// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DPLTracingService.h
/// \brief DPL service that wraps the Tracing library and instruments each
///        processing invocation with an OpenTelemetry span.
///
/// When O2_WITH_DPL_TRACING is defined (i.e. AliceO2::Tracing was found at
/// CMake time) the service:
///   - Creates a Tracer from the --tracing-backend CLI option on init.
///   - Extracts a W3C TraceContextHeader from the first input message.
///   - Opens a child span in preProcessing and ends it in postProcessing.
///   - Injects the resulting SpanContext back into outgoing messages via
///     the DataAllocator header stack (handled by the DPL output machinery).
///
/// Without the library the service is a zero-cost stub.

#ifndef O2_FRAMEWORK_DPLTRACINGSERVICE_H_
#define O2_FRAMEWORK_DPLTRACINGSERVICE_H_

#ifdef O2_WITH_DPL_TRACING
#include "Tracing/Tracer.h"
#include "Tracing/Span.h"
#include "Tracing/SpanContext.h"
#include "Tracing/TracingFactory.h"
#include "Tracing/Tags.h"
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>
#endif

#include "Framework/Signpost.h"

#include "Headers/TraceContextHeader.h"
#include "Framework/DataRef.h"
#include "Framework/DataRefUtils.h"
#include "Framework/InputRecord.h"
#include "Framework/ProcessingContext.h"
#include "Framework/ServiceRegistryRef.h"
#include "Framework/TimingInfo.h"
#include "Framework/DeviceSpec.h"

namespace o2::framework
{

#ifdef O2_WITH_DPL_TRACING
/// Thread-local context of the current device-level span, updated by
/// beginSpan / endSpan so that signpost hooks can create proper child spans.
inline thread_local o2::tracing::SpanContext tDPLCurrentSpanCtx{};

/// Process-global map from signpost id → active sub-span.
/// Protected by a mutex; reads on the hot-path are fast (span create/end
/// is rare compared to actual data processing).
struct SignpostSpanBridge {
  static SignpostSpanBridge& instance()
  {
    static SignpostSpanBridge s;
    return s;
  }

  void setTracer(o2::tracing::Tracer* t) { mTracer = t; }

  void startSpan(const char* name, int64_t id)
  {
    if (!mTracer || !tDPLCurrentSpanCtx.valid()) {
      // No active device span yet (e.g. run_callback fires before preProcessing).
      // Skip to avoid orphaned root spans that duplicate the Phase 3 device span.
      return;
    }
    auto span = mTracer->startSpan(name, tDPLCurrentSpanCtx);
    std::lock_guard<std::mutex> lock(mMutex);
    mSpans.emplace(id, std::move(span));
  }

  void endSpan(const char* /*name*/, int64_t id)
  {
    std::unique_ptr<o2::tracing::Span> span;
    {
      std::lock_guard<std::mutex> lock(mMutex);
      auto it = mSpans.find(id);
      if (it == mSpans.end()) {
        return;
      }
      span = std::move(it->second);
      mSpans.erase(it);
    }
    span->end();
  }

 private:
  o2::tracing::Tracer* mTracer{nullptr};
  std::mutex mMutex;
  std::unordered_map<int64_t, std::unique_ptr<o2::tracing::Span>> mSpans;
};
#endif // O2_WITH_DPL_TRACING

struct DPLTracingService {
#ifdef O2_WITH_DPL_TRACING
  std::unique_ptr<o2::tracing::Tracer> tracer;
  std::unique_ptr<o2::tracing::Span> currentSpan;
#endif
  int processingCount{0};

  /// Initialise from the --tracing-backend option value.
  /// No-op if tracing is not compiled in.
  void init(const std::string& backendUrl,
            const std::string& deviceName,
            uint32_t runNumber)
  {
#ifdef O2_WITH_DPL_TRACING
    tracer = o2::tracing::TracingFactory::Get(backendUrl);
    tracer->addGlobalAttribute(o2::tracing::tags::kServiceName, deviceName);
    if (runNumber != static_cast<uint32_t>(-1)) {
      tracer->addGlobalAttribute(o2::tracing::tags::kRunNumber,
                                 std::to_string(runNumber));
    }
    // Register the signpost bridge and install global hooks so that
    // O2_SIGNPOST_START / O2_SIGNPOST_END calls emit child spans.
    SignpostSpanBridge::instance().setTracer(tracer.get());
    o2_signpost_start_hook.store(
      [](const char* name, int64_t id) {
        SignpostSpanBridge::instance().startSpan(name, id);
      },
      std::memory_order_relaxed);
    o2_signpost_end_hook.store(
      [](const char* name, int64_t id) {
        SignpostSpanBridge::instance().endSpan(name, id);
      },
      std::memory_order_relaxed);
#endif
  }

  /// Called in preProcessing: extract parent context from headers, start span.
  void beginSpan(ProcessingContext& ctx)
  {
    processingCount++;
#ifdef O2_WITH_DPL_TRACING
    if (!tracer) {
      return;
    }
    // Scan all inputs for TraceContextHeaders.
    // For fan-in devices (multiple upstream services) OTEL only supports one
    // parent span, so we pick the LAST valid context found — this tends to be
    // the deepest hop in a chain (e.g. C rather than A for a diamond A→C→D).
    // All other upstream contexts are stored as span links so the full fan-in
    // is visible in the trace view even though the service map shows one edge.
    o2::tracing::SpanContext parentCtx{};
    std::vector<o2::tracing::SpanContext> linkCtxs;
    auto& inputs = ctx.inputs();
    for (int i = 0; i < inputs.size(); ++i) {
      auto ref = inputs.getByPos(i);
      if (ref.header == nullptr) {
        continue;
      }
      auto* tch = o2::header::get<o2::header::TraceContextHeader*>(ref.header);
      if (tch && tch->valid()) {
        auto newCtx = o2::tracing::SpanContext::fromW3C(tch->traceparent);
        if (parentCtx.valid()) {
          linkCtxs.push_back(parentCtx); // demote previous to link
        }
        parentCtx = newCtx;
      }
    }

    auto& timing = ctx.services().get<TimingInfo>();
    // SERVER = receiving upstream data; CLIENT = originating (no upstream context).
    auto kind = parentCtx.valid() ? o2::tracing::SpanKind::Server
                                  : o2::tracing::SpanKind::Client;
    currentSpan = tracer->startSpan("dpl/process", parentCtx, kind);
    for (auto& lctx : linkCtxs) {
      currentSpan->addLink(lctx);
    }
    currentSpan->setAttribute(o2::tracing::tags::kTimeslice,
                              static_cast<int64_t>(timing.timeslice));
    if (timing.runNumber != static_cast<uint32_t>(-1)) {
      currentSpan->setAttribute(o2::tracing::tags::kRunNumber,
                                static_cast<int64_t>(timing.runNumber));
    }
    // Expose this span's context thread-locally so O2_SIGNPOST_START hooks
    // can create proper child spans.
    tDPLCurrentSpanCtx = currentSpan->context();
#endif
  }

  /// Returns the current span context as a TraceContextHeader for injection into
  /// outgoing message headers. Call this while a span is active (between beginSpan
  /// and endSpan). Returns an empty/invalid header when no span is active.
  ///
  /// A fresh CLIENT child span (dpl/send) is created and immediately closed for
  /// every call. Service-map processors require a CLIENT→SERVER pair across a
  /// service boundary to draw a topology edge. Creating one span per output call
  /// (rather than a shared singleton) ensures that fan-out devices sending to
  /// multiple downstream services each produce a distinct CLIENT span, allowing
  /// the service-map processor to independently correlate each edge.
  o2::header::TraceContextHeader currentOutgoingContext() const
  {
#ifdef O2_WITH_DPL_TRACING
    if (currentSpan && tracer) {
      auto sendSpan = tracer->startSpan("dpl/send", currentSpan->context(),
                                        o2::tracing::SpanKind::Client);
      auto ctx = sendSpan->context();
      sendSpan->end();
      if (ctx.valid()) {
        auto w3c = ctx.toW3C();
        return o2::header::TraceContextHeader(w3c.c_str());
      }
    }
#endif
    return o2::header::TraceContextHeader{};
  }

  /// Called by DataRelayer when a timeslice slot is dropped because it can never
  /// be completed (missing upstream inputs). Emits a SERVER span with
  /// SpanStatus::Error so that APM error-rate (RED) metrics capture the failure.
  /// partialTch may be default-constructed (invalid) when none of the expected
  /// inputs arrived; in that case the span has no upstream parent.
  void emitDroppedSlot(uint64_t timeslice,
                       const o2::header::TraceContextHeader& partialTch,
                       std::string_view missingInputs)
  {
#ifdef O2_WITH_DPL_TRACING
    if (!tracer) {
      return;
    }
    o2::tracing::SpanContext parentCtx{};
    if (partialTch.valid()) {
      parentCtx = o2::tracing::SpanContext::fromW3C(partialTch.traceparent);
    }
    auto span = tracer->startSpan("dpl/process", parentCtx, o2::tracing::SpanKind::Server);
    span->setAttribute(o2::tracing::tags::kTimeslice, static_cast<int64_t>(timeslice));
    if (!missingInputs.empty()) {
      span->setAttribute("dpl.missing_inputs", std::string(missingInputs));
    }
    span->setStatus(o2::tracing::SpanStatus::Error, "incomplete slot dropped");
    span->end();
#endif
  }

  /// Called in postProcessing: end the span and return the propagatable context.
  /// Returns an empty/invalid context when tracing is not compiled in.
  o2::header::TraceContextHeader endSpan()
  {
    o2::header::TraceContextHeader tch{};
#ifdef O2_WITH_DPL_TRACING
    if (currentSpan) {
      auto ctx = currentSpan->context();
      currentSpan->end();
      currentSpan.reset();
      tDPLCurrentSpanCtx = {}; // clear thread-local so stray signpost hooks don't parent to a closed span
      if (ctx.valid()) {
        auto w3c = ctx.toW3C();
        tch = o2::header::TraceContextHeader(w3c.c_str());
      }
    }
#endif
    return tch;
  }
};

} // namespace o2::framework

#endif // O2_FRAMEWORK_DPLTRACINGSERVICE_H_
