# DPL Distributed Tracing — Technical Specification

**Branch:** `otel-tracing`  
**Commits covered:** `3c5fed6f51` → `f27865f1f0`  
**Date:** 2026-04-21

---

## 1. Overview

This feature integrates OpenTelemetry-compatible distributed tracing into the
Data Processing Layer (DPL). Every processing invocation performed by a DPL
device emits an OpenTelemetry span. Spans are propagated between devices via
the O2 header stack, forming a connected trace that mirrors the DPL topology.
Fault conditions — timeslice slots that can never be completed because upstream
inputs are missing — are additionally represented as error spans so that
RED-metric (Rate / Errors / Duration) dashboards capture them.

The implementation is entirely opt-in at build time: when the optional
`AliceO2::Tracing` CMake target is absent, every tracing path compiles away to
nothing and there is no runtime overhead.

---

## 2. DPL Background

The following DPL constructs are referenced throughout this document.

### 2.1 Device

A DPL *device* (`DataProcessorSpec`) is an independent OS process that
performs one stage of the data-processing pipeline. Devices communicate
exclusively through typed, framed FairMQ messages. The DPL driver spawns,
monitors, and connects devices according to a topology described at
workflow-configuration time.

### 2.2 ServiceSpec and ServiceRegistry

Services are singleton objects owned by a device and accessible to algorithms
and other services via the `ServiceRegistry`. A `ServiceSpec` is a descriptor
that defines:

| Field | Purpose |
|---|---|
| `name` | Unique identifier used for lookup |
| `init` | Factory called once when the device starts; returns a `ServiceHandle` |
| `configure` | Called after option parsing |
| `preProcessing` | Called before each user algorithm invocation |
| `postProcessing` | Called after each user algorithm invocation |
| `exit` | Called on device shutdown; responsible for resource cleanup |
| `kind` | Concurrency model; `ServiceKind::Serial` means one instance per device |

`CommonServices::defaultServices()` assembles the list of services that are
unconditionally present in every device.

### 2.3 ProcessingContext

The `ProcessingContext` passed to `preProcessing` / `postProcessing` callbacks
(and to user algorithms) provides access to:

- `inputs()` — the `InputRecord` holding all FairMQ messages for the current
  timeslice.
- `outputs()` — the `DataAllocator` for constructing outgoing messages.
- `services()` — the `ServiceRegistryRef` for retrieving other services.

### 2.4 O2 Header Stack

Every DPL message carries a *header stack*: a contiguous, variable-length
sequence of fixed-size `BaseHeader`-derived structs prepended to the payload.
The stack is built with `o2::header::Stack` and transmitted in the FairMQ
header region. Standard entries are `DataHeader` (data origin / description /
subspecification) and `DataProcessingHeader` (timeslice, start time). Any
number of additional headers can be appended; `o2::header::get<T*>()` performs
a typed walk of the stack to locate a header by its 8-byte `sHeaderType`
descriptor.

### 2.5 DataRelayer

The `DataRelayer` is the message routing core inside a device. It maintains a
two-dimensional cache indexed by `(TimesliceSlot, InputIndex)`. Incoming
messages are placed into slots keyed by their timeslice value. When all
expected inputs for a slot have arrived, the slot is marked complete and
dispatched to the user algorithm. `setOldestPossibleInput()` is called whenever
downstream devices advance the oldest-possible timeslice watermark; it scans
for slots that can now never be completed and prunes them, logging warnings
about the missing inputs.

### 2.6 DataAllocator

`DataAllocator` is the output API available to user algorithms. Its internal
`headerMessageFromOutput()` method assembles the `o2::header::Stack` for each
outgoing message and allocates the FairMQ header buffer from the
transport-specific memory pool.

### 2.7 Signpost

`O2_SIGNPOST_START(log, id, name, ...)` and `O2_SIGNPOST_END(log, id, name, ...)`
are lightweight instrumentation macros defined in `Framework/Foundation/include/Framework/Signpost.h`.
On macOS they call into the system `os_signpost` API; on other platforms they
fall through to a minimal in-process interval tracker used for stack-trace
capture. They are already sprinkled throughout the DPL internals (e.g. around
`DataRelayer` slot operations and the `DeviceState` event loop).

### 2.8 TimingInfo

`TimingInfo` is a per-stream service that exposes the current timeslice value
and run number to algorithm and service callbacks. It is populated from the
incoming `DataProcessingHeader` immediately before `preProcessing` runs.

---

## 3. New Constructs

### 3.1 `TraceContextHeader` (`DataFormats/Headers`)

```
DataFormats/Headers/include/Headers/TraceContextHeader.h
DataFormats/Headers/src/TraceContextHeader.cxx
```

A new `BaseHeader`-derived struct that carries one W3C `traceparent` string
across a device boundary. The wire format of `traceparent` is:

```
00-<32 hex traceId>-<16 hex spanId>-<02 hex flags>
```

This is exactly 55 characters plus a NUL terminator, so the field is allocated
as a 56-byte fixed-size array, making the header trivially copyable and safe to
embed in any FairMQ memory region.

**Static descriptor fields:**

| Field | Value |
|---|---|
| `sHeaderType` | `"TrcCtx  "` (8-byte descriptor) |
| `sSerializationMethod` | `gSerializationMethodNone` |
| `sVersion` | `1` |

The `valid()` predicate tests whether `traceparent[0] != '\0'`, allowing
downstream code to distinguish an injected header from a default-constructed
one without inspecting the W3C string.

The header is added to the `O2::Headers` CMake library in
`DataFormats/Headers/CMakeLists.txt`.

---

### 3.2 `DPLTracingService` (`Framework/Core`)

```
Framework/Core/src/DPLTracingService.h
```

A header-only service class (`struct DPLTracingService`) implementing the full
tracing lifecycle for one device process. All methods are guarded by
`#ifdef O2_WITH_DPL_TRACING`; without the Tracing library the struct contains
only a `processingCount` counter and every method is a no-op.

#### 3.2.1 `init(backendUrl, deviceName, runNumber)`

Creates an `o2::tracing::Tracer` via `TracingFactory::Get(backendUrl)` and
attaches two global attributes — `service.name` (set to the DPL device name)
and optionally `run.number`. It then registers the `SignpostSpanBridge` (see
§3.3) and installs two `std::atomic` function pointers,
`o2_signpost_start_hook` and `o2_signpost_end_hook`, that are read by the
`O2_SIGNPOST_START` / `O2_SIGNPOST_END` macros on every invocation.

The `--tracing-backend` CLI option selects the backend. The built-in values
understood by `TracingFactory` are:

| URI | Behaviour |
|---|---|
| `no-op://` | Default; tracer discards all spans |
| `stdout://` | Prints spans to standard output (development) |
| `otlp-grpc://host:port` | Exports via OTLP gRPC to a collector |

#### 3.2.2 `beginSpan(ProcessingContext&)` — called from `preProcessing`

1. Increments `processingCount`.
2. Iterates over all inputs in `InputRecord`. For each input that carries a
   `TraceContextHeader`, parses the W3C `traceparent` into a
   `SpanContext`. If more than one valid context is found, the first becomes a
   *span link* rather than the parent; only the last valid context is used as
   the direct parent. This models fan-in: a single OTEL parent edge is drawn
   while all contributing upstream spans remain visible via links.
3. Creates a span named `"dpl/process"`. The span kind is:
   - `SpanKind::Server` when an upstream `traceparent` was found (the device is
     receiving work from an upstream service).
   - `SpanKind::Client` when no upstream context exists (the device is the
     origin of a new trace, e.g. a raw-data reader).
4. Attaches the timeslice and run number from `TimingInfo` as span attributes.
5. Writes the resulting `SpanContext` to the thread-local
   `tDPLCurrentSpanCtx`, making it available to `SignpostSpanBridge` hooks
   fired from the same thread.

#### 3.2.3 `currentOutgoingContext()` — called from `DataAllocator`

Returns a `TraceContextHeader` to be injected into an outgoing message. For
each call a short-lived `CLIENT`-kind span named `"dpl/send"` is created as a
child of the current processing span, its context is serialised to W3C format,
and the span is immediately closed. Creating a separate CLIENT span per output
call (rather than reusing the SERVER span) is required by OTEL service-map
processors: they identify a service-to-service edge by finding a CLIENT span
whose `traceparent` context is the parent of a downstream SERVER span. For a
fan-out device sending to *N* downstream devices, *N* distinct CLIENT spans
are created, allowing independent correlation of each edge.

#### 3.2.4 `emitDroppedSlot(timeslice, partialTch, missingInputs)`

Called by `DataRelayer::setOldestPossibleInput()` when a timeslice slot is
pruned because it will never be completed. Emits a standalone `"dpl/process"`
span with `SpanStatus::Error` and the attribute `dpl.missing_inputs` (a
comma-separated list of `DataSpecUtils::describe()` strings for each absent
input). If any partial input arrived before the slot was dropped, its
`TraceContextHeader` is used to parent the error span into the upstream trace
tree; otherwise the span is a root with no parent.

#### 3.2.5 `endSpan()` — called from `postProcessing`

Closes the current span, clears `tDPLCurrentSpanCtx` to prevent stray
signpost hooks from parenting sub-spans to an already-closed span, and returns
the span's context serialised as a `TraceContextHeader`.

---

### 3.3 `SignpostSpanBridge`

A process-global singleton (in `DPLTracingService.h`) that maps signpost
interval identifiers to active `o2::tracing::Span` objects. It holds a pointer
to the device `Tracer` (set by `DPLTracingService::init`) and an
`std::unordered_map<int64_t, std::unique_ptr<Span>>` protected by a mutex.

When `startSpan(name, id)` is called it checks that `tDPLCurrentSpanCtx` is
valid — i.e. that `beginSpan` has run and the device-level span is open. If
not, the call is silently dropped to avoid orphaned root spans from
`O2_SIGNPOST_START` invocations that fire outside of processing (e.g. during
device initialisation). When `endSpan(name, id)` is called the corresponding
span is moved out of the map and closed. This creates a two-level span
hierarchy: the device-level `"dpl/process"` span as root, signpost intervals
as immediate children.

---

### 3.4 `Signpost.h` — hook extension

Two `std::atomic<o2_signpost_hook_fn>` globals are added:

```cpp
extern std::atomic<o2_signpost_hook_fn> o2_signpost_start_hook;
extern std::atomic<o2_signpost_hook_fn> o2_signpost_end_hook;
```

where `o2_signpost_hook_fn` is `void (*)(const char* name, int64_t id)`.

The `O2_SIGNPOST_START` and `O2_SIGNPOST_END` macros are extended with a
trailing `do { ... } while(0)` block that performs a relaxed atomic load and,
only if the pointer is non-null (`O2_BUILTIN_UNLIKELY`), calls through to the
hook. The `unlikely` branch hint keeps the unconditional overhead to one atomic
load and a not-taken branch per signpost site when no tracing backend is
configured — cost on the order of a cache-line read.

Definitions of the two atomics live in `Signpost.h` behind
`#ifdef O2_SIGNPOST_IMPLEMENTATION`, which is only triggered in one translation
unit.

---

## 4. Integration Points

### 4.1 `CommonServices::tracingSpec()`

`Framework/Core/src/CommonServices.cxx`

The pre-existing stub `TracingInfrastructure` (a struct with a single counter
and no-op callbacks) is replaced with a fully functional `ServiceSpec` backed
by `DPLTracingService`. The updated spec:

- **`init`**: constructs a `DPLTracingService`, reads the `--tracing-backend`
  option value, and calls `svc->init()` unless the value is `"no-op://"`.
  Skipping `init` when the backend is `no-op://` avoids loading any tracing
  library code at all.
- **`preProcessing`**: calls `svc->beginSpan(ctx)`.
- **`postProcessing`**: calls `svc->endSpan()`. The returned
  `TraceContextHeader` is currently discarded (`[[maybe_unused]]`); injecting
  it into the output is handled by `DataAllocator` (§4.2).
- **`exit`**: nullifies both signpost hook atomics before deleting the service,
  preventing any late-firing signpost from accessing freed memory.

`tracingSpec()` is appended to the `defaultServices()` vector, making it
present in every DPL device.

A thread-local `tDPLCurrentSpanCtx` of type `o2::tracing::SpanContext` is
defined in `CommonServices.cxx` (behind `#ifdef O2_WITH_DPL_TRACING`) and
declared `extern` in `DPLTracingService.h` so that both `Framework::Core` and
`Framework::DataTakingSupport` can link to the same instance.

### 4.2 `DataAllocator::headerMessageFromOutput()`

`Framework/Core/src/DataAllocator.cxx`

Before constructing the `o2::header::Stack` for an outgoing message, the
method now queries the tracing service:

```cpp
auto tch = mRegistry.get<DPLTracingService>().currentOutgoingContext();
if (tch.valid()) {
    return o2::pmr::getMessage(o2::header::Stack{channelAlloc, dh, dph, spec.metaHeader, tch});
}
return o2::pmr::getMessage(o2::header::Stack{channelAlloc, dh, dph, spec.metaHeader});
```

When a span is active and the tracing backend is configured, the
`TraceContextHeader` is appended to the stack. Every call to `make<T>()`,
`snapshot()`, or `adopt()` that goes through `headerMessageFromOutput()` will
carry context, covering the common output paths without requiring any change to
user algorithm code.

### 4.3 `DataRelayer::setOldestPossibleInput()`

`Framework/Core/src/DataRelayer.cxx`

Two additions are made inside the loop that prunes incomplete slots:

1. **Context extraction**: when an input element is found present in the cache,
   the code attempts to read a `TraceContextHeader` from its header message.
   The first valid header found is stored as `partialTch`.

2. **Error span emission**: after the loop, if `didDrop` is true, the list of
   absent inputs is accumulated as a comma-separated `missingInputs` string and
   `mContext.get<DPLTracingService>().emitDroppedSlot(...)` is called.

This produces an error span for every dropped slot, regardless of whether
partial data arrived.

### 4.4 `DeviceSpecHelpers` and `runDataProcessing`

`--tracing-backend` is registered as a forwarded device option in
`DeviceSpecHelpers::getForwardedDeviceOptions()` and given the default value
`"no-op://"` in `runDataProcessing.cxx`'s `doChild()`. This means:

- The driver propagates the option to all child devices without each
  `DataProcessorSpec` needing to declare it.
- The default value disables tracing entirely, preserving the zero-cost
  behaviour for deployments that do not configure a backend.

### 4.5 InfoLogger trace correlation (`DataTakingSupport/Plugin.cxx`)

`Framework/DataTakingSupport/src/Plugin.cxx`

The InfoLogger sink (created by `createInfoLoggerSinkHelper`) is extended to
inject the current W3C trace identifiers into every log line emitted during
processing:

```cpp
msgCtx.setField(InfoLoggerContext::FieldName::TraceId, w3c.substr(3, 32));
msgCtx.setField(InfoLoggerContext::FieldName::SpanId,  w3c.substr(36, 16));
```

The substrings are derived from the W3C `traceparent` format
(`00-<traceId32>-<spanId16>-<flags2>`): `traceId` begins at offset 3 and is 32
characters; `spanId` begins at offset 36 and is 16 characters. Correlation is
only applied when `tDPLCurrentSpanCtx.valid()` is true, i.e. when a
`beginSpan` call has set up an active span on the current thread. This allows
log messages to be joined with trace spans in backends that support unified
telemetry (e.g. OpenSearch with the OTel plugin, or Grafana Tempo + Loki).

---

## 5. Build System

### 5.1 Optional dependency

`dependencies/O2Dependencies.cmake` adds:

```cmake
find_package(Tracing CONFIG)
set_package_properties(Tracing PROPERTIES TYPE OPTIONAL ...)
```

The package is `OPTIONAL`; a build without it succeeds and produces a fully
functional binary with all tracing paths compiled to no-ops.

### 5.2 Compile definition

`Framework/Core/CMakeLists.txt` links `AliceO2::Tracing` as a generator
expression and propagates `O2_WITH_DPL_TRACING` as a `PRIVATE` compile
definition to `Framework`. `Framework/DataTakingSupport/CMakeLists.txt`
applies the same pattern for `FrameworkDataTakingSupport`.

Using `PRIVATE` ensures the define is not visible to downstream consumers of
`O2::Framework`, avoiding accidental compilation of tracing code in targets
that do not link the Tracing library.

---

## 6. Span Hierarchy and Propagation Model

```
┌──────────────────────────────────────────────────────────┐
│  Device A (SpanKind::Client — no upstream context)        │
│                                                          │
│  dpl/process  ─────────────────────────────────────────  │
│    └── O2_SIGNPOST_START("foo")  (child via bridge)      │
│    └── O2_SIGNPOST_START("bar")  (child via bridge)      │
│                                                          │
│  DataAllocator::headerMessageFromOutput()                │
│    └── dpl/send  [CLIENT, closed immediately]            │
│         traceparent injected into outgoing message       │
└──────────────────────────────────────────────────────────┘
                           │  FairMQ message + TraceContextHeader
                           ▼
┌──────────────────────────────────────────────────────────┐
│  Device B (SpanKind::Server — upstream context present)   │
│                                                          │
│  dpl/process  ← parent = Device A's dpl/send span        │
│    └── sub-spans from O2_SIGNPOST_START calls            │
└──────────────────────────────────────────────────────────┘
```

- **CLIENT/SERVER pairing**: OTEL service maps require a CLIENT span in device
  A whose context is the parent of a SERVER span in device B. The `dpl/send`
  span plays the CLIENT role; the `dpl/process` span in the receiving device
  plays the SERVER role.

- **Fan-in**: if device C receives messages from both A and B, `beginSpan`
  picks the last valid context as the single OTEL parent and records all other
  contexts as span links. The service map draws one primary edge while the
  trace view shows all contributing upstream spans.

- **Fan-out**: each call to `headerMessageFromOutput()` independently creates
  and closes a `dpl/send` span. The N downstream devices each receive a
  distinct `traceparent`, resulting in N independent edges in the service map.

- **Error spans**: dropped timeslice slots emit a `dpl/process` SERVER span
  with `SpanStatus::Error`, immediately closed with no sub-spans. APM RED
  metrics will capture these as errors against the affected device's service
  name.

---

## 7. Known Limitations and Open Items

- **`endSpan` return value unused in `postProcessing`**: the `TraceContextHeader`
  returned by `endSpan()` is marked `[[maybe_unused]]` in the `postProcessing`
  lambda. Context injection is handled by `currentOutgoingContext()` inside
  `DataAllocator`, so this is not a functional gap, but the discarded return
  value suggests the postProcessing path may have originally been intended to
  carry context for a different injection mechanism.

- **`traceparent` links not yet working**: commit `f27865f1f0` ("Add links (not
  working?)") explicitly notes that span links — the mechanism used for fan-in
  devices to record non-primary upstream contexts — are not yet producing
  correct output in the trace backend. The data structures are populated
  correctly in `beginSpan`; the issue is likely in how `TracingFactory` or the
  underlying OTLP exporter serialises link records.

- **Run number unavailable at init**: `DPLTracingService::init` receives
  `static_cast<uint32_t>(-1)` as the run number because `TimingInfo` is
  per-stream and not valid during service initialisation. The run number is
  attached as a span attribute in `beginSpan` instead, once `TimingInfo` is
  populated for the current timeslice. The global `addGlobalAttribute` call for
  `kRunNumber` in `init` is therefore skipped.

- **No backpressure on `SignpostSpanBridge`**: the unordered map has no bound.
  If an `O2_SIGNPOST_START` call is made but the corresponding
  `O2_SIGNPOST_END` is never executed (e.g. due to an exception), the entry
  will remain in the map for the lifetime of the process. This is unlikely in
  practice since signpost intervals are always paired, but there is no
  defensive cleanup.
