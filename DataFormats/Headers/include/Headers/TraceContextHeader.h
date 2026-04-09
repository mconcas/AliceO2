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

/// \file TraceContextHeader.h
/// \brief O2 header carrying a W3C Trace Context (traceparent) for distributed tracing.
///
/// Injected into the O2 header stack of outgoing DPL messages by the TracingService.
/// Extracted on the receiving side to create child spans, enabling a connected
/// topology graph across DPL devices.

#ifndef O2_HEADERS_TRACECONTEXTHEADER_H
#define O2_HEADERS_TRACECONTEXTHEADER_H

#include "Headers/DataHeader.h"
#include <cstring>

namespace o2::header
{

/// \struct TraceContextHeader
/// \brief Carries one W3C traceparent string (55 bytes) in the O2 header stack.
///
/// Wire format of the traceparent field:
///   "00-<32 hex traceId>-<16 hex spanId>-<02 hex flags>\0"
///
/// The struct is fixed-size and trivially copyable so it is safe to
/// embed directly in a FairMQ message header region.
struct TraceContextHeader : public BaseHeader {
  static const uint32_t sVersion;
  static const o2::header::HeaderType sHeaderType;
  static const o2::header::SerializationMethod sSerializationMethod;

  /// W3C traceparent value, null-terminated.
  /// "00-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-XXXXXXXXXXXXXXXX-XX"
  ///   2  +1+ 32                          +1+ 16           +1+2  = 55 chars + '\0'
  static constexpr size_t kTraceparentMaxLen = 56; // 55 chars + NUL

  char traceparent[kTraceparentMaxLen];

  TraceContextHeader()
    : BaseHeader(sizeof(TraceContextHeader), sHeaderType, sSerializationMethod, sVersion)
  {
    std::memset(traceparent, 0, sizeof(traceparent));
  }

  explicit TraceContextHeader(const char* tp)
    : BaseHeader(sizeof(TraceContextHeader), sHeaderType, sSerializationMethod, sVersion)
  {
    std::memset(traceparent, 0, sizeof(traceparent));
    if (tp) {
      std::strncpy(traceparent, tp, kTraceparentMaxLen - 1);
    }
  }

  /// Returns true if this header holds a non-empty traceparent.
  bool valid() const { return traceparent[0] != '\0'; }
};

} // namespace o2::header

#endif // O2_HEADERS_TRACECONTEXTHEADER_H
