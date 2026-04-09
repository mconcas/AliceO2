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

#include "Headers/TraceContextHeader.h"

// Out-of-line definitions for BaseHeader static members.
const uint32_t o2::header::TraceContextHeader::sVersion = 1;
const o2::header::HeaderType o2::header::TraceContextHeader::sHeaderType = "TrcCtx  ";
const o2::header::SerializationMethod o2::header::TraceContextHeader::sSerializationMethod = o2::header::gSerializationMethodNone;
