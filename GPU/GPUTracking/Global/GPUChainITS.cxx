// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUChainITS.cxx
/// \author David Rohr

#include "GPUChainITS.h"
#include "GPUReconstructionIncludesITS.h"
#include "DataFormatsITS/TrackITS.h"
#include <algorithm>

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::its;

GPUChainITS::~GPUChainITS()
{
  mITSTrackerTraits.reset();
  mITSVertexerTraits.reset();
}

GPUChainITS::GPUChainITS(GPUReconstruction* rec, unsigned int maxTracks) : GPUChain(rec), mMaxTracks(maxTracks) {}

void GPUChainITS::RegisterPermanentMemoryAndProcessors() { mRec->RegisterGPUProcessor(&processors()->itsFitter, GetRecoStepsGPU() & RecoStep::ITSTracking); }

void GPUChainITS::RegisterGPUProcessors()
{
  if (GetRecoStepsGPU() & RecoStep::ITSTracking) {
    mRec->RegisterGPUDeviceProcessor(&processorsShadow()->itsFitter, &processors()->itsFitter);
  }
}

void GPUChainITS::MemorySize(size_t& gpuMem, size_t& pageLockedHostMem)
{
  gpuMem = mMaxTracks * sizeof(GPUITSTrack) + GPUCA_MEMALIGN;
  pageLockedHostMem = gpuMem;
}

int GPUChainITS::Init() { return 0; }

TrackerTraits* GPUChainITS::GetITSTrackerTraits()
{
  if (mITSTrackerTraits == nullptr) {
    mRec->GetITSTraits(&mITSTrackerTraits, nullptr, nullptr);
    mITSTrackerTraits->SetRecoChain(this);
  }
  return mITSTrackerTraits.get();
}

VertexerTraits* GPUChainITS::GetITSVertexerTraits()
{
  if (mITSVertexerTraits == nullptr) {
    mRec->GetITSTraits(nullptr, &mITSVertexerTraits, nullptr);
  }
  return mITSVertexerTraits.get();
}

TimeFrame* GPUChainITS::GetITSTimeframe()
{
  if (mITSTimeFrame == nullptr) {
    mRec->GetITSTraits(nullptr, nullptr, &mITSTimeFrame);
  }
  LOGP(info, "Setting timeFrame allocator to external");
  mITSTimeFrame->setExtAllocator(true); // meaningful only when GPU
  mITSTimeFrame->setChain(this);
  return mITSTimeFrame.get();
}

int GPUChainITS::PrepareEvent() { return 0; }

int GPUChainITS::Finalize() { return 0; }

int GPUChainITS::RunChain() { return 0; }