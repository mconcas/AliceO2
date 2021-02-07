// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DetectorsVertexing/SVertexer.h"
#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsITS/TrackITS.h"

// move eventually
#include "ITStrackingCUDA/Vector.h"
#include "ITStrackingCUDA/Array.h"

//Macro for checking GPU errors following a cuda launch or api call
#define gpuCheckError()                                                               \
  {                                                                                   \
    cudaError_t e = cudaGetLastError();                                               \
    if (e != cudaSuccess) {                                                           \
      printf("GPU failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
      exit(0);                                                                        \
    }                                                                                 \
  }

namespace o2
{
namespace vertexing
{

struct SVParamsGPU {
  // Struct to pass configuration for SVertexer as static symbol on GPU
  unsigned char useAbsDCA = false;
  unsigned char propagateToPCA = false;
  float maxRIni = 0.f;
  float minParamChange = 0.f;
  float minRelChi2Change = 0.f;
  float maxDZIni = 0.f;
  float maxChi2 = 0.f;
  float bz = 0.f;
};

__constant__ SVParamsGPU gpuStaticConf;

namespace kernels
{
GPUg() void testFitterKernel()
{
  o2::vertexing::DCAFitterN<2> mFitter2Prong;

  printf("End.\n");
}
} // namespace kernels

using PID = o2::track::PID;
using TrackTPCITS = o2::dataformats::TrackTPCITS;
using TrackITS = o2::its::TrackITS;
using TrackTPC = o2::tpc::TrackTPC;

void SVertexerCUDA::init()
{
  mSVParams = &SVertexerParams::Instance();
  auto bz = o2::base::Propagator::Instance()->getNominalBz();

  SVParamsGPU parsHost;

  parsHost.useAbsDCA = mSVParams->useAbsDCA;
  parsHost.maxRIni = mSVParams->maxRIni;
  parsHost.minParamChange = mSVParams->minParamChange;
  parsHost.minRelChi2Change = mSVParams->minRelChi2Change;
  parsHost.maxDZIni = mSVParams->maxDZIni;
  parsHost.maxChi2 = mSVParams->maxChi2;
  parsHost.propagateToPCA = false;
  parsHost.bz = bz;

  // precalculated selection cuts
  mMinR2ToMeanVertex = mSVParams->minRfromMeanVertex * mSVParams->minRfromMeanVertex;
  mMaxDCAXY2ToMeanVertex = mSVParams->maxDCAXYfromMeanVertex * mSVParams->maxDCAXYfromMeanVertex;
  mMinCosPointingAngle = mSVParams->minCosPointingAngle;

  mV0Hyps[SVertexerParams::Photon].set(PID::Photon, PID::Electron, PID::Electron, mSVParams->pidCutsPhoton, bz);
  mV0Hyps[SVertexerParams::K0].set(PID::K0, PID::Pion, PID::Pion, mSVParams->pidCutsK0, bz);
  mV0Hyps[SVertexerParams::Lambda].set(PID::Lambda, PID::Proton, PID::Pion, mSVParams->pidCutsLambda, bz);
  mV0Hyps[SVertexerParams::AntiLambda].set(PID::Lambda, PID::Pion, PID::Proton, mSVParams->pidCutsLambda, bz);
  mV0Hyps[SVertexerParams::HyperTriton].set(PID::HyperTriton, PID::Helium3, PID::Pion, mSVParams->pidCutsHTriton, bz);
  mV0Hyps[SVertexerParams::AntiHyperTriton].set(PID::HyperTriton, PID::Pion, PID::Helium3, mSVParams->pidCutsHTriton, bz);

  // Load configuration to constant static symbol on GPU
  cudaMemcpyToSymbol(gpuStaticConf, &parsHost, sizeof(parsHost));
  printf("Calling kernel\n");
  kernels::testFitterKernel<<<1, 1>>>();
  cudaDeviceSynchronize();
  gpuCheckError();

  // NB: no DCA-fitter has been configured yet.
}

void SVertexerCUDA::process(const gsl::span<const PVertex>& vertices,   // primary vertices
                            const gsl::span<const GIndex>& trackIndex,  // Global ID's for associated tracks
                            const gsl::span<const VRef>& vtxRefs,       // references from vertex to these track IDs
                            const o2d::GlobalTrackAccessor& tracksPool, // accessor to various tracks
                            std::vector<V0>& v0s,                       // found V0s
                            std::vector<RRef>& vtx2V0refs               // references from PVertex to V0
)
{
  std::unordered_map<uint64_t, int> cache; // cache for tested combinations, the value >0 will give the entry of prevalidated V0 in the v0sTmp
  std::vector<V0> v0sTmp(1);               // 1st one is dummy!
  std::vector<int> v0sIdx;                 // id's in v0sTmp used attached to p.vertices
  std::vector<RRef> pv2v0sRefs;            // p.vertex to v0 index references
  std::vector<char> selQ(trackIndex.size(), 0);

  kernels::testFitterKernel<<<1, 1>>>();
  gpuCheckError();
}

void SVertexerCUDA::testDCAFitterGPU(std::vector<o2::track::TrackParCov>& trks)
{
  o2::its::gpu::Vector<o2::track::TrackParCov> tr;
}

} // namespace vertexing
} // namespace o2
