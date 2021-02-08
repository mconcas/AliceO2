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
#include "TrackParametrization.cxx"
#include "TrackParametrizationWithError.cxx"
// #include "DetectorsBase/Propagator.h"
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

// struct SVParamsGPU {
//   // Struct to pass configuration for SVertexer as static symbol on GPU
//   unsigned char useAbsDCA = false;
//   unsigned char propagateToPCA = false;
//   float maxRIni = 0.f;
//   float minParamChange = 0.f;
//   float minRelChi2Change = 0.f;
//   float maxDZIni = 0.f;
//   float maxChi2 = 0.f;
//   float bz = 0.f;
// };

// __constant__ SVParamsGPU gpuStaticConf;

namespace kernels
{
GPUg() void testFitterKernel(o2::vertexing::DCAFitterN<2> fitter2Prong, o2::track::TrackParCov track0, o2::track::TrackParCov track1)
{
  fitter2Prong.setBz(5.0);
  fitter2Prong.setPropagateToPCA(true);  // After finding the vertex, propagate tracks to the DCA. This is default anyway
  fitter2Prong.setMaxR(200);             // do not consider V0 seeds with 2D circles crossing above this R. This is default anyway
  fitter2Prong.setMaxDZIni(4);           // do not consider V0 seeds with tracks Z-distance exceeding this. This is default anyway
  fitter2Prong.setMinParamChange(1e-3);  // stop iterations if max correction is below this value. This is default anyway
  fitter2Prong.setMinRelChi2Change(0.9); // stop iterations if chi2 improves by less that this factor
  fitter2Prong.setUseAbsDCA(true);

  fitter2Prong.print();

  int ncA = fitter2Prong.process(track0, track1); // HERE WE FIT THE VERTICES
  // track1.propagateParamTo(0.1, 5.0);

  printf("End.\n");
}
} // namespace kernels

using PID = o2::track::PID;
using TrackTPCITS = o2::dataformats::TrackTPCITS;
using TrackITS = o2::its::TrackITS;
using TrackTPC = o2::tpc::TrackTPC;

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
}

void SVertexerCUDA::testDCAFitterGPU(std::vector<o2::track::TrackParCov>& trks)
{
  // init and load clusters on card
  // o2::its::gpu::Vector<o2::track::TrackParCov> tr{2, 2};
  // tr.reset(trks.data(), static_cast<int>(trks.size()));
  // auto tr0 = trks[0];
  // auto tr1 = trks[1];

  // call gpu kernel for dca fitter instance
  kernels::testFitterKernel<<<1, 1>>>(mFitter2Prong, trks[0], trks[1]);
  cudaDeviceSynchronize();
  gpuCheckError();
}

} // namespace vertexing
} // namespace o2
