// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Kernels.cu
/// \author: mconcas@cern.ch

#include "../Shared/Kernels.h"
#include "../Shared/Common.h"
#include <stdio.h>

double bytesToKB(size_t s) { return (double)s / (1024.0); }
double bytesToGB(size_t s) { return (double)s / (1024.0 * 1024.0 * 1024.0); }

#define GPUCHECK(error)                                                                        \
  if (error != cudaSuccess) {                                                                  \
    printf("%serror: '%s'(%d) at %s:%d%s\n", KRED, cudaGetErrorString(error), error, __FILE__, \
           __LINE__, KNRM);                                                                    \
    failed("API returned error code.");                                                        \
  }

namespace o2
{
namespace benchmark
{
namespace gpu
{

///////////////////
/// Kernels go here

template <class buffer_type>
GPUg() void readerKernel(
  size_t Ntimes,
  buffer_type* result,
  buffer_type* buffer,
  size_t bufferSize)
{
  for (size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < bufferSize; i += blockDim.x * gridDim.x) {
    buffer_type result{0};
    result += buffer[i];
  }
}
///////////////////

} // namespace gpu

void printDeviceProp(int deviceId)
{
  const int w1 = 34;
  std::cout << std::left;
  std::cout << std::setw(w1)
            << "--------------------------------------------------------------------------------"
            << std::endl;
  std::cout << std::setw(w1) << "device#" << deviceId << std::endl;

  cudaDeviceProp props;
  GPUCHECK(cudaGetDeviceProperties(&props, deviceId));

  std::cout << std::setw(w1) << "Name: " << props.name << std::endl;
  std::cout << std::setw(w1) << "pciBusID: " << props.pciBusID << std::endl;
  std::cout << std::setw(w1) << "pciDeviceID: " << props.pciDeviceID << std::endl;
  std::cout << std::setw(w1) << "pciDomainID: " << props.pciDomainID << std::endl;
  std::cout << std::setw(w1) << "multiProcessorCount: " << props.multiProcessorCount << std::endl;
  std::cout << std::setw(w1) << "maxThreadsPerMultiProcessor: " << props.maxThreadsPerMultiProcessor
            << std::endl;
  std::cout << std::setw(w1) << "isMultiGpuBoard: " << props.isMultiGpuBoard << std::endl;
  std::cout << std::setw(w1) << "clockRate: " << (float)props.clockRate / 1000.0 << " Mhz" << std::endl;
  std::cout << std::setw(w1) << "memoryClockRate: " << (float)props.memoryClockRate / 1000.0 << " Mhz"
            << std::endl;
  std::cout << std::setw(w1) << "memoryBusWidth: " << props.memoryBusWidth << std::endl;
  std::cout << std::setw(w1) << "clockInstructionRate: " << (float)props.clockRate / 1000.0
            << " Mhz" << std::endl;
  std::cout << std::setw(w1) << "totalGlobalMem: " << std::fixed << std::setprecision(2)
            << bytesToGB(props.totalGlobalMem) << " GB" << std::endl;
#if !defined(__CUDACC__)
  std::cout << std::setw(w1) << "maxSharedMemoryPerMultiProcessor: " << std::fixed << std::setprecision(2)
            << bytesToKB(props.sharedMemPerMultiprocessor) << " KB" << std::endl;
#endif
#if defined(__HIPCC__)
  std::cout << std::setw(w1) << "maxSharedMemoryPerMultiProcessor: " << std::fixed << std::setprecision(2)
            << bytesToKB(props.maxSharedMemoryPerMultiProcessor) << " KB" << std::endl;
#endif
  std::cout << std::setw(w1) << "totalConstMem: " << props.totalConstMem << std::endl;
  std::cout << std::setw(w1) << "sharedMemPerBlock: " << (float)props.sharedMemPerBlock / 1024.0 << " KB"
            << std::endl;
  std::cout << std::setw(w1) << "canMapHostMemory: " << props.canMapHostMemory << std::endl;
  std::cout << std::setw(w1) << "regsPerBlock: " << props.regsPerBlock << std::endl;
  std::cout << std::setw(w1) << "warpSize: " << props.warpSize << std::endl;
  std::cout << std::setw(w1) << "l2CacheSize: " << props.l2CacheSize << std::endl;
  std::cout << std::setw(w1) << "computeMode: " << props.computeMode << std::endl;
  std::cout << std::setw(w1) << "maxThreadsPerBlock: " << props.maxThreadsPerBlock << std::endl;
  std::cout << std::setw(w1) << "maxThreadsDim.x: " << props.maxThreadsDim[0] << std::endl;
  std::cout << std::setw(w1) << "maxThreadsDim.y: " << props.maxThreadsDim[1] << std::endl;
  std::cout << std::setw(w1) << "maxThreadsDim.z: " << props.maxThreadsDim[2] << std::endl;
  std::cout << std::setw(w1) << "maxGridSize.x: " << props.maxGridSize[0] << std::endl;
  std::cout << std::setw(w1) << "maxGridSize.y: " << props.maxGridSize[1] << std::endl;
  std::cout << std::setw(w1) << "maxGridSize.z: " << props.maxGridSize[2] << std::endl;
  std::cout << std::setw(w1) << "major: " << props.major << std::endl;
  std::cout << std::setw(w1) << "minor: " << props.minor << std::endl;
  std::cout << std::setw(w1) << "concurrentKernels: " << props.concurrentKernels << std::endl;
  std::cout << std::setw(w1) << "cooperativeLaunch: " << props.cooperativeLaunch << std::endl;
  std::cout << std::setw(w1) << "cooperativeMultiDeviceLaunch: " << props.cooperativeMultiDeviceLaunch << std::endl;
#if defined(__HIPCC__)
  std::cout << std::setw(w1) << "arch.hasGlobalInt32Atomics: " << props.arch.hasGlobalInt32Atomics << std::endl;
  std::cout << std::setw(w1) << "arch.hasGlobalFloatAtomicExch: " << props.arch.hasGlobalFloatAtomicExch
            << std::endl;
  std::cout << std::setw(w1) << "arch.hasSharedInt32Atomics: " << props.arch.hasSharedInt32Atomics << std::endl;
  std::cout << std::setw(w1) << "arch.hasSharedFloatAtomicExch: " << props.arch.hasSharedFloatAtomicExch
            << std::endl;
  std::cout << std::setw(w1) << "arch.hasFloatAtomicAdd: " << props.arch.hasFloatAtomicAdd << std::endl;
  std::cout << std::setw(w1) << "arch.hasGlobalInt64Atomics: " << props.arch.hasGlobalInt64Atomics << std::endl;
  std::cout << std::setw(w1) << "arch.hasSharedInt64Atomics: " << props.arch.hasSharedInt64Atomics << std::endl;
  std::cout << std::setw(w1) << "arch.hasDoubles: " << props.arch.hasDoubles << std::endl;
  std::cout << std::setw(w1) << "arch.hasWarpVote: " << props.arch.hasWarpVote << std::endl;
  std::cout << std::setw(w1) << "arch.hasWarpBallot: " << props.arch.hasWarpBallot << std::endl;
  std::cout << std::setw(w1) << "arch.hasWarpShuffle: " << props.arch.hasWarpShuffle << std::endl;
  std::cout << std::setw(w1) << "arch.hasFunnelShift: " << props.arch.hasFunnelShift << std::endl;
  std::cout << std::setw(w1) << "arch.hasThreadFenceSystem: " << props.arch.hasThreadFenceSystem << std::endl;
  std::cout << std::setw(w1) << "arch.hasSyncThreadsExt: " << props.arch.hasSyncThreadsExt << std::endl;
  std::cout << std::setw(w1) << "arch.hasSurfaceFuncs: " << props.arch.hasSurfaceFuncs << std::endl;
  std::cout << std::setw(w1) << "arch.has3dGrid: " << props.arch.has3dGrid << std::endl;
  std::cout << std::setw(w1) << "arch.hasDynamicParallelism: " << props.arch.hasDynamicParallelism << std::endl;
  std::cout << std::setw(w1) << "gcnArchName: " << props.gcnArchName << std::endl;
#endif
  std::cout << std::setw(w1) << "isIntegrated: " << props.integrated << std::endl;
  std::cout << std::setw(w1) << "maxTexture1D: " << props.maxTexture1D << std::endl;
  std::cout << std::setw(w1) << "maxTexture2D.width: " << props.maxTexture2D[0] << std::endl;
  std::cout << std::setw(w1) << "maxTexture2D.height: " << props.maxTexture2D[1] << std::endl;
  std::cout << std::setw(w1) << "maxTexture3D.width: " << props.maxTexture3D[0] << std::endl;
  std::cout << std::setw(w1) << "maxTexture3D.height: " << props.maxTexture3D[1] << std::endl;
  std::cout << std::setw(w1) << "maxTexture3D.depth: " << props.maxTexture3D[2] << std::endl;
#if defined(__HIPCC__)
  std::cout << std::setw(w1) << "isLargeBar: " << props.isLargeBar << std::endl;
  std::cout << std::setw(w1) << "asicRevision: " << props.asicRevision << std::endl;
#endif

  int deviceCnt;
  GPUCHECK(cudaGetDeviceCount(&deviceCnt));
  std::cout << std::setw(w1) << "peers: ";
  for (int i = 0; i < deviceCnt; i++) {
    int isPeer;
    GPUCHECK(cudaDeviceCanAccessPeer(&isPeer, i, deviceId));
    if (isPeer) {
      std::cout << "device#" << i << " ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(w1) << "non-peers: ";
  for (int i = 0; i < deviceCnt; i++) {
    int isPeer;
    GPUCHECK(cudaDeviceCanAccessPeer(&isPeer, i, deviceId));
    if (!isPeer) {
      std::cout << "device#" << i << " ";
    }
  }
  std::cout << std::endl;

  size_t free, total;
  GPUCHECK(cudaMemGetInfo(&free, &total));

  std::cout << std::fixed << std::setprecision(2);
  std::cout << std::setw(w1) << "memInfo.total: " << bytesToGB(total) << " GB" << std::endl;
  std::cout << std::setw(w1) << "memInfo.free:  " << bytesToGB(free) << " GB (" << std::setprecision(0)
            << (float)free / total * 100.0 << "%)" << std::endl;
}

template <class buffer_type>
template <typename... T>
float GPUbenchmark<buffer_type>::measure(void (GPUbenchmark<buffer_type>::*task)(T...), const char* taskName, T&&... args)
{
  float diff{0.f};
  std::cout << std::setw(2) << ">>> " << taskName;
  auto start = std::chrono::high_resolution_clock::now();
  (this->*task)(std::forward<T>(args)...);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> diff_t{end - start};
  diff = diff_t.count();
  std::cout << std::setw(2) << " completed in: \x1B[32m" << diff << " ms\x1B[0m" << std::endl;
  return diff;
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::printDevices()
{
  int deviceCnt;
  GPUCHECK(cudaGetDeviceCount(&deviceCnt));

  for (int i = 0; i < deviceCnt; i++) {
    GPUCHECK(cudaSetDevice(i));
    printDeviceProp(i);
  }
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::generalInit(const int deviceId)
{
  cudaDeviceProp props;
  size_t free;

  // Fetch and store features
  GPUCHECK(cudaGetDeviceProperties(&props, deviceId));
  GPUCHECK(cudaMemGetInfo(&free, &mState.totalMemory));

  mState.nMultiprocessors = props.multiProcessorCount;
  mState.nMaxThreadsPerBlock = props.maxThreadsPerMultiProcessor;
  mState.nMaxThreadsPerDimension = props.maxThreadsDim[0];
  mState.scratchSize = static_cast<long int>(FREE_MEMORY_FRACTION_TO_ALLOCATE * free);
  std::cout << ">>> Running on: " << props.name << std::endl;
  // Allocate scratch on GPU
  GPUCHECK(cudaMalloc(reinterpret_cast<void**>(&mState.scratchPtr), mState.scratchSize));
  mState.computeScratchPtrs();

  // Initialize corresponding buffers on host and copy content on GPU
  mState.getHostBuffers().resize(mState.getScratchPtrs().size());
  for (size_t iScratchPart{0}; iScratchPart < mState.getScratchPtrs().size(); ++iScratchPart) {
    mState.getHostBuffers()[iScratchPart].resize(gpuState<buffer_type>::getArraySize());
    GPUCHECK(cudaMemcpy(mState.getScratchPtrs()[iScratchPart], mState.getHostBuffers()[iScratchPart].data(), gpuState<buffer_type>::getArraySize() * sizeof(buffer_type), cudaMemcpyHostToDevice));
  }
  std::cout << "    ├ Allocated: " << std::setprecision(2) << bytesToGB(mState.scratchSize) << "/" << std::setprecision(2) << bytesToGB(mState.totalMemory)
            << "(GB) [" << std::setprecision(3) << (100.f) * (mState.scratchSize / (float)mState.totalMemory) << "%]\n"
            << "    ├ Number of scratch partitions: " << mState.getMaxSegments() << " of " << PARTITION_SIZE_GB << "GB each\n"
            << "    ├ Size of arrays in segments: " << gpuState<buffer_type>::getArraySize() << " elements" << std::endl
            << "    └ Memory buffers copied from host to device"
            << std::endl;
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::readingInit()
{
  mState.deviceReadingResultsPtrs.resize(mState.getMaxSegments());
  mState.hostReadingResultsVector.resize(mState.getMaxSegments());
  for (size_t iScratchPart{0}; iScratchPart < mState.getMaxSegments(); ++iScratchPart) {
    GPUCHECK(cudaMalloc(reinterpret_cast<void**>(&(mState.deviceReadingResultsPtrs[iScratchPart])), sizeof(buffer_type)));
  }
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::readingBenchmark()
{
  dim3 nBlocks{static_cast<uint32_t>(mState.nMultiprocessors / mState.getMaxSegments())};
  dim3 nThreads{static_cast<uint32_t>(std::min(mState.nMaxThreadsPerDimension, mState.nMaxThreadsPerBlock))};

  for (size_t iScratchPart{0}; iScratchPart < mState.getMaxSegments(); ++iScratchPart) {
    gpu::readerKernel<buffer_type><<<nBlocks, nThreads>>>(1, mState.deviceReadingResultsPtrs[iScratchPart], mState.getScratchPtrs()[iScratchPart], gpuState<buffer_type>::getArraySize());
  }
  GPUCHECK(cudaDeviceSynchronize());
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::readingFinalize()
{
  for (size_t iScratchPart{0}; iScratchPart < mState.getMaxSegments(); ++iScratchPart) {
    GPUCHECK(cudaMemcpy(&mState.hostReadingResultsVector[iScratchPart], mState.deviceReadingResultsPtrs[iScratchPart], sizeof(buffer_type), cudaMemcpyDeviceToHost));
    std::cout << "result " << iScratchPart << ": " << mState.hostReadingResultsVector[iScratchPart] << std::endl;
  }
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::generalFinalize()
{
  GPUCHECK(cudaFree(mState.scratchPtr));
}

template <class buffer_type>
void GPUbenchmark<buffer_type>::run()
{
  // printDevices();
  generalInit(0);

  readingInit();
  measure(&GPUbenchmark<buffer_type>::readingBenchmark, "Reading benchmark");
  GPUbenchmark<buffer_type>::generalFinalize();
}

template class GPUbenchmark<char>;
// template class GPUbenchmark<uint4>;
template class GPUbenchmark<size_t>;
template class GPUbenchmark<int>;

} // namespace benchmark
} // namespace o2