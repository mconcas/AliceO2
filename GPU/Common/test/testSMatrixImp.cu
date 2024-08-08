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

/// \file testGPUSMatrixImp.cu
/// \author Matteo Concas, Maksym KIzitskyi

#define BOOST_TEST_MODULE Test GPUSMatrixImpl
#ifdef __HIPCC__
#define GPUPLATFORM "HIP"
#include "hip/hip_runtime.h"
#else
#define GPUPLATFORM "CUDA"
#include <cuda.h>
#endif

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <MathUtils/SMatrixGPU.h>
#include <Math/SMatrix.h>
#include <random>

using MatSym3DGPU = o2::math_utils::SMatrixGPU<float, 3, 3, o2::math_utils::MatRepSymGPU<float, 3>>;
using MatSym3D = ROOT::Math::SMatrix<float, 3, 3, ROOT::Math::MatRepSym<float, 3>>;
using Mat3DGPU = o2::math_utils::SMatrixGPU<float, 3, 3, o2::math_utils::MatRepStdGPU<float, 3, 3>>;
using Mat3D = ROOT::Math::SMatrix<float, 3, 3, ROOT::Math::MatRepStd<float, 3, 3>>;

// Macro for checking CUDA errors
#define GPU_CHECK(call)                                                                      \
  do {                                                                                       \
    cudaError_t error = call;                                                                \
    if (error != cudaSuccess) {                                                              \
      fprintf(stderr, "CUDA Error: %s (error code %d)\n", cudaGetErrorString(error), error); \
      return;                                                                                \
    }                                                                                        \
  } while (0)

namespace gpu
{
enum PrintMode {
  Decimal,
  Binary,
  Hexadecimal
};

__device__ void floatToBinaryString(float number, char* buffer)
{
  unsigned char* bytePointer = reinterpret_cast<unsigned char*>(&number);
  for (int byteIndex = 3; byteIndex >= 0; --byteIndex) {
    unsigned char byte = bytePointer[byteIndex];
    for (int bitIndex = 7; bitIndex >= 0; --bitIndex) {
      buffer[(3 - byteIndex) * 8 + (7 - bitIndex)] = (byte & (1 << bitIndex)) ? '1' : '0';
    }
  }
  buffer[32] = '\0'; // Null terminator
}

template <typename MatrixType>
GPUd() void printMatrix(const MatrixType& matrix, const char* name, const PrintMode mode)
{
  if (mode == PrintMode::Binary) {
    char buffer[33];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        floatToBinaryString(matrix(i, j), buffer);
        printf("%s(%d,%d) = %s\n", name, i, j, buffer);
      }
    }
  }
  if (mode == PrintMode::Decimal) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        printf("%s(%i,%i) = %f\n", name, i, j, matrix(i, j));
      }
    }
  }
  if (mode == PrintMode::Hexadecimal) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        printf("%s(%d,%d) = %x\n", name, i, j, o2::gpu::CAMath::Float2UIntReint(matrix(i, j)));
      }
    }
  }
}

// Invert test for symmetric matrix
template <typename T, int D>
__global__ void invertMatrixKernel(T* matrix)
{
  matrix->Invert();
}
} // namespace gpu

// Function to compare two matrices element-wise with a specified tolerance
template <typename MatrixType>
void compareMatricesElementWise(const MatrixType& mat1, const MatrixType& mat2, float tolerance)
{
  auto tol = boost::test_tools::tolerance(tolerance);

  for (unsigned int i = 0; i < mat1.kRows; ++i) {
    for (unsigned int j = 0; j < mat1.kCols; ++j) {
      BOOST_TEST(mat1(i, j) == mat2(i, j), tol);
    }
  }
}

// RAII class for CUDA resources
class GPUMemory
{
 public:
  GPUMemory(size_t size)
  {
    GPU_CHECK(cudaMalloc(&device_ptr, size));
  }
  ~GPUMemory()
  {
    GPU_CHECK(cudaFree(device_ptr));
  }
  void* get() const { return device_ptr; }

 private:
  void* device_ptr;
};

class GPUBenchmark
{
 public:
  GPUBenchmark()
  {
    GPU_CHECK(cudaEventCreate(&startEvent));
    GPU_CHECK(cudaEventCreate(&stopEvent));
  }

  ~GPUBenchmark()
  {
    GPU_CHECK(cudaEventDestroy(startEvent));
    GPU_CHECK(cudaEventDestroy(stopEvent));
  }

  void start()
  {
    GPU_CHECK(cudaEventRecord(startEvent));
  }

  void stop()
  {
    GPU_CHECK(cudaEventRecord(stopEvent));
    GPU_CHECK(cudaEventSynchronize(stopEvent));
    GPU_CHECK(cudaEventElapsedTime(&duration, startEvent, stopEvent));
  }

  float getDuration() const { return duration; }
  void printDuration() const
  {
    std::cout << "Kernel execution time: " << duration << " ms" << std::endl;
  }

 private:
  cudaEvent_t startEvent, stopEvent;
  float duration;
};

template <typename T>
void discardResult(const T&)
{
}

void prologue()
{
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess || !deviceCount) {
    std::cerr << "No " << GPUPLATFORM << " devices found" << std::endl;
    return;
  }

  for (int iDevice = 0; iDevice < deviceCount; ++iDevice) {
    cudaDeviceProp deviceProp;
    discardResult(cudaGetDeviceProperties(&deviceProp, iDevice));
    printf("%s Device %d: %s\n", GPUPLATFORM, iDevice, deviceProp.name);
  }
}

struct GPUSMatrixImplFixtureSolo {
  GPUSMatrixImplFixtureSolo() : SMatrixSym_d(sizeof(MatSym3DGPU)), SMatrixSym_h(), SMatrix_d(sizeof(Mat3DGPU)), SMatrix_h()
  {
    prologue();
    initializeMatrices();
    printMatrixSizes();
  }

  ~GPUSMatrixImplFixtureSolo() = default;
  void initializeMatrices()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(1.0, 10.0);

    // Initialize host matrices with random values
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
        SMatrixSym_h(i, j) = dis(gen);
        SMatrix_h(i, j) = dis(gen);
      }
    }
    SMatrixSym_original_h = SMatrixSym_h;
    SMatrix_original_h = SMatrix_h;

    // Copy host matrices to device
    GPU_CHECK(cudaMemcpy(SMatrixSym_d.get(), &SMatrixSym_h, sizeof(MatSym3DGPU), cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(SMatrix_d.get(), &SMatrix_h, sizeof(Mat3DGPU), cudaMemcpyHostToDevice));
  }

  void printMatrixSizes() const
  {
    printf("sizeof(MatSym3DGPU) = %zu bytes\n", sizeof(MatSym3DGPU));
    printf("sizeof(MatSym3D) = %zu bytes\n", sizeof(MatSym3D));
    printf("sizeof(Mat3DGPU) = %zu bytes\n", sizeof(Mat3DGPU));
    printf("sizeof(Mat3D) = %zu bytes\n", sizeof(Mat3D));
  }

  GPUMemory SMatrixSym_d;
  MatSym3D SMatrixSym_h;
  MatSym3D SMatrixSym_original_h;
  GPUMemory SMatrix_d;
  Mat3D SMatrix_h;
  Mat3D SMatrix_original_h;
};

BOOST_FIXTURE_TEST_CASE(MatrixInversion, GPUSMatrixImplFixtureSolo)
{
  float tolerance = 0.00001f;

  GPUBenchmark benchmark;
  benchmark.start();
  gpu::invertMatrixKernel<MatSym3DGPU, 3><<<1, 1>>>(static_cast<MatSym3DGPU*>(SMatrixSym_d.get()));
  benchmark.stop();
  benchmark.printDuration();
  discardResult(cudaDeviceSynchronize());
  GPU_CHECK(cudaGetLastError());
  GPU_CHECK(cudaMemcpy(&SMatrixSym_h, SMatrixSym_d.get(), sizeof(MatSym3DGPU), cudaMemcpyDeviceToHost));

  MatSym3D identitySym;
  identitySym(0, 0) = 1;
  identitySym(1, 1) = 1;
  identitySym(2, 2) = 1;
  auto operationSym = SMatrixSym_h * SMatrixSym_original_h;
  MatSym3D resultSym;
  ROOT::Math::AssignSym::Evaluate(resultSym, operationSym);
  compareMatricesElementWise(resultSym, identitySym, tolerance);

  benchmark.start();
  gpu::invertMatrixKernel<Mat3DGPU, 3><<<1, 1>>>(static_cast<Mat3DGPU*>(SMatrix_d.get()));
  benchmark.stop();
  benchmark.printDuration();
  discardResult(cudaDeviceSynchronize());
  GPU_CHECK(cudaGetLastError());
  GPU_CHECK(cudaMemcpy(&SMatrix_h, SMatrix_d.get(), sizeof(Mat3DGPU), cudaMemcpyDeviceToHost));

  Mat3D identity;
  identity(0, 0) = 1;
  identity(1, 1) = 1;
  identity(2, 2) = 1;
  auto operation = SMatrix_h * SMatrix_original_h;
  Mat3D result;
  ROOT::Math::Assign<float, 3, 3, decltype(operation), ROOT::Math::MatRepStd<float, 3, 3>, ROOT::Math::MatRepStd<float, 3, 3>>::Evaluate(result, operation);
  compareMatricesElementWise(result, identity, tolerance);
}

// struct GPUSMatrixImplFixtureDuo {
//   GPUSMatrixImplFixtureDuo() : i(3), SMatrixSym_d_A(sizeof(MatSym3DGPU)), SMatrixSym_h_A(), SMatrix_d_A(sizeof(Mat3DGPU)), SMatrix_h_A(), SMatrixSym_d_B(sizeof(MatSym3DGPU)), SMatrixSym_h_B(), SMatrix_d_B(sizeof(Mat3DGPU)), SMatrix_h_B()
//   {
//     prologue();
//     initializeMatrices();
//     printMatrixSizes();
//   }

//   ~GPUSMatrixImplFixtureDuo() = default;

//   void initializeMatrices()
//   {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<float> dis(1.0, 10.0);

//     // Initialize host matrices with random values
//     for (int i = 0; i < 3; ++i) {
//       for (int j = i; j < 3; ++j) {
//         SMatrixSym_h_A(i, j) = dis(gen);
//         SMatrix_h_A(i, j) = dis(gen);

//         SMatrixSym_h_B(i, j) = dis(gen);
//         SMatrix_h_B(i, j) = dis(gen);
//       }
//     }

//     // Copy host matrices to device
//     GPU_CHECK(cudaMemcpy(SMatrixSym_d_A.get(), &SMatrixSym_h_A, sizeof(MatSym3DGPU), cudaMemcpyHostToDevice));
//     GPU_CHECK(cudaMemcpy(SMatrix_d_A.get(), &SMatrix_h_A, sizeof(Mat3DGPU), cudaMemcpyHostToDevice));

//     GPU_CHECK(cudaMemcpy(SMatrixSym_d_B.get(), &SMatrixSym_h_B, sizeof(MatSym3DGPU), cudaMemcpyHostToDevice));
//     GPU_CHECK(cudaMemcpy(SMatrix_d_B.get(), &SMatrix_h_B, sizeof(Mat3DGPU), cudaMemcpyHostToDevice));
//   }

//   void printMatrixSizes() const
//   {
//     printf("sizeof(MatSym3DGPU) = %zu\n", sizeof(MatSym3DGPU));
//     printf("sizeof(MatSym3D) = %zu\n", sizeof(MatSym3D));
//     printf("sizeof(Mat3DGPU) = %zu\n", sizeof(Mat3DGPU));
//     printf("sizeof(Mat3D) = %zu\n", sizeof(Mat3D));
//   }

//   int i;
//   GPUMemory SMatrixSym_d_A;
//   MatSym3D SMatrixSym_h_A;

//   GPUMemory SMatrixSym_d_B;
//   MatSym3D SMatrixSym_h_B;

//   GPUMemory SMatrix_d_A;
//   Mat3D SMatrix_h_A;

//   GPUMemory SMatrix_d_B;
//   Mat3D SMatrix_h_B;
// };

// // Copy test for symmetric matrix
// template <typename T>
// __global__ void copySymMatrixKernel(
//   MatSym3DGPU* srcMatrix,
//   MatSym3DGPU* dstMatrix,
//   const PrintMode mode = PrintMode::Decimal)
// {
//   printf("\nStart copying general matrix\n");
//   printMatrix(*dstMatrix, "Before copying: ", mode);
//   printf("\nCopied values:\n");
//   printMatrix(*srcMatrix, "Copied values: ", mode);
//   printf("\nResult:\n");
//   *dstMatrix = *srcMatrix;
//   printMatrix(*dstMatrix, "After copying: ", mode);
//   printf("\n-------------------------------------------------------\n");
// }

// // Copy test for general matrix
// template <typename T>
// __global__ void copyMatrixKernel(
//   Mat3DGPU* srcMatrix,
//   Mat3DGPU* dstMatrix,
//   const PrintMode mode = PrintMode::Decimal)
// {
//   printf("\nStart copying general matrix\n");
//   printMatrix(*dstMatrix, "Before copying: ", mode);
//   printf("\nCopied values:\n");
//   printMatrix(*srcMatrix, "Copied values: ", mode);
//   printf("\nResult:\n");
//   *dstMatrix = *srcMatrix;
//   printMatrix(*dstMatrix, "After copying: ", mode);
//   printf("\n-------------------------------------------------------\n");
// }

// BOOST_FIXTURE_TEST_CASE(TestMatrixCopyingAndComparison, GPUSMatrixImplFixtureDuo)
// {
//   copySymMatrixKernel<float><<<1, 1>>>(static_cast<MatSym3DGPU*>(SMatrixSym_d_A.get()), static_cast<MatSym3DGPU*>(SMatrixSym_d_B.get()));
//   discardResult(cudaDeviceSynchronize());
//   GPU_CHECK(cudaGetLastError());

//   GPU_CHECK(cudaMemcpy(&SMatrixSym_h_B, SMatrixSym_d_B.get(), sizeof(MatSym3DGPU), cudaMemcpyDeviceToHost));

//   compareMatricesElementWise(SMatrixSym_h_A, SMatrixSym_h_B, 0.0);

//   copyMatrixKernel<float><<<1, 1>>>(static_cast<Mat3DGPU*>(SMatrix_d_A.get()), static_cast<Mat3DGPU*>(SMatrix_d_B.get()));
//   discardResult(cudaDeviceSynchronize());
//   GPU_CHECK(cudaGetLastError());

//   GPU_CHECK(cudaMemcpy(&SMatrix_h_B, SMatrix_d_B.get(), sizeof(Mat3DGPU), cudaMemcpyDeviceToHost));

//   compareMatricesElementWise(SMatrix_h_A, SMatrix_h_B, 0.0);
// }