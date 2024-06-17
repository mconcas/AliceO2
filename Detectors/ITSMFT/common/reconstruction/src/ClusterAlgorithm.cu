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

/// \file ClusterAlgorithm.cu
/// \brief Implementation of the Playne CCL Algorithm for Clustering on the GPU
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]

#include <cuda.h>
#include <cuda_runtime.h>

#include <chrono>

#include "ITSMFTReconstruction/ClusterAlgorithm.h"
#include "ITSMFTReconstruction/BoundingBox.h"

using namespace o2::itsmft;

#define CHECK_CUDA_ERROR(err)                                       \
  do {                                                              \
    if (err != cudaSuccess) {                                       \
      fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err)); \
    }                                                               \
  } while (0)

namespace o2::itsmft::gpu
{
using Point = o2::itsmft::Point;
using MinimalistBoundingBox = o2::itsmft::MinimalistBoundingBox;
using BoundingBox = o2::itsmft::BoundingBox;

__device__ void print(int* data, int* labels, int* parent, int* rank, int nrow, int ncol)
{
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      // Print elements of the Data Matrix
      printf("%d ", data[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Labels Matrix
      printf("%d ", labels[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Parent Matrix
      printf("%d ", parent[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Rank Matrix
      printf("%d ", rank[i * ncol + j]);
    }
    printf("\n");
  }
  printf("\n");
}

// recursively find root of cluster
__device__ int find(int x, int* parent)
{
  if (x != parent[x]) {
    parent[x] = find(parent[x], parent);
  }
  return parent[x];
}

// merge two clusters
__device__ int unify(int x, int y, int* parent, int* rank)
{
  int rootX = find(x, parent);
  int rootY = find(y, parent);

  if (rootX == rootY)
    return;

  if (rank[rootX] < rank[rootY]) {
    parent[rootX] = rootY;
    return rootY;
  } else {
    parent[rootY] = rootX;
    if (rank[rootX] == rank[rootY]) {
      rank[rootX]++;
    }
    return rootX;
  }
}

__global__ void ccl_kernel(int N, int* data, int* coordinates, int* regionSizes, int* regionHeights, int* regionWidths,
                           int* startIndices, int* labels, int* parent, int* rank, int* numClusters,
                           int* clusterSizes, Point* scratchMemory, int* scratchMemIndex, int* stridedSizes, 
                           int* stridedPositions, int stride, MinimalistBoundingBox* boxes, MinimalistBoundingBox* stridedBoxes, 
                           int* chipIds, int* stridedChipIds)
{
  int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadIndex >= N) {
    return;
  }

  for (int regionIdx = threadIndex; regionIdx < N; regionIdx += blockDim.x * gridDim.x) {
    int startIndex = startIndices[regionIdx];
    int width = regionWidths[regionIdx];
    int height = regionHeights[regionIdx];
    int regionStartRow = coordinates[2 * regionIdx];
    int regionStartCol = coordinates[2 * regionIdx + 1];

    int currentLabel = 1;

    // First pass
    for (int r = 0; r < height; ++r) {
      for (int c = 0; c < width; ++c) {
        int idx = startIndex + r * width + c;
        if (data[idx] == 1) {
          int leftIdx = (c == 0 ? -1 : idx - 1);
          int topIdx = (r == 0 ? -1 : idx - width);

          bool connectLeft = leftIdx != -1 && data[leftIdx] == 1;
          bool connectTop = topIdx != -1 && data[topIdx] == 1;

          if (!connectLeft && !connectTop) {
            labels[idx] = currentLabel++;
            parent[idx] = idx;
            rank[idx] = 0;
            numClusters[regionIdx]++;
            clusterSizes[idx] = 1;
            boxes[idx].min_r = r + regionStartRow;
            boxes[idx].min_c = c + regionStartCol;
            boxes[idx].max_r = r + regionStartRow;
            boxes[idx].max_c = c + regionStartCol;
          } else if (connectLeft && connectTop) {
            int topRootIdx = find(topIdx, parent);
            int leftRootIdx = find(leftIdx, parent);

            if (topRootIdx == leftRootIdx) {
              labels[idx] = labels[topRootIdx];
              parent[idx] = topRootIdx;
              rank[idx] = 0;

              // in this case, the bounding box already encompasses the newly added pixel
            } else {
              int newRoot = unify(leftRootIdx, topRootIdx, parent, rank);
              labels[idx] = labels[newRoot];
              parent[idx] = newRoot;
              numClusters[regionIdx]--;
              clusterSizes[newRoot] = clusterSizes[leftRootIdx] + clusterSizes[topRootIdx];

              boxes[newRoot].min_r = min(boxes[topRootIdx].min_r, boxes[leftRootIdx].min_r);
              boxes[newRoot].min_c = min(boxes[topRootIdx].min_c, boxes[leftRootIdx].min_c);
              boxes[newRoot].max_r = max(boxes[topRootIdx].max_r, boxes[leftRootIdx].max_r);
              boxes[newRoot].max_c = max(boxes[topRootIdx].max_c, boxes[leftRootIdx].max_c);
            }
          } else {
            int rootIdx;
            if (connectLeft) {
              rootIdx = find(leftIdx, parent);
            } else if (connectTop) {
              rootIdx = find(topIdx, parent);
            }
            labels[idx] = labels[rootIdx];
            parent[idx] = rootIdx;
            rank[idx] = 0;
            clusterSizes[rootIdx]++;

            boxes[rootIdx].min_r = min(boxes[rootIdx].min_r, r + regionStartRow);
            boxes[rootIdx].min_c = min(boxes[rootIdx].min_c, c + regionStartCol);
            boxes[rootIdx].max_r = max(boxes[rootIdx].max_r, r + regionStartRow);
            boxes[rootIdx].max_c = max(boxes[rootIdx].max_c, c + regionStartCol);
          }
        }
      }
    }

    int localClusterIndex = 0;

    // Second pass
    for (int r = 0; r < height; ++r) {
      for (int c = 0; c < width; ++c) {
        int idx = startIndex + r * width + c;
        if (data[idx] != 0) {
          int rootIdx = find(idx, parent);
          labels[idx] = labels[rootIdx];

          if (idx == rootIdx) {
            int currentClusterSize = clusterSizes[idx];
            MinimalistBoundingBox currentClusterBox = boxes[idx];

            // global scratch memory index
            int reservedIndex = atomicAdd(scratchMemIndex, currentClusterSize);

            // pseudo local strided index
            int stridedIndex = stride * threadIndex + localClusterIndex;
            stridedSizes[stridedIndex] = currentClusterSize;
            stridedBoxes[stridedIndex] = currentClusterBox;
            stridedPositions[stridedIndex] = reservedIndex;
            stridedChipIds[stridedIndex] = chipIds[regionIdx];

            int localPixelIndex = 0;
            for (int r2 = 0; r2 < height; ++r2) {
              for (int c2 = 0; c2 < width; ++c2) {
                int idx2 = startIndex + r2 * width + c2;
                if (data[idx2] != 0 && find(idx2, parent) == rootIdx) {
                  scratchMemory[reservedIndex + localPixelIndex].r = r2 + regionStartRow;
                  scratchMemory[reservedIndex + localPixelIndex].c = c2 + regionStartCol;
                  localPixelIndex++;
                }
              }
            }
          }
          localClusterIndex++;
        }
      }
    }
  }
}
} // namespace o2::itsmft::gpu

std::vector<int> flatten(const std::vector<std::vector<std::vector<int>>>& data)
{
  std::vector<int> flatData;
  for (const auto& region : data)
    for (const auto& row : region)
      for (const auto& pixel : row)
        flatData.push_back(pixel);  
  return flatData;
}

std::vector<int> flattenCoordinates(const std::vector<std::pair<int,int>>& coordinates)
{
  std::vector<int> flatCoordinates;
  for (const auto& p : coordinates)
  {
    flatCoordinates.push_back(p.first);  
    flatCoordinates.push_back(p.second); 
  } 
  return flatCoordinates;
}

void ClusterAlgorithm::clusterize(const std::vector<std::vector<std::vector<int>>>& data, const std::vector<int>& chipIds, const std::vector<std::pair<int,int>>& coordinates, 
                                  std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  auto start_all = std::chrono::high_resolution_clock::now();
  auto start_pre = std::chrono::high_resolution_clock::now();

  std::vector<int> flatData = flatten(data);
  std::vector<int> flatCoordinates = flattenCoordinates(coordinates); 

  int totalNumPixels = flatData.size();
  int numRegions = data.size();
  int scratchMemLength = 100 * numRegions; // assuming that there will not be more than 100 pixels in clusters per region on average
  int stride = 8;

  std::vector<int> regionSizes;
  std::vector<int> regionHeights;
  std::vector<int> regionWidths;
  std::vector<int> startIndices;
  int* parent = new int[totalNumPixels]();

  Point* scratchMemory = new Point[scratchMemLength]();
  int* stridedSizes = new int[numRegions * stride]();
  int* stridedPositions = new int[numRegions * stride]();
  int scratchMemIndex = 0;
  MinimalistBoundingBox* boxes = new MinimalistBoundingBox[totalNumPixels]();
  MinimalistBoundingBox* stridedBoxes = new MinimalistBoundingBox[numRegions * stride]();
  int* stridedChipIds = new int[numRegions * stride]();

  auto start_dataprep = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < totalNumPixels; i++) {
    parent[i] = i;
  }

  startIndices.push_back(0);

  for (const auto& region : data) {
    regionSizes.push_back(region.size() * (region.empty() ? 0 : region[0].size()));
    regionHeights.push_back(region.size());
    regionWidths.push_back((region.empty() ? 0 : region[0].size()));
    startIndices.push_back(startIndices.back() + regionSizes.back());
  }

  // Remove the last start index (beyond the end of the data)
  startIndices.pop_back();

  auto stop_dataprep = std::chrono::high_resolution_clock::now();

  int* deviceData;
  int* deviceChipIds;
  int* deviceSizes;
  int* deviceHeights;
  int* deviceWidths;
  int* deviceStartIndices;
  int* deviceLabels;
  int* deviceParents;
  int* deviceRanks;
  int* deviceClusterSizes;
  int* deviceNumClusters;
  int* deviceCoordinates;

  Point* deviceScratchMemory;
  int* deviceScratchMemIndex;
  int* deviceStridedSizes;
  int* deviceStridedPositions;
  int* deviceStridedChipIds;
  MinimalistBoundingBox* deviceBoxes;
  MinimalistBoundingBox* deviceStridedBoxes;

  auto start_malloc = std::chrono::high_resolution_clock::now();

  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceData, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceCoordinates, flatCoordinates.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceChipIds, numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceSizes, regionSizes.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceHeights, regionHeights.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceWidths, regionWidths.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStartIndices, startIndices.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceLabels, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceParents, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceRanks, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceClusterSizes, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceNumClusters, numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceScratchMemory, scratchMemLength * sizeof(Point)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceScratchMemIndex, sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedSizes, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedPositions, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedChipIds, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceBoxes, totalNumPixels * sizeof(MinimalistBoundingBox)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedBoxes, stride * numRegions * sizeof(MinimalistBoundingBox)));

  auto stop_malloc = std::chrono::high_resolution_clock::now();

  // Add corresponding amount of bytes for each of the cudaMallocs from above
  int totalMallocBytes = 0;
  totalMallocBytes += flatCoordinates.size() * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(int);
  totalMallocBytes += numRegions * sizeof(int);
  totalMallocBytes += regionSizes.size() * sizeof(int);
  totalMallocBytes += regionHeights.size() * sizeof(int);
  totalMallocBytes += regionWidths.size() * sizeof(int);
  totalMallocBytes += startIndices.size() * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(int);
  totalMallocBytes += numRegions * sizeof(int);
  totalMallocBytes += scratchMemLength * sizeof(Point);
  totalMallocBytes += sizeof(int);
  totalMallocBytes += stride * numRegions * sizeof(int);
  totalMallocBytes += stride * numRegions * sizeof(int);
  totalMallocBytes += stride * numRegions * sizeof(int);
  totalMallocBytes += totalNumPixels * sizeof(MinimalistBoundingBox);
  totalMallocBytes += stride * numRegions * sizeof(MinimalistBoundingBox);

  auto start_memcpy = std::chrono::high_resolution_clock::now();

  CHECK_CUDA_ERROR(cudaMemcpy(deviceData, flatData.data(), totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceCoordinates, flatCoordinates.data(), flatCoordinates.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceChipIds, chipIds.data(), numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceSizes, regionSizes.data(), regionSizes.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceWidths, regionWidths.data(), regionWidths.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceHeights, regionHeights.data(), regionHeights.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStartIndices, startIndices.data(), startIndices.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceParents, parent, totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));

  auto stop_memcpy = std::chrono::high_resolution_clock::now();

  // Add corresponding amount of bytes for each of the cudaMemcpys from above
  int totalMemcpyBytes = 0;
  totalMemcpyBytes += flatCoordinates.size() * sizeof(int);
  totalMemcpyBytes += totalNumPixels * sizeof(int);
  totalMemcpyBytes += numRegions * sizeof(int);
  totalMemcpyBytes += regionSizes.size() * sizeof(int);
  totalMemcpyBytes += regionHeights.size() * sizeof(int);
  totalMemcpyBytes += regionWidths.size() * sizeof(int);
  totalMemcpyBytes += startIndices.size() * sizeof(int);
  totalMemcpyBytes += totalNumPixels * sizeof(int);

  auto start_memset = std::chrono::high_resolution_clock::now();

  CHECK_CUDA_ERROR(cudaMemset(deviceLabels, 0, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceRanks, 0, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceClusterSizes, 0, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceNumClusters, 0, numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceScratchMemory, 0, scratchMemLength * sizeof(Point)));
  CHECK_CUDA_ERROR(cudaMemset(deviceScratchMemIndex, 0, sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceStridedSizes, 0, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceStridedPositions, 0, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceStridedChipIds, 0, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMemset(deviceBoxes, 0, totalNumPixels * sizeof(MinimalistBoundingBox)));
  CHECK_CUDA_ERROR(cudaMemset(deviceStridedBoxes, 0, stride * numRegions * sizeof(MinimalistBoundingBox)));

  auto stop_memset = std::chrono::high_resolution_clock::now();

  // Add corresponding amount of bytes for each of the cudaMemsets from above
  int totalMemsetBytes = 0;
  totalMemsetBytes += totalNumPixels * sizeof(int);
  totalMemsetBytes += totalNumPixels * sizeof(int);
  totalMemsetBytes += totalNumPixels * sizeof(int);
  totalMemsetBytes += numRegions * sizeof(int);
  totalMemsetBytes += scratchMemLength * sizeof(Point);
  totalMemsetBytes += sizeof(int);
  totalMemsetBytes += stride * numRegions * sizeof(int);
  totalMemsetBytes += stride * numRegions * sizeof(int);
  totalMemsetBytes += stride * numRegions * sizeof(int);
  totalMemsetBytes += totalNumPixels * sizeof(MinimalistBoundingBox);
  totalMemsetBytes += stride * numRegions * sizeof(MinimalistBoundingBox);

  const int threadsPerBlock = 256;
  const int numBlocks = (numRegions + threadsPerBlock - 1) / threadsPerBlock;

  auto end_pre = std::chrono::high_resolution_clock::now();

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  gpu::ccl_kernel<<<numBlocks, threadsPerBlock>>>(numRegions, deviceData, deviceCoordinates, deviceSizes, deviceHeights,
                                                  deviceWidths, deviceStartIndices, deviceLabels, deviceParents, deviceRanks,
                                                  deviceNumClusters, deviceClusterSizes, deviceScratchMemory, deviceScratchMemIndex,
                                                  deviceStridedSizes, deviceStridedPositions, stride, deviceBoxes, deviceStridedBoxes, 
                                                  deviceChipIds, deviceStridedChipIds);
  cudaEventRecord(stop);

  auto start_post = std::chrono::high_resolution_clock::now();

  CHECK_CUDA_ERROR(cudaMemcpy(scratchMemory, deviceScratchMemory, scratchMemLength * sizeof(int), cudaMemcpyDeviceToHost)); // copying can be shortened for tiny optimization
  CHECK_CUDA_ERROR(cudaMemcpy(stridedSizes, deviceStridedSizes, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedPositions, deviceStridedPositions, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedChipIds, deviceStridedChipIds, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedBoxes, deviceStridedBoxes, stride * numRegions * sizeof(MinimalistBoundingBox), cudaMemcpyDeviceToHost));

  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  std::cout << "GPU KERNEL TOOK " << milliseconds << "ms" << std::endl;

  for (int clusterIdx = 0; clusterIdx < stride * numRegions; ++clusterIdx) {
    if (stridedSizes[clusterIdx] == 0)
      continue;

    clusterPixels.emplace_back();
    int clusterStartPosition = stridedPositions[clusterIdx];
    int chipId = stridedChipIds[clusterIdx];

    const MinimalistBoundingBox& sBox = stridedBoxes[clusterIdx];
    BoundingBox bBox(chipId);
    bBox.rowMin = static_cast<uint16_t>(sBox.min_r);
    bBox.colMin = static_cast<uint16_t>(sBox.min_c);
    bBox.rowMax = static_cast<uint16_t>(sBox.max_r);
    bBox.colMax = static_cast<uint16_t>(sBox.max_c);
    clusterBBoxes.push_back(bBox);

    for (int pixelIdx = 0; pixelIdx < stridedSizes[clusterIdx]; ++pixelIdx) {
      clusterPixels.back().push_back(PixelData(scratchMemory[clusterStartPosition + pixelIdx].r, scratchMemory[clusterStartPosition + pixelIdx].c));
    }
  }

  delete[] parent;
  delete[] scratchMemory;
  delete[] stridedSizes;
  delete[] stridedPositions;
  delete[] stridedChipIds;
  delete[] boxes;
  delete[] stridedBoxes;

  CHECK_CUDA_ERROR(cudaFree(deviceData));
  CHECK_CUDA_ERROR(cudaFree(deviceCoordinates));
  CHECK_CUDA_ERROR(cudaFree(deviceChipIds));
  CHECK_CUDA_ERROR(cudaFree(deviceSizes));
  CHECK_CUDA_ERROR(cudaFree(deviceHeights));
  CHECK_CUDA_ERROR(cudaFree(deviceWidths));
  CHECK_CUDA_ERROR(cudaFree(deviceStartIndices));
  CHECK_CUDA_ERROR(cudaFree(deviceLabels));
  CHECK_CUDA_ERROR(cudaFree(deviceParents));
  CHECK_CUDA_ERROR(cudaFree(deviceRanks));
  CHECK_CUDA_ERROR(cudaFree(deviceClusterSizes));
  CHECK_CUDA_ERROR(cudaFree(deviceNumClusters));
  CHECK_CUDA_ERROR(cudaFree(deviceScratchMemory));
  CHECK_CUDA_ERROR(cudaFree(deviceScratchMemIndex));
  CHECK_CUDA_ERROR(cudaFree(deviceStridedSizes));
  CHECK_CUDA_ERROR(cudaFree(deviceStridedPositions));
  CHECK_CUDA_ERROR(cudaFree(deviceStridedChipIds));
  CHECK_CUDA_ERROR(cudaFree(deviceBoxes));
  CHECK_CUDA_ERROR(cudaFree(deviceStridedBoxes));

  auto end_post = std::chrono::high_resolution_clock::now();
  auto end_all = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_all - start_all);
  std::cout << "Time taken for all: " << duration.count() / 1000 << " ms" << std::endl;

  duration = std::chrono::duration_cast<std::chrono::microseconds>(end_pre - start_pre);
  std::cout << "Time taken for pre: " << duration.count() / 1000 << " ms" << std::endl;

  duration = std::chrono::duration_cast<std::chrono::microseconds>(end_post - start_post);
  std::cout << "Time taken for post: " << duration.count() / 1000 << " ms" << std::endl;

  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_dataprep - start_dataprep);
  std::cout << "Time taken for dataprep: " << duration.count() / 1000 << " ms" << std::endl;

  std::cout << std::endl;

  auto malloctime = std::chrono::duration_cast<std::chrono::microseconds>(stop_malloc - start_malloc);
  std::cout << "Time taken for MALLOC: " << malloctime.count() << " us" << std::endl;
  std::cout << "Bytes allocated using MALLOC: " << totalMallocBytes << std::endl;
  double throughput_malloc = static_cast<double>(totalMallocBytes) / (malloctime.count() / 1e6) / 1e9;
  std::cout << "Resulting Throughput (GB/s): " << throughput_malloc << std::endl;
  std::cout << std::endl;

  auto memcpytime = std::chrono::duration_cast<std::chrono::microseconds>(stop_memcpy - start_memcpy);
  std::cout << "Time taken for MEMCPY: " << memcpytime.count() << " us" << std::endl;
  std::cout << "Bytes set using MEMCPY: " << totalMemcpyBytes << std::endl;
  double throughput_memcpy = static_cast<double>(totalMemcpyBytes) / (memcpytime.count() / 1e6) / 1e9;
  std::cout << "Resulting Throughput (GB/s): " << throughput_memcpy << std::endl;
  std::cout << std::endl;

  auto memsettime = std::chrono::duration_cast<std::chrono::microseconds>(stop_memset - start_memset);
  std::cout << "Time taken for MEMSET: " << memsettime.count() << " us" << std::endl;
  std::cout << "Bytes set using MEMSET: " << totalMemsetBytes << std::endl;
  double throughput_memset = static_cast<double>(totalMemsetBytes) / (memsettime.count() / 1e6) / 1e9;
  std::cout << "Resulting Throughput (GB/s): " << throughput_memset << std::endl;
  std::cout << std::endl;
}