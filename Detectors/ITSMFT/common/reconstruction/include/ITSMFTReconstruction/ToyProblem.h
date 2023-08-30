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

/// \file ToyProblem.h
/// \brief Preparation for GPU algorithm testing suitef
#ifndef ALICEO2_ITS_TOYPROBLEM_H
#define ALICEO2_ITS_TOYPROBLEM_H

#define _PERFORM_TIMING_

// MAX DISTANCE FOR REGION EXTRACTION IN ROWS
#define MAX_DIST_X 3
// MAX DISTANCE FOR REGION EXTRACTION IN COLUMNS
#define MAX_DIST_Y 3

#include <utility>
#include <vector>
#include <cstring>
#include <memory>
#include <chrono>
#include <iostream>
#include <functional>
#include "ITSMFTReconstruction/PixelData.h"
#include "ITSMFTReconstruction/RegionExtractor.h"
// #include "ITSMFTReconstruction/ClusterAlgorithm.h"
#include "ITSMFTReconstruction/BoundingBox.h"

#ifdef _PERFORM_TIMING_
#include <TStopwatch.h>
#endif

namespace o2
{
namespace itsmft
{
class ToyProblem
{
 public:
  ToyProblem(std::unique_ptr<RegionExtractor> regionExtractor);
  ~ToyProblem() = default;
  ToyProblem(const ToyProblem&) = delete;

  void setRegionExtractor(std::unique_ptr<RegionExtractor> newRegionExtractor);
  // void setClusterAlgorithm(std::unique_ptr<ClusterAlgorithm> newClusterAlgorithm);

  void addChip(ChipPixelData* chipData);
  void addChipAsync(ChipPixelData* chipData);
  void performClustering(std::vector<BoundingBox2>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);
  void executeToyProblem(std::vector<BoundingBox2>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);

 private:
  // stores the "preclustered" regions that will be used for further computation with some algorithm of choice
  std::vector<std::vector<std::vector<int>>> extractedRegions;
  // store coordinate within chip for each region
  std::vector<std::pair<int,int>> coordinates;
  // store chip ID for each region
  std::vector<int> chipIds;
  // region extraction strategy
  std::unique_ptr<RegionExtractor> regionExtractor;
  // clustering strategy
  // std::unique_ptr<ClusterAlgorithm> clusterAlgorithm;
  // stores function calls for addChip operation to measure overhead of region extraction
  std::vector<std::function<void()>> extractionTasks;

  void executeExtractionAsync();
  void postProcess();
};

class Timer
{
 public:
  Timer(const std::string& name)
    : name(name), start(std::chrono::high_resolution_clock::now())
  {
  }

  ~Timer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << name << " took " << elapsed.count() << " seconds\n";
  }

 private:
  std::string name;
  std::chrono::high_resolution_clock::time_point start;
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_TOYPROBLEM_H */
