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
/// \brief Creates environment for testing GPU implementation of clusterer in a (for now) minimal and controlled manner.
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]
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

  /**
   * \brief Constructor that initializes ToyProblem with a specific region extractor implementation.
   * \param regionExtractor A unique pointer to the region extraction algorithm to be used.
   */
  ToyProblem(std::unique_ptr<RegionExtractor> regionExtractor);

  ~ToyProblem() = default;
  ToyProblem(const ToyProblem&) = delete;

  /**
   * \brief Sets a new region extractor.
   * \param newRegionExtractor A unique pointer to the new region extractor.
   */
  void setRegionExtractor(std::unique_ptr<RegionExtractor> newRegionExtractor);

  /**
   * \brief Adds chip data to the ToyProblem environment and initiates preprocessing.
   * 
   * In contrast to async version, preprocessing is immediately executed.
   * 
   * \param chipData Pointer to the chip data to be added.
   */
  void addChip(ChipPixelData* chipData);

    /**
   * \brief Adds chip data to the ToyProblem environment.
   * 
   * In contrast to regular version, this version does not initiate preprocessing.
   * 
   * \param chipData Pointer to the chip data to be added.
   */
  void addChipAsync(ChipPixelData* chipData);

  /**
   * \brief Performs clustering on the fully preprocessed chip data.
   * 
   * \param clusterBBoxes Output vector to be populated with bounding boxes of identified clusters.
   * \param clusterPixels Output vector to be populated with pixel data of identified clusters.
   */
  void performClustering(std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);

  /**
   * \brief Executes the toy problem. This includes region extraction/preprocessing and clustering.
   * 
   * Any chips added using addChipAsync will be preprocessed 
   * before the actual clustering is performed using performClustering.
   * 
   * \param clusterBBoxes Output vector to be populated with bounding boxes of identified clusters.
   * \param clusterPixels Output vector to be populated with pixel data of identified clusters.
   */
  void executeToyProblem(std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);

 private:
  std::vector<std::vector<std::vector<int>>> extractedRegions; ///< Stores the preclustered regions for further computation.
  std::vector<std::pair<int,int>> coordinates; ///< Coordinates within chip for each region.
  std::vector<int> chipIds; ///< Chip ID for each region.
  std::unique_ptr<RegionExtractor> regionExtractor; ///< Region extraction strategy.
  std::vector<std::function<void()>> extractionTasks; ///< Stores function calls for addChip operation to measure overhead of region extraction.

  /**
   * \brief Performs preprocessing on chips added using addChipAsync.
   */
  void executeExtractionAsync();

  /**
   * \brief Any postprocessing steps to be executed after the main ToyProblem execution.
   */
  void postProcess();
};

class Timer
{
 public:
  /**
   * \brief Constructor that starts the timer.
   * \param name Name of the timer.
   */
  Timer(const std::string& name)
    : name(name), start(std::chrono::high_resolution_clock::now())
  {
  }

  /**
   * \brief Destructor that stops the timer and prints the elapsed time.
   */
  ~Timer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << name << " took " << elapsed.count() << " seconds\n";
  }

 private:
  std::string name; ///< Name of the timer.
  std::chrono::high_resolution_clock::time_point start; ///< Start time of the timer.
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_TOYPROBLEM_H */
