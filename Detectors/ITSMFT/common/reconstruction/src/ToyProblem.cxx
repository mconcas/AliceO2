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

/// \file ToyProblem.cxx
/// \brief Implementation of the ITS cluster finder toy problem
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]

#include "ITSMFTReconstruction/RegionExtractor.h"
#include "ITSMFTReconstruction/BoundingBox.h"
#include "ITSMFTReconstruction/ToyProblem.h"
#include "ITSMFTReconstruction/ClusterAlgorithm.h"

using namespace o2::itsmft;

ToyProblem::ToyProblem(std::unique_ptr<RegionExtractor> regionExtractor)
  : regionExtractor(std::move(regionExtractor))
{
}

void ToyProblem::addChip(ChipPixelData* chipData)
{
  if (!regionExtractor) {
    throw std::runtime_error("Clustering algorithm is not set");
  }
  ExtractedRegionsWrapper chipRegionsWrapped = regionExtractor->preprocess(chipData, MAX_DIST_X, MAX_DIST_Y);
  std::vector<std::vector<std::vector<int>>> chipRegions = chipRegionsWrapped.regions;
  std::vector<std::pair<int, int>> chipRegionCoordinates = chipRegionsWrapped.coordinates;
  extractedRegions.insert(extractedRegions.end(), chipRegions.begin(), chipRegions.end());
  chipIds.push_back(chipData->getChipID());
  coordinates.insert(coordinates.end(), chipRegionCoordinates.begin(), chipRegionCoordinates.end());
}

void ToyProblem::addChipAsync(ChipPixelData* chipData)
{
  extractionTasks.push_back([this, chipData] { addChip(chipData); });
}

void ToyProblem::executeExtractionAsync()
{
  if (extractionTasks.empty())
    return;
  std::cout << "Number of tasks: " << extractionTasks.size() << std::endl;
  Timer timer("executeExtractionAsync");
  for (auto& task : extractionTasks) {
    task();
  }
}

void ToyProblem::executeToyProblem(std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  Timer timer("executeToyProblem");
  executeExtractionAsync();
  performClustering(clusterBBoxes, clusterPixels);
}

void ToyProblem::performClustering(std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  Timer timer("performClustering");
  ClusterAlgorithm algo;
  algo.clusterize(extractedRegions, chipIds, coordinates, clusterBBoxes, clusterPixels);
}

void ToyProblem::postProcess()
{
  Timer timer("postProcess");
}