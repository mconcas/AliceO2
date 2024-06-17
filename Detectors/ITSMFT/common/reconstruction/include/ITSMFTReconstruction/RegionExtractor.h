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

/// \file RegionExtractor.h
/// \brief Abstract strategy class for region extraction
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]
#ifndef ALICEO2_ITS_REGIONEXTRACTOR_H
#define ALICEO2_ITS_REGIONEXTRACTOR_H

#include <utility>
#include <vector>
#include <cstring>
#include <memory>
#include <iostream>
#include "ITSMFTReconstruction/PixelData.h"

namespace o2
{
namespace itsmft
{

/**
 * \brief Wrapper struct to hold extracted regions and their coordinates.
 */
struct ExtractedRegionsWrapper {
  std::vector<std::vector<std::vector<int>>> regions;
  std::vector<std::pair<int, int>> coordinates;
};

class RegionExtractor
{
 public:
  using ChipPixelData = o2::itsmft::ChipPixelData;

  virtual ~RegionExtractor() = default;

  /**
   * \brief Pure virtual function to preprocess the chip data and extract subregions.
   *
   * This function needs to be implemented by derived classes to preprocess the given chip data and extract regions
   * based on the specified maximum distances in x and y directions.
   * This step is used to reduce load on the CUDA kernel.
   *
   * \param chipData Pointer to the chip data to be processed.
   * \param maxdist_x Meta information for preprocessing: Maximum allowed distance in the x direction for region expansion.
   * \param maxdist_y Meta information for preprocessing: Maximum allowed distance in the y direction for region expansion.
   * \return ExtractedRegionsWrapper containing the extracted regions.
   */
  virtual ExtractedRegionsWrapper preprocess(const ChipPixelData* chipData, const int maxdist_x, const int maxdist_y) = 0;
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_REGIONEXTRACTOR_H */