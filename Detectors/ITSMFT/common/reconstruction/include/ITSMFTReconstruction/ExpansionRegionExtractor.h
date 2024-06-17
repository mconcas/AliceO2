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

/// \file ExpansionRegionExtractor.h
/// \brief concrete strategy class for region extraction using expansion method
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]
#ifndef ALICEO2_ITS_EXPANSIONREGIONEXTRACTOR_H
#define ALICEO2_ITS_EXPANSIONREGIONEXTRACTOR_H

#include <utility>
#include <vector>
#include <cstring>
#include <memory>
#include <iostream>
#include <set>
#include <algorithm>
#include "ITSMFTReconstruction/RegionExtractor.h"

namespace o2
{
namespace itsmft
{

class ExpansionRegionExtractor : public RegionExtractor
{
  using PixelData = o2::itsmft::PixelData;

 public:
  /**
   * \brief Preprocess the chip data by extracting subregions using the expansion method.
   *
   * This function preprocesses the given chip data to identify subregions containing pixel data using an expansion method.
   * It scans the maximum allowed distances in x and y directions to form a set of subregions.
   * This step is used to reduce load on the CUDA kernel.
   *
   * \param chipData Pointer to the chip data to be processed.
   * \param maxdist_x Maximum allowed distance in the x direction for region expansion.
   * \param maxdist_y Maximum allowed distance in the y direction for region expansion.
   * \return ExtractedRegionsWrapper containing the extracted regions.
   */
  ExtractedRegionsWrapper preprocess(const ChipPixelData* chipData, const int maxdist_x, const int maxdist_y) override;

 private:
  struct Region {
    int row;    ///< Row-coordinate of the top-left corner
    int col;    ///< Column-coordinate of the top-left corner
    int width;  ///< Width of the region
    int height; ///< Height of the region
  };

  /**
   * \brief Expand a region based on the maximum allowed distances.
   *
   * This function expands the given region based on the maximum allowed distances in x and y directions.
   * It checks whether there is a pixel within the maximum distances around the current region and swallows
   * any pixel within that range.
   *
   * \param fullRegion 2D vector representing the full region.
   * \param regionInfo Struct containing the subregion information.
   * \param maxdist_x Maximum allowed distance scanned in the x direction for expansion.
   * \param maxdist_y Maximum allowed distance scanned in the y direction for expansion.
   * \param pixelSet Set of PixelData to be updated during the expansion.
   * \param pixelDataPointers 2D vector of pointers to PixelData.
   * \return True if the region was successfully expanded, false otherwise.
   */
  bool expandRegion(std::vector<std::vector<int>>& fullRegion,
                    Region& regionInfo,
                    const int maxdist_x,
                    const int maxdist_y,
                    std::set<PixelData>& pixelSet,
                    std::vector<std::vector<PixelData*>>& pixelDataPointers);

  /**
   * \brief Convert sparse pixel data to a grid representation.
   *
   * \param pixels Vector of sparse PixelData representation to be converted.
   * \return 2D vector representing the pixel data in grid form.
   */
  std::vector<std::vector<int>> convertSparsePixelsToGrid(const std::vector<PixelData> pixels);
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_EXPANSIONREGIONEXTRACTOR_H */