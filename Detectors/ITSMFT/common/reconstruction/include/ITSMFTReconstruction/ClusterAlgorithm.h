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

/// \file ClusterAlgorithm.h
/// \brief Definition of a class for the Playne CCL Algorithm for Clustering on the GPU [DOI: 10.1109/TPDS.2018.2799216]
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]
#ifndef ALICEO2_ITS_ClusterAlgorithm_H
#define ALICEO2_ITS_ClusterAlgorithm_H

#include <utility>
#include <vector>
#include <cstring>
#include <memory>
#include <iostream>
#include "ITSMFTReconstruction/PixelData.h"
#include "ITSMFTReconstruction/Clusterer.h"
#include "ITSMFTReconstruction/BoundingBox.h"

namespace o2
{
namespace itsmft
{
using PixelData = o2::itsmft::PixelData;
using BoundingBox = o2::itsmft::BoundingBox;

struct Point {
  int r;
  int c;
};

class ClusterAlgorithm
{
 public:
  ClusterAlgorithm() = default;

  /**
   * \brief Identifies clusters within specified regions of pixel data.
   *
   * This function processes the given regions of pixel data, along with corresponding chip IDs and region coordinates,
   * to identify clusters. The identified clusters are represented by bounding boxes and their corresponding pixel data.
   *
   * \param regions A 3D vector representing the pixel data in different precomputed regions. Each element is an integer representing the pixel value.
   * \param chipIds A vector of integers representing the chip IDs for each region.
   * \param regionCoordinates A vector of (row, column) pairs representing the top-left corner coordinates of each region.
   * \param clusterBBoxes Output vector to be populated with the bounding boxes of identified clusters.
   * \param clusterPixels Output vector to be populated with the pixel data of identified clusters.
   */
  void clusterize(const std::vector<std::vector<std::vector<int>>>& regions, const std::vector<int>& chipIds, const std::vector<std::pair<int, int>>& regionCoordinates,
                  std::vector<BoundingBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_ClusterAlgorithm_H */
