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
/// \brief concrete strategy class for clustering using Playne GPU algorithm [DOI: 10.1109/TPDS.2018.2799216]
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
using BoundingBox2 = o2::itsmft::BoundingBox2;

struct Point {
  int r;
  int c;
};

class ClusterAlgorithm
{
 public:
  ClusterAlgorithm() = default;
  void clusterize(const std::vector<std::vector<std::vector<int>>>& regions, const std::vector<int>& chipIds, const std::vector<std::pair<int,int>>& regionCoordinates,
                  std::vector<BoundingBox2>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels);
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_ClusterAlgorithm_H */
