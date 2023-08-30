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
    ExtractedRegionsWrapper preprocess(const ChipPixelData* chipData, const int maxdist_x, const int maxdist_y) override;

private:

    struct Region
    {
        int row;  // row-coordinate of the top-left corner
        int col;  // col-coordinate of the top-left corner
        int width;
        int height;
    };

    bool expandRegion(std::vector<std::vector<int>>& fullRegion,
                                                Region& regionInfo, 
                                                const int maxdist_x,
                                                const int maxdist_y,
                                                std::set<PixelData>& pixelSet, 
                                                std::vector<std::vector<PixelData*>>& pixelDataPointers);

    std::vector<std::vector<int>> convertSparsePixelsToGrid(const std::vector<PixelData> pixels);
};

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_EXPANSIONREGIONEXTRACTOR_H */