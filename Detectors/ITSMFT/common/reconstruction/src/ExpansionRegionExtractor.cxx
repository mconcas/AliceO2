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

/// \file ExpansionRegionExtractor.cxx
/// \brief Implementation of the region extraction using expanding regions.

#include "ITSMFTReconstruction/Clusterer.h"
#include "ITSMFTReconstruction/ExpansionRegionExtractor.h"

using namespace o2::itsmft;

ExtractedRegionsWrapper ExpansionRegionExtractor::preprocess(const ChipPixelData* chipData, const int maxdist_x, const int maxdist_y)
{
  std::vector<std::vector<std::vector<int>>> extractedRegions;
  std::vector<PixelData> pixelData = chipData->getData();

  // could be replaced by unordered_set for better runtime complexity
  // requires hash function for PixelData, however
  std::set<PixelData> pixelSet(pixelData.begin(), pixelData.end());

  std::vector<std::vector<int>> fullRegion = convertSparsePixelsToGrid(pixelData);
  
  int numRows = fullRegion.size();
  int numCols = 0;
  if (!fullRegion.empty()) {
    numCols = fullRegion[0].size();
  }

  std::vector<std::vector<PixelData*>> pixelDataPointers(numRows, std::vector<PixelData*>(numCols, nullptr));
  for (const PixelData& pixel : pixelSet) {
    const PixelData* ptr = &pixel;
    pixelDataPointers[pixel.getRow()][pixel.getCol()] = const_cast<PixelData*>(ptr);
  }

  std::vector<std::pair<int, int>> regionCoordinates;

  while (!pixelSet.empty()) {
    PixelData nextPixel = *pixelSet.begin();
    int nextRow = nextPixel.getRow();
    int nextCol = nextPixel.getCol();
    Region regionInfo = {nextRow, nextCol, 1, 1};

    while (expandRegion(fullRegion, regionInfo, maxdist_x, maxdist_y, pixelSet, pixelDataPointers)) { }

    std::vector<std::vector<int>> region(regionInfo.height, std::vector<int>(regionInfo.width, 0));

    for (int i = 0; i < regionInfo.height; ++i) {
      for (int j = 0; j < regionInfo.width; ++j) {
        region[i][j] = fullRegion[i + regionInfo.row][j + regionInfo.col];
      }
    }

    for (int row = regionInfo.row; row < regionInfo.row + regionInfo.height; ++row) {
      for (int col = regionInfo.col; col < regionInfo.col + regionInfo.width; ++col) {

        if (row >= fullRegion.size() || col >= fullRegion[row].size()) {
          std::cout << "ERROR: Trying to access out of bounds index in fullRegion" << std::endl;
        }

        if (fullRegion[row][col] != 1) continue;

        fullRegion[row][col] = 0;

        if (row >= pixelDataPointers.size() || col >= pixelDataPointers[row].size()) {
          std::cout << "ERROR: Trying to access out of bounds index in pixelDataPointers" << std::endl;
        }

        PixelData* pixelData = pixelDataPointers[row][col];
        if (pixelData) {
          if (pixelSet.find(*pixelData) == pixelSet.end()) {
            std::cout << "ERROR: Trying to erase an object from pixelSet that doesn't exist" << std::endl;
          }

          pixelSet.erase(*pixelData);
          pixelDataPointers[row][col] = nullptr;
        }
      }
    }
    regionCoordinates.push_back(std::make_pair(regionInfo.row, regionInfo.col));
    extractedRegions.push_back(region);
  }
  ExtractedRegionsWrapper wrappedRegions = { extractedRegions, regionCoordinates };
  return wrappedRegions;
}

bool ExpansionRegionExtractor::expandRegion(std::vector<std::vector<int>>& fullRegion,
                                            Region& regionInfo,
                                            const int maxdist_x,
                                            const int maxdist_y,
                                            std::set<PixelData>& pixelSet,
                                            std::vector<std::vector<PixelData*>>& pixelDataPointers)
{
  int fullRegionHeight = fullRegion.size();
  int fullRegionWidth = fullRegion[0].size();

  int startRow = std::max(regionInfo.row - maxdist_x, 0);
  int startCol = std::max(regionInfo.col - maxdist_y, 0);
  int endRow = std::min(regionInfo.row + regionInfo.height + maxdist_x, fullRegionHeight);
  int endCol = std::min(regionInfo.col + regionInfo.width + maxdist_y, fullRegionWidth);

  bool regionExpanded = false;

  for (int row = startRow; row < endRow; ++row) {
    for (int col = startCol; col < endCol; ++col) {
      if (row >= regionInfo.row && row < regionInfo.row + regionInfo.height &&
          col >= regionInfo.col && col < regionInfo.col + regionInfo.width) {
        continue;
      }

      if (fullRegion[row][col] == 1) {
        if (col < regionInfo.col) {
          regionInfo.col = col;
          regionInfo.width += (regionInfo.col - col + 1);
        } else if (col >= regionInfo.col + regionInfo.width) {
          regionInfo.width += (col - regionInfo.col - regionInfo.width + 1);
        }

        if (row < regionInfo.row) {
          regionInfo.row = row;
          regionInfo.height += (regionInfo.row - row + 1);
        } else if (row >= regionInfo.row + regionInfo.height) {
          regionInfo.height += (row - regionInfo.row - regionInfo.height + 1);
        }
        regionExpanded = true;
      }
    }
  }
  return regionExpanded;
}

std::vector<std::vector<int>> ExpansionRegionExtractor::convertSparsePixelsToGrid(const std::vector<PixelData> pixels)
{
  uint16_t maxRow = 0, maxCol = 0;
  for (const auto& pixel : pixels) {
    maxRow = std::max(maxRow, pixel.getRow());
    maxCol = std::max(maxCol, pixel.getCol());
  }

  std::vector<std::vector<int>> grid(maxRow + 1, std::vector<int>(maxCol + 1, 0));

  for (const auto& pixel : pixels) {
    grid[pixel.getRow()][pixel.getCol()] = 1;
  }

  return grid;
}