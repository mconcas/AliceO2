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

/// \file BoundingBox.h
/// \brief Todo
#ifndef ALICEO2_ITS_BOUNDINGBOX_H
#define ALICEO2_ITS_BOUNDINGBOX_H

namespace o2
{
namespace itsmft
{

struct BoundingBox {
  int min_r;
  int min_c;
  int max_r;
  int max_c;
};

struct BoundingBox2 {
    uint16_t chipID = 0xffff;
    uint16_t rowMin = 0xffff;
    uint16_t colMin = 0xffff;
    uint16_t rowMax = 0;
    uint16_t colMax = 0;
    BoundingBox2(uint16_t c) : chipID(c) {}
    bool isInside(uint16_t row, uint16_t col) const { return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax; }
    auto rowSpan() const { return rowMax - rowMin + 1; }
    auto colSpan() const { return colMax - colMin + 1; }
    void clear()
    {
      rowMin = colMin = 0xffff;
      rowMax = colMax = 0;
    }
    void adjust(uint16_t row, uint16_t col)
    {
      if (row < rowMin) {
        rowMin = row;
      }
      if (row > rowMax) {
        rowMax = row;
      }
      if (col < colMin) {
        colMin = col;
      }
      if (col > colMax) {
        colMax = col;
      }
    }
  };

} // namespace itsmft
} // namespace o2
#endif /* ALICEO2_ITS_BOUNDINGBOX_H */
