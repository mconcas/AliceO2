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
/// \brief Definitions and implementation of BoundingBox classes used in the ToyProblem
/// \author Nikolaus Draeger [https://cds.cern.ch/record/2879828]

#ifndef ALICEO2_ITS_BOUNDINGBOX_H
#define ALICEO2_ITS_BOUNDINGBOX_H

namespace o2
{
namespace itsmft
{

/**
 * \brief MinimalistBoundingBox is a more abstract version of a bounding box used mainly in the context of CUDA code.
 */
struct MinimalistBoundingBox {
    int min_r;
    int min_c;
    int max_r;
    int max_c;
};

/**
 * \brief BoundingBox represents a box around clusters and includes functionality used in their computation
 *        as well as meta information such as the chip ID.
 */
struct BoundingBox {

    uint16_t chipID = 0xffff;  ///< Chip ID associated with this bounding box
    uint16_t rowMin = 0xffff;  ///< Minimum row index
    uint16_t colMin = 0xffff;  ///< Minimum column index
    uint16_t rowMax = 0;       ///< Maximum row index
    uint16_t colMax = 0;       ///< Maximum column index

    BoundingBox(uint16_t c) : chipID(c) {}

    /**
     * \brief Check if a point (row, col) is inside the bounding box.
     * \param row Row index of the point.
     * \param col Column index of the point.
     * \return True if the point is inside the bounding box, false otherwise.
     */
    bool isInside(uint16_t row, uint16_t col) const {
      return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
    }

    /**
     * \brief Get the span of rows covered by the bounding box.
     * \return Number of rows covered by the bounding box.
     */
    auto rowSpan() const { return rowMax - rowMin + 1; }

    /**
     * \brief Get the span of columns covered by the bounding box.
     * \return Number of columns covered by the bounding box.
     */
    auto colSpan() const { return colMax - colMin + 1; }

    /**
     * \brief Clear the bounding box, resetting it to its initial state.
     */
    void clear()
    {
      rowMin = colMin = 0xffff;
      rowMax = colMax = 0;
    }

    /**
     * \brief Adjust the bounding box such that the given point (row, col) is within bounds.
     * \param row Row index of the point.
     * \param col Column index of the point.
     */
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
