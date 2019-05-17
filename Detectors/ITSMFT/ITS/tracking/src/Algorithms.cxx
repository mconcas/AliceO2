// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Algorithms.cxx
/// \brief
///

#include <iterator>
#include "ITStracking/Algorithms.h"

namespace o2
{
namespace ITS
{

Centroid::Centroid(int indices[2], float position[3])
{
  for (int i{ 0 }; i < 2; ++i) {
    mIndices[i] = indices[i];
  }
  for (int i{ 0 }; i < 3; ++i) {
    mPosition[i] = position[i];
  }
}

float Centroid::ComputeDistance(const Centroid& c1, const Centroid& c2)
{
  return MATH_SQRT((c1.mPosition[0] - c2.mPosition[0]) * (c1.mPosition[0] - c2.mPosition[0]) +
                   (c1.mPosition[1] - c2.mPosition[1]) * (c1.mPosition[1] - c2.mPosition[1]) +
                   (c1.mPosition[2] - c2.mPosition[2]) * (c1.mPosition[2] - c2.mPosition[2]));
}

} // namespace ITS
} // namespace o2