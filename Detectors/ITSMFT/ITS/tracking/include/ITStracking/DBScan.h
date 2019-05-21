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
/// \file DBScan.h
/// \brief
///

// #include "ITStracking/Graph.h"

struct Centroid final {
  Centroid() = default;
  Centroid(int indices[2], float position[3]);
  static float ComputeDistance(const Centroid& c1, const Centroid& c2);
  int mIndices[2];
  float mPosition[3];
};

class DBScan3D
{
};