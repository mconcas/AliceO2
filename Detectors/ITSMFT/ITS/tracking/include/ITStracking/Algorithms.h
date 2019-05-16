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
/// \file Algorithms.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_ALGORITHMS_H_
#define TRACKINGITSU_INCLUDE_ALGORITHMS_H_

#include <array>
#include <utility>
#include <vector>

namespace o2
{
namespace ITS
{

struct Centroid {
    Centroid() = default;
    GPU_HOST_DEVICE Centroid(int indices[2], float position[3]);

    int mIndices[2];
    float mPosition[3];
}

inline GPU_HOST_DEVICE Centroid::Centroid(int indices[2], float position[3]) {
  for (int i{0}; i < 2; ++i) {
    mIndices[i] = indices[i];
  }
  for (int i{0}; i < 3; ++i) {
      mPosition[i] = position[i];
  }
}

class CentroidGraph {
  public:
  CentroidGraph() = delete;
  GPU_HOST CentroidGraph(std::vector<Centroids>& centroids);
  GPU_HOST_DEVICE Centroid(Centroids, const int size);
  private:
  std::vector<Centroid> mCentroids;
  typedef std::pair<int, int> Edge;
  // LookupTable to browse edges, Indices are vertices
  std::vector<std::pair<int, int>> mCompactAdjList; 
  std::vector<Edge> mEdges;
};

class DBScan3D {
  
};

}
}
#endif