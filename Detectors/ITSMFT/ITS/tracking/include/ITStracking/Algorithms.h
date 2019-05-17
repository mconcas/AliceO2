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
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include "ITStracking/Definitions.h"

namespace o2
{
namespace ITS
{

typedef std::pair<int, int> Edge;

struct Centroid final {
  Centroid() = default;
  Centroid(int indices[2], float position[3]);
  static float ComputeDistance(const Centroid& c1, const Centroid& c2);
  int mIndices[2];
  float mPosition[3];
};

template <typename T>
class Graph
{
 public:
  Graph() = delete;
  template <typename... Args>
  Graph(const std::vector<T>& vertices, bool (*linkFunction)(T&, T&, Args...), Args&&... args);

 private:
  std::vector<GraphVertex<Centroid>> mGraphVertices;
  std::vector<std::pair<int, int>> mVerticesEdgeData; //
  std::vector<Edge> mEdges;                           // Adjacency list
};

template <typename T>
template <typename... Args>
Graph<T>::Graph(const std::vector<T>& vertices, bool (*linkFunction)(T& v1, T& v2, Args...), Args&&... args)
{
  int tot_nedges{ 0 };
  for (int iVertex1{ 0 }; iVertex1 < vertices.size() - 1; ++iVertex1) {
    int nedges{ 0 };
    for (int iVertex2{ iVertex1 }; iVertex2 < vertices.size(); ++iVertex2) {
      if ((this->*linkFunction)(vertices[iVertex1], vertices[iVertex2], std::forward<Args>(args)...)) {
        mEdges.emplace_back(std::make_pair(iVertex1, iVertex2));
        ++nedges;
      }
    }
    tot_nedges += nedges;
    mVerticesEdgeData[i] = std::make_pair(nedges, totedges);
  }
}

// template <typename T>
// class GraphVertex
// {
//  public:
//   GraphVertex() = delete;
//   T GetPayload() const { return mPayload; }
//   template <typename... Args>
//   GraphVertex(Args&&... args) : mPayload{ std::forward<Args>(args)... }
//   {
//     mVisited = false;
//   }

//  private:
//   T mPayload;
//   char mVisited;
//   // three possible statuses: core, border, noise
// };
class DBScan3D
{
};

} // namespace ITS
} // namespace o2
#endif