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
#include <iostream>
#include <utility>
#include <vector>



namespace o2
{
namespace its
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
  Graph(const std::vector<T>& vertices, std::function<bool(const T&, const T&)> linkFunction);

 private:
  std::vector<std::pair<int, int>> mVerticesEdgeData; //
  std::vector<Edge> mEdges;                           // Adjacency list
};

template <typename T>
Graph<T>::Graph(const std::vector<T>& vertices, std::function<bool(const T& v1, const T& v2)> linkFunction)
{
  mVerticesEdgeData.resize(vertices.size());
  int tot_nedges{ 0 };
  for (size_t iVertex1{ 0 }; iVertex1 < vertices.size(); ++iVertex1) {
    int nedges{ 0 };
    for (size_t iVertex2{ 0 }; iVertex2 < vertices.size(); ++iVertex2) {
      if (iVertex1 != iVertex2 && linkFunction(vertices[iVertex1], vertices[iVertex2])) {
        mEdges.emplace_back(std::make_pair(iVertex1, iVertex2));
        ++nedges;
      }
    }
    tot_nedges += nedges;
    mVerticesEdgeData[iVertex1] = std::make_pair(nedges, tot_nedges);
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

} // namespace its
} // namespace o2
#endif