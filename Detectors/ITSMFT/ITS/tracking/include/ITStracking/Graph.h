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
/// \file Graph.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_ALGORITHMS_H_
#define TRACKINGITSU_INCLUDE_ALGORITHMS_H_

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

namespace o2
{
namespace its
{

typedef std::pair<int, int> Edge;

template <typename T>
class Graph
{
 public:
  Graph() = delete;
  explicit Graph(const size_t nThreads = 1);
  void init(std::vector<T>&);
  std::vector<int> getCluster(const int);
  void computeEdges(std::function<bool(const T& v1, const T& v2)>);
  std::vector<std::vector<Edge>> getEdges() const { return mEdges; }
  char isMultiThreading() const { return mIsMultiThread; }

  std::vector<T>* mVertices = nullptr; // Observer pointer
  std::vector<std::thread> mExecutors; // Difficult to pass

 private:
  void findVertexEdges(std::vector<Edge>& localEdges, const T& vertex, const size_t vId, const size_t size);

  // Multithread block
  size_t mNThreads;
  char mIsMultiThread;

  // Common data members
  std::function<bool(const T&, const T&)> mLinkFunction;
  std::vector<std::vector<Edge>> mEdges;
  std::vector<char> mVisited;
};

template <typename T>
Graph<T>::Graph(const size_t nThreads) : mNThreads{ nThreads }
{
  mIsMultiThread = nThreads > 1 ? true : false;
}

template <typename T>
void Graph<T>::init(std::vector<T>& vertices)
{

  // Graph initialization
  mVertices = &vertices;
  if (mIsMultiThread) {
    mNThreads = std::min(static_cast<const size_t>(std::thread::hardware_concurrency()), mNThreads);
    mExecutors.resize(mNThreads);
  }

  mEdges.resize(vertices.size());
  mVisited.resize(vertices.size(), false);
}

template <typename T>
void Graph<T>::computeEdges(std::function<bool(const T& v1, const T& v2)> linkFunction)
{
  mLinkFunction = linkFunction;
  int tot_nedges = 0;
  const size_t size = { mVertices->size() };
  if (!mIsMultiThread) {
    std::cout << "\tSingle thread implementation" << std::endl;
    for (size_t iVertex{ 0 }; iVertex < size; ++iVertex) {
      findVertexEdges(mEdges[iVertex], (*mVertices)[iVertex], iVertex, size);
      tot_nedges += static_cast<int>(mEdges[iVertex].size());
    }
  } else {
    std::cout << "\tMultithread implementation" << std::endl;
    mNThreads = std::min(static_cast<const size_t>(std::thread::hardware_concurrency()), mNThreads);
    mExecutors.resize(mNThreads);
    const size_t stride{ static_cast<size_t>(std::ceil(mVertices->size() / static_cast<size_t>(mExecutors.size()))) };
    for (size_t iExecutor{ 0 }; iExecutor < mExecutors.size(); ++iExecutor) {
      // We cannot pass a template function to std::thread(), using lambda instead
      mExecutors[iExecutor] = std::thread(
        [iExecutor, stride, this](const auto& linkFunction) {
          for (size_t iVertex1{ iExecutor * stride }; iVertex1 < stride * (iExecutor + 1) && iVertex1 < mVertices->size(); ++iVertex1) {
            for (size_t iVertex2{ 0 }; iVertex2 < mVertices->size(); ++iVertex2) {
              if (iVertex1 != iVertex2 && linkFunction((*mVertices)[iVertex1], (*mVertices)[iVertex2])) {
                mEdges[iVertex1].emplace_back(iVertex1, iVertex2);
              }
            }
          }
        },
        mLinkFunction);
    }
  }
  for (auto&& thread : mExecutors) {
    thread.join();
  }
}
template <typename T>
void Graph<T>::findVertexEdges(std::vector<Edge>& localEdges, const T& vertex, const size_t vId, const size_t size)
{
  for (size_t iVertex2{ 0 }; iVertex2 < size; ++iVertex2) {
    if (vId != iVertex2 && mLinkFunction(vertex, (*mVertices)[iVertex2])) {
      localEdges.emplace_back(std::make_pair(vId, iVertex2));
    }
  }
}

template <typename T>
std::vector<int> Graph<T>::getCluster(const int vertexId)
{
  std::queue<int> idQueue;
  std::vector<int> indices;
  idQueue.emplace(vertexId);
  mVisited[vertexId] = true;
  // indices.emplace_back(vertexId);
  if (!mIsMultiThread) {
    // std::cout << "\tSingle thread clustering" << std::endl;
    // Consume the queue
    while (!idQueue.empty()) {
      const int id = idQueue.front();
      idQueue.pop();
      
      for (Edge edge : mEdges[id]) {
        if (!mVisited[edge.second]) {
          idQueue.emplace(edge.second);
          indices.emplace_back(edge.second);
          mVisited[edge.second] = true;
        }
      }
    }
  }
  std::fill(mVisited.begin(), mVisited.end(), false);
  return indices;
}

} // namespace its
} // namespace o2

#endif