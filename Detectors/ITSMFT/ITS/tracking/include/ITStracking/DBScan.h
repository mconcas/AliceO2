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

#ifndef O2_ITS_TRACKING_DBSCAN_H_
#define O2_ITS_TRACKING_DBSCAN_H_

#include "ITStracking/Graph.h"

namespace o2
{
namespace its
{

template <typename T>
class DBScan : Graph<T>
{
 public:
  DBScan() = delete;
  explicit DBScan(const size_t nThreads);
  void Init(std::vector<T>&);
  void ClassifyVertices(std::function<unsigned char(std::vector<Edge>&)> classFunction);

 private:
  std::vector<std::vector<unsigned char>> mStates;
  std::function<unsigned char(std::vector<std::vector<Edge>>&)> mClassFunction;
};

template <typename T>
DBScan<T>::DBScan(const size_t nThreads) : Graph<T>(nThreads)
{
}

template <typename T>
void DBScan<T>::Init(std::vector<T>& vertices)
{
  this->Graph<T>::Init(vertices);
}

template <typename T>
void DBScan<T>::ClassifyVertices(std::function<unsigned char(std::vector<Edge>& edges)> classFunction)
{
  mClassFunction = classFunction;
  const size_t size = { this->mVertices.size() };
  mStates.resize(size);

  if (!this->mIsMultiThread) {
    std::cout << "\tClassifying elements" << std::endl;
    for (size_t iVertex{ 0 }; iVertex < size; ++iVertex) {
      mStates[iVertex] = classFunction(this->mEdges[iVertex]);
    }
  } else {
    std::cout << "\tClassifying elements in parallel" << std::endl;
    const size_t stride{ static_cast<size_t>(std::ceil(this->mVertices.size() / static_cast<size_t>(this->mExecutors.size()))) };
    for (size_t iExecutor{ 0 }; iExecutor < this->mExecutors.size(); ++iExecutor) {
      // We cannot pass a template function to std::thread(), using lambda instead
      this->mExecutors[iExecutor] = std::thread(
        [iExecutor, stride, this](const auto& classFunction) {
          for (size_t iVertex{ iExecutor * stride }; iVertex < stride * (iExecutor + 1) && iVertex < this->mVertices.size(); ++iVertex) {
            for (size_t iEdge{ 0 }; iEdge < this->mEdges[iVertex].size(); ++iEdge) {
              classFunction(this->mEdges[iEdge]);
            }
          }
        },
        mClassFunction);
    }
  }
  for (auto&& thread : this->mExecutors) {
    thread.join();
  }
}

struct Centroid final {
  Centroid() = default;
  Centroid(int* indices, float* position);
  void Init();
  static float ComputeDistance(const Centroid& c1, const Centroid& c2);

  int mIndices[2];
  float mPosition[3];
};

} // namespace its
} // namespace o2
#endif