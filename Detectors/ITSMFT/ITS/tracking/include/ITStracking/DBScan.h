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

#include "ITStracking/Graph.h"

namespace o2
{
namespace its
{

struct Centroid final {
  Centroid() = default;
  Centroid(int* indices, float* position);
  static float ComputeDistance(const Centroid& c1, const Centroid& c2);

  int mIndices[2];
  float mPosition[3];
};

template <typename T>
class DBScan : Graph<T>
{
 public:
  DBScan() = delete;
  explicit DBScan(const std::vector<T>& vertices, const size_t nThreads);
  void ClassifyVertices(std::function<unsigned char(std::vector<std::vector<Edge>>&)> classFunction);

 private:
  std::vector<std::vector<unsigned char>> mStates;
  std::function<unsigned char(std::vector<std::vector<Edge>>&)> mClassFunction;
};

template <typename T>
DBScan<T>::DBScan(const std::vector<T>& vertices, const size_t nThreads) : Graph<T>(vertices, nThreads)
{
}

template <typename T>
void DBScan<T>::ClassifyVertices(std::function<unsigned char(std::vector<std::vector<Edge>>& edges)> classFunction)
{
  mClassFunction = classFunction;
  const size_t size = { this->mVertices.size() };
  mStates.resize(size);

  if (!this->mIsMultiThread) {
    std::cout << "Single thread implementation" << std::endl;
    for (size_t iVertex{ 0 }; iVertex < size; ++iVertex) {
      mStates[iVertex] = classFunction(this->mEdges);
    }
  } // else {
  //   std::cout << "Multithread implementation" << std::endl;
  //   mNThreads = std::min(static_cast<const size_t>(std::thread::hardware_concurrency()), mNThreads);
  //   mExecutors.resize(mNThreads);
  //   const size_t stride{ static_cast<size_t>(std::ceil(mVertices.size() / static_cast<size_t>(mExecutors.size()))) };
  //   for (size_t iExecutor{ 0 }; iExecutor < mExecutors.size(); ++iExecutor) {
  //     // We cannot pass a template function to std::thread, using lambda instead
  //     mExecutors[iExecutor] = std::thread(
  //       [iExecutor, stride, this](const auto& linkFunction) {
  //         for (size_t iVertex1{ iExecutor * stride }; iVertex1 < stride * (iExecutor + 1) && iVertex1 < mVertices.size(); ++iVertex1) {
  //           for (size_t iVertex2{ 0 }; iVertex2 < mVertices.size(); ++iVertex2) {
  //             if (iVertex1 != iVertex2 && linkFunction(mVertices[iVertex1], mVertices[iVertex2])) {
  //               mEdges[iVertex1].emplace_back(iVertex1, iVertex2);
  //             }
  //           }
  //         }
  //       },
  //       mLinkFunction);
  //   }
  // }
  // for (auto&& thread : mExecutors) {
  //   thread.join();
  // }
}

} // namespace its
} // namespace o2