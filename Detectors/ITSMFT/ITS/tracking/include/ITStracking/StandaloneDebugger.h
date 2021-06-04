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
/// \fileStandaloneDebugger.h
/// \brief separate TreeStreamerRedirector class to be used with GPU
/// \author matteo.concas@cern.ch

#ifndef O2_ITS_STANDALONE_DEBUGGER_H_
#define O2_ITS_STANDALONE_DEBUGGER_H_

#include <string>
#include <iterator>

// Tracker
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#include "ITStracking/ROframe.h"
#endif

#include "DataFormatsITS/TrackITS.h"
#include "ITStracking/PrimaryVertexContext.h"
#include "ITStracking/Road.h"

namespace o2
{

namespace utils
{
class TreeStreamRedirector;
}

namespace its
{
class Tracklet;
class Line;
class ROframe;
class ClusterLines;

using constants::its::UnusedIndex;

#if !defined(__CUDACC__) && !defined(__HIPCC__)
template <int numClusters = TrackITSExt::MaxClusters>
struct FakeTrackInfo {
 public:
  FakeTrackInfo() = default;
  FakeTrackInfo(PrimaryVertexContext* pvc, const ROframe& event, TrackITSExt& track, bool storeClusters) : isFake{false},
                                                                                                           isAmbiguousId{false},
                                                                                                           track{track},
                                                                                                           mainLabel{UnusedIndex, UnusedIndex, UnusedIndex, false}
  {
    occurrences.clear();
    mcLabels.clear();
    mcLabels.resize(numClusters);
    clusters.clear();
    clusters.resize(numClusters);
    trackingFrameInfos.clear();
    trackingFrameInfos.resize(numClusters);

    for (auto& c : clusStatuses) {
      c = -1;
    }
    for (size_t iCluster{0}; iCluster < numClusters; ++iCluster) {
      int extIndex = track.getClusterIndex(iCluster);
      if (extIndex == -1) {
        continue;
      }
      o2::MCCompLabel mcLabel = event.getClusterLabels(iCluster, extIndex);
      mcLabels[iCluster] = mcLabel;

      bool found = false;
      for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
        std::pair<o2::MCCompLabel, int>& occurrence = occurrences[iOcc];
        if (mcLabel == occurrence.first) {
          ++occurrence.second;
          found = true;
          break;
        }
      }
      if (!found) {
        occurrences.emplace_back(mcLabel, 1);
      }
    }

    if (occurrences.size() > 1) {
      isFake = true;
    }
    std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) {
      return e1.second > e2.second;
    });
    mainLabel = occurrences[0].first;

    for (size_t iOcc{1}; iOcc < occurrences.size(); ++iOcc) {
      if (occurrences[iOcc].second == occurrences[0].second) {
        isAmbiguousId = true;
        break;
      }
    }

    for (size_t iCluster{0}; iCluster < numClusters; ++iCluster) {
      int extIndex = track.getClusterIndex(iCluster);
      if (extIndex == -1) {
        continue; // clusStatuses[iCluster] => -1 // unset
      }
      o2::MCCompLabel lbl = event.getClusterLabels(iCluster, extIndex);
      if (lbl == mainLabel && occurrences[0].second > 1 && !lbl.isNoise()) { // if we have MaxClusters fake clusters -> occurrences[0].second = 1
        clusStatuses[iCluster] = 1;
      } else {
        clusStatuses[iCluster] = 0;
        firstFakeIndex = (firstFakeIndex == -1) ? iCluster : firstFakeIndex;
        ++nFakeClusters;
      }
    }
    if (storeClusters) {
      for (auto iCluster{0}; iCluster < track.getNumberOfClusters(); ++iCluster) {
        const int index = track.getClusterIndex(iCluster);
        if (index != constants::its::UnusedIndex) {
          clusters[iCluster] = pvc->getClusters()[iCluster][index];
          trackingFrameInfos[iCluster] = event.getTrackingFrameInfoOnLayer(iCluster).at(index);
        }
      }
    }
  }

  FakeTrackInfo(PrimaryVertexContext* pvc, const ROframe& event, o2::its::Road& road)
  {
    occurrences.clear();
    mcLabels.clear();
    mcLabels.resize(numClusters);
    clusters.clear();
    clusters.resize(numClusters);
    trackingFrameInfos.clear();
    trackingFrameInfos.resize(numClusters);

    for (auto& c : clusStatuses) {
      c = -1;
    }

    std::vector<int> roadClusters(numClusters, constants::its::UnusedIndex);
    int lastCellLevel = constants::its::UnusedIndex;
    for (int iCell{0}; iCell < numClusters - 2; ++iCell) { // cellsperroad = numClusters -2
      const int cellIndex = road[iCell];
      if (cellIndex == constants::its::UnusedIndex) {
        continue;
      } else {
        roadClusters[iCell] = pvc->getCells()[iCell][cellIndex].getFirstClusterIndex();
        roadClusters[iCell + 1] = pvc->getCells()[iCell][cellIndex].getSecondClusterIndex();
        roadClusters[iCell + 2] = pvc->getCells()[iCell][cellIndex].getThirdClusterIndex();
        assert(roadClusters[iCell] != constants::its::UnusedIndex &&
               roadClusters[iCell + 1] != constants::its::UnusedIndex &&
               roadClusters[iCell + 2] != constants::its::UnusedIndex);
        lastCellLevel = iCell;
      }
    }

    /// From primary vertex context index to event index (== the one used as input of the tracking code)
    for (int iC{0}; iC < roadClusters.size(); iC++) {
      if (roadClusters[iC] != constants::its::UnusedIndex) {
        roadClusters[iC] = pvc->getClusters()[iC][roadClusters[iC]].clusterId;
      }
    }

    for (size_t iCluster{0}; iCluster < numClusters; ++iCluster) {
      int extIndex = roadClusters[iCluster];
      if (extIndex == constants::its::UnusedIndex) {
        continue;
      }
      o2::MCCompLabel mcLabel = event.getClusterLabels(iCluster, extIndex);
      mcLabels[iCluster] = mcLabel;

      bool found = false;
      for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
        std::pair<o2::MCCompLabel, int>& occurrence = occurrences[iOcc];
        if (mcLabel == occurrence.first) {
          ++occurrence.second;
          found = true;
          break;
        }
      }
      if (!found) {
        occurrences.emplace_back(mcLabel, 1);
      }
    }

    if (occurrences.size() > 1) {
      isFake = true;
    }
    std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) {
      return e1.second > e2.second;
    });
    mainLabel = occurrences[0].first;

    for (size_t iOcc{1}; iOcc < occurrences.size(); ++iOcc) {
      if (occurrences[iOcc].second == occurrences[0].second) {
        isAmbiguousId = true;
        break;
      }
    }

    for (size_t iCluster{0}; iCluster < numClusters; ++iCluster) {
      int extIndex = roadClusters[iCluster];
      if (extIndex == constants::its::UnusedIndex) {
        continue;
      }
      o2::MCCompLabel lbl = event.getClusterLabels(iCluster, extIndex);
      if (lbl == mainLabel && occurrences[0].second > 1 && !lbl.isNoise()) { // if we have MaxClusters fake roadClusters -> occurrences[0].second = 1
        clusStatuses[iCluster] = 1;
      } else {
        clusStatuses[iCluster] = 0;
        firstFakeIndex = (firstFakeIndex == -1) ? iCluster : firstFakeIndex;
        ++nFakeClusters;
      }
    }
  }

  // Data
  std::vector<MCCompLabel> mcLabels;
  std::vector<std::pair<MCCompLabel, int>>
    occurrences;
  MCCompLabel mainLabel;
  std::array<int, numClusters> clusStatuses;
  std::vector<o2::its::Cluster> clusters;
  std::vector<o2::its::TrackingFrameInfo> trackingFrameInfos;
  o2::its::TrackITSExt track;

  bool isFake;
  int firstFakeIndex = -1;
  bool isAmbiguousId;
  int nFakeClusters = 0;
  ClassDefNV(FakeTrackInfo, 1);
}; // namespace its
#endif

class StandaloneDebugger
{
 public:
  explicit StandaloneDebugger(const std::string debugTreeFileName = "dbg_ITS.root");
  ~StandaloneDebugger();
  void setDebugTreeFileName(std::string);

  // Monte carlo oracle
  int getEventId(const int firstClusterId, const int secondClusterId, ROframe* frame);

  // Tree part
  const std::string& getDebugTreeFileName() const { return mDebugTreeFileName; }
  void fillCombinatoricsTree(std::array<std::vector<Cluster>, constants::its::LayersNumberVertexer>&,
                             std::vector<Tracklet>,
                             std::vector<Tracklet>,
                             const ROframe*);
  void fillCombinatoricsMCTree(std::vector<Tracklet>, std::vector<Tracklet>);
  void fillTrackletSelectionTree(std::array<std::vector<Cluster>, constants::its::LayersNumberVertexer>&,
                                 std::vector<Tracklet> comb01,
                                 std::vector<Tracklet> comb12,
                                 std::vector<std::array<int, 2>>,
                                 const ROframe*);
  void fillLinesSummaryTree(std::vector<Line>, const ROframe*);
  void fillPairsInfoTree(std::vector<Line>, const ROframe*);
  void fillLineClustersTree(std::vector<ClusterLines> clusters, const ROframe* event);
  void fillXYZHistogramTree(std::array<std::vector<int>, 3>, const std::array<int, 3>);
  void fillVerticesInfoTree(float x, float y, float z, int size, int rId, int eId, float pur);

  // Tracker debug utilities
  void dumpTrackToBranchWithInfo(std::string branchName, int layer, int iteration, o2::its::TrackITSExt track, const ROframe event, PrimaryVertexContext* pvc, const bool dumpClusters = false);
  void dumpTmpTrackToBranchWithInfo(std::string branchName, int layer, int iteration, o2::its::TrackITSExt track, const ROframe event, PrimaryVertexContext* pvc, float pChi2, const bool dumpClusters = false);
  void dumpTrkChi2(float chiFake, float chiTrue);
  void dumpLayerFake(int l, bool hascorrect);
  static int getBinIndex(const float, const int, const float, const float);

  // Smoother debug utilities
  void dumpSmootherChi2(int layer, float original, float smoothed, int startLayer, int trackLength);

 private:
  std::string mDebugTreeFileName = "dbg_ITS.root"; // output filename
  o2::utils::TreeStreamRedirector* mTreeStream;    // observer
};

inline void StandaloneDebugger::setDebugTreeFileName(const std::string name)
{
  if (!name.empty()) {
    mDebugTreeFileName = name;
  }
}

} // namespace its
} // namespace o2

#endif /*O2_ITS_STANDALONE_DEBUGGER_H_*/