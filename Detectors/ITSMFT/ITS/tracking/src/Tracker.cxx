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
/// \file Tracker.cxx
/// \brief
///

#include "ITStracking/Tracker.h"

#include "ITStracking/Cell.h"
#include "ITStracking/Constants.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITStracking/Smoother.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/TrackerTraits.h"
#include "ITStracking/TrackerTraitsCPU.h"
#include "ITStracking/TrackingConfigParam.h"
#include "DetectorsBase/GeometryManager.h"

#include "ReconstructionDataFormats/Track.h"
#include <cassert>
#include <iostream>
#include <utility>
#include <dlfcn.h>
#include <cstdlib>
#include <string>

#include <random>

namespace o2
{
namespace its
{

using constants::its::UnusedIndex;
const MCCompLabel unusedMCLabel = {UnusedIndex, UnusedIndex,
                                   UnusedIndex, false};

Tracker::Tracker(o2::its::TrackerTraits* traits)
{
  /// Initialise standard configuration with 1 iteration
  mTrkParams.resize(1);
  mMemParams.resize(1);
  mTraits = traits;
  mPrimaryVertexContext = mTraits->getPrimaryVertexContext();
#ifdef CA_STANDALONE_DEBUGGER
  mDebugger = new StandaloneDebugger("dbg_ITSTrackerCPU.root");
#endif
}

#if defined(CA_DEBUG) || defined(CA_STANDALONE_DEBUGGER)
Tracker::~Tracker()
{
  delete mDebugger;
}
#else
Tracker::~Tracker() = default;
#endif

void Tracker::clustersToTracks(const ROframe& event, std::ostream& timeBenchmarkOutputStream)
{
  const int verticesNum = event.getPrimaryVerticesNum();
  mTracks.clear();
  mTrackLabels.clear();
  mDebugRoad.clear();

  for (int iVertex = 0; iVertex < verticesNum; ++iVertex) {

    float total{0.f};

    for (int iteration = 0; iteration < mTrkParams.size(); ++iteration) {

      int numCls = 0;
      for (unsigned int iLayer{0}; iLayer < mTrkParams[iteration].NLayers; ++iLayer) {
        numCls += event.getClusters()[iLayer].size();
      }
      if (numCls < mTrkParams[iteration].MinTrackLength) {
        continue;
      }

      mTraits->UpdateTrackingParameters(mTrkParams[iteration]);
      /// Ugly hack -> Unifiy float3 definition in CPU and CUDA/HIP code
      int pass = iteration + iVertex; /// Do not reinitialise the context if we analyse pile-up events
      std::array<float, 3> pV = {event.getPrimaryVertex(iVertex).x, event.getPrimaryVertex(iVertex).y, event.getPrimaryVertex(iVertex).z};
      total += evaluateTask(&Tracker::initialisePrimaryVertexContext, "Context initialisation",
                            timeBenchmarkOutputStream, mMemParams[iteration], mTrkParams[iteration], event.getClusters(), pV, pass);
      total += evaluateTask(&Tracker::computeTracklets, "Tracklet finding", timeBenchmarkOutputStream);
      total += evaluateTask(&Tracker::computeCells, "Cell finding", timeBenchmarkOutputStream);
      total += evaluateTask(&Tracker::findCellsNeighbours, "Neighbour finding", timeBenchmarkOutputStream, iteration);
      total += evaluateTask(&Tracker::findRoads, "Road finding", timeBenchmarkOutputStream, iteration);
      total += evaluateTask(&Tracker::findTracks, "Track finding", timeBenchmarkOutputStream, event);
    }

    if (constants::DoTimeBenchmarks && fair::Logger::Logging(fair::Severity::info)) {
      timeBenchmarkOutputStream << std::setw(2) << " - "
                                << "Vertex processing completed in: " << total << "ms" << std::endl;
    }
  }
  if (verticesNum > 0) {
    float smoothElapsed = evaluateTask(&Tracker::smoothTracks, "Track smoothing", timeBenchmarkOutputStream, event);
  }
  if (event.hasMCinformation()) {
    computeTracksMClabels(event);
  } else {
    rectifyClusterIndices(event);
  }
}

void Tracker::computeTracklets()
{
  mTraits->computeLayerTracklets();
}

void Tracker::computeCells()
{
  mTraits->computeLayerCells();
}

void Tracker::findCellsNeighbours(int& iteration)
{
  for (int iLayer{0}; iLayer < mTrkParams[iteration].CellsPerRoad() - 1; ++iLayer) {

    if (mPrimaryVertexContext->getCells()[iLayer + 1].empty() ||
        mPrimaryVertexContext->getCellsLookupTable()[iLayer].empty()) {
      continue;
    }

    int layerCellsNum{static_cast<int>(mPrimaryVertexContext->getCells()[iLayer].size())};
    const int nextLayerCellsNum{static_cast<int>(mPrimaryVertexContext->getCells()[iLayer + 1].size())};
    mPrimaryVertexContext->getCellsNeighbours()[iLayer].resize(nextLayerCellsNum);

    for (int iCell{0}; iCell < layerCellsNum; ++iCell) {

      const Cell& currentCell{mPrimaryVertexContext->getCells()[iLayer][iCell]};
      const int nextLayerTrackletIndex{currentCell.getSecondTrackletIndex()};
      const int nextLayerFirstCellIndex{mPrimaryVertexContext->getCellsLookupTable()[iLayer][nextLayerTrackletIndex]};
      if (nextLayerFirstCellIndex != constants::its::UnusedIndex &&
          mPrimaryVertexContext->getCells()[iLayer + 1][nextLayerFirstCellIndex].getFirstTrackletIndex() ==
            nextLayerTrackletIndex) {

        for (int iNextLayerCell{nextLayerFirstCellIndex}; iNextLayerCell < nextLayerCellsNum; ++iNextLayerCell) {

          Cell& nextCell{mPrimaryVertexContext->getCells()[iLayer + 1][iNextLayerCell]};
          if (nextCell.getFirstTrackletIndex() != nextLayerTrackletIndex) {
            break;
          }

          const float3 currentCellNormalVector{currentCell.getNormalVectorCoordinates()};
          const float3 nextCellNormalVector{nextCell.getNormalVectorCoordinates()};
          const float3 normalVectorsDeltaVector{currentCellNormalVector.x - nextCellNormalVector.x,
                                                currentCellNormalVector.y - nextCellNormalVector.y,
                                                currentCellNormalVector.z - nextCellNormalVector.z};

          const float deltaNormalVectorsModulus{(normalVectorsDeltaVector.x * normalVectorsDeltaVector.x) +
                                                (normalVectorsDeltaVector.y * normalVectorsDeltaVector.y) +
                                                (normalVectorsDeltaVector.z * normalVectorsDeltaVector.z)};
          const float deltaCurvature{std::abs(currentCell.getCurvature() - nextCell.getCurvature())};

          if (deltaNormalVectorsModulus < mTrkParams[iteration].NeighbourMaxDeltaN[iLayer] &&
              deltaCurvature < mTrkParams[iteration].NeighbourMaxDeltaCurvature[iLayer]) {

            mPrimaryVertexContext->getCellsNeighbours()[iLayer][iNextLayerCell].push_back(iCell);

            const int currentCellLevel{currentCell.getLevel()};

            if (currentCellLevel >= nextCell.getLevel()) {

              nextCell.setLevel(currentCellLevel + 1);
            }
          }
        }
      }
    }
  }
}

void Tracker::findRoads(int& iteration)
{
  for (int iLevel{mTrkParams[iteration].CellsPerRoad()}; iLevel >= mTrkParams[iteration].CellMinimumLevel(); --iLevel) {
    CA_DEBUGGER(int nRoads = -mPrimaryVertexContext->getRoads().size());
    const int minimumLevel{iLevel - 1};

    for (int iLayer{mTrkParams[iteration].CellsPerRoad() - 1}; iLayer >= minimumLevel; --iLayer) {

      const int levelCellsNum{static_cast<int>(mPrimaryVertexContext->getCells()[iLayer].size())};

      for (int iCell{0}; iCell < levelCellsNum; ++iCell) {

        Cell& currentCell{mPrimaryVertexContext->getCells()[iLayer][iCell]};

        if (currentCell.getLevel() != iLevel) {
          continue;
        }

        mPrimaryVertexContext->getRoads().emplace_back(iLayer, iCell);

        /// For 3 clusters roads (useful for cascades and hypertriton) we just store the single cell
        /// and we do not do the candidate tree traversal
        if (iLevel == 1) {
          continue;
        }

        const int cellNeighboursNum{static_cast<int>(
          mPrimaryVertexContext->getCellsNeighbours()[iLayer - 1][iCell].size())};
        bool isFirstValidNeighbour = true;

        for (int iNeighbourCell{0}; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

          const int neighbourCellId = mPrimaryVertexContext->getCellsNeighbours()[iLayer - 1][iCell][iNeighbourCell];
          const Cell& neighbourCell = mPrimaryVertexContext->getCells()[iLayer - 1][neighbourCellId];

          if (iLevel - 1 != neighbourCell.getLevel()) {
            continue;
          }

          if (isFirstValidNeighbour) {

            isFirstValidNeighbour = false;

          } else {

            mPrimaryVertexContext->getRoads().emplace_back(iLayer, iCell);
          }

          traverseCellsTree(neighbourCellId, iLayer - 1);
        }

        // TODO: crosscheck for short track iterations
        // currentCell.setLevel(0);
      }
    }
#ifdef CA_DEBUG
    nRoads += mPrimaryVertexContext->getRoads().size();
    std::cout << "+++ Roads with " << iLevel + 2 << " clusters: " << nRoads << " / " << mPrimaryVertexContext->getRoads().size() << std::endl;
#endif
  }
}

void Tracker::findTracks(const ROframe& event)
{
  mTracks.reserve(mTracks.capacity() + mPrimaryVertexContext->getRoads().size());
  std::vector<TrackITSExt> tracks;
  std::vector<intermediateChi2> chi2s;
  tracks.reserve(mPrimaryVertexContext->getRoads().size());
  chi2s.reserve(mPrimaryVertexContext->getRoads().size());

#ifdef CA_DEBUG
  std::vector<int> roadCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> fitCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> backpropagatedCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> refitCounters(mTrkParams[0].NLayers - 3, 0);
  std::vector<int> nonsharingCounters(mTrkParams[0].NLayers - 3, 0);
#endif

  for (auto& road : mPrimaryVertexContext->getRoads()) {
    std::vector<int> clusters(mTrkParams[0].NLayers, constants::its::UnusedIndex);
    int lastCellLevel = constants::its::UnusedIndex;
    CA_DEBUGGER(int nClusters = 2);

    for (int iCell{0}; iCell < mTrkParams[0].CellsPerRoad(); ++iCell) {
      const int cellIndex = road[iCell];
      if (cellIndex == constants::its::UnusedIndex) {
        continue;
      } else {
        clusters[iCell] = mPrimaryVertexContext->getCells()[iCell][cellIndex].getFirstClusterIndex();
        clusters[iCell + 1] = mPrimaryVertexContext->getCells()[iCell][cellIndex].getSecondClusterIndex();
        clusters[iCell + 2] = mPrimaryVertexContext->getCells()[iCell][cellIndex].getThirdClusterIndex();
        assert(clusters[iCell] != constants::its::UnusedIndex &&
               clusters[iCell + 1] != constants::its::UnusedIndex &&
               clusters[iCell + 2] != constants::its::UnusedIndex);
        lastCellLevel = iCell;
        CA_DEBUGGER(nClusters++);
      }
    }

    CA_DEBUGGER(assert(nClusters >= mTrkParams[0].MinTrackLength));
    CA_DEBUGGER(roadCounters[nClusters - 4]++);

    if (lastCellLevel == constants::its::UnusedIndex) {
      continue;
    }

    /// From primary vertex context index to event index (== the one used as input of the tracking code)
    for (int iC{0}; iC < clusters.size(); iC++) {
      if (clusters[iC] != constants::its::UnusedIndex) {
        clusters[iC] = mPrimaryVertexContext->getClusters()[iC][clusters[iC]].clusterId;
      }
    }

    /// Track seed preparation. Clusters are numbered progressively from the outermost to the innermost.
    const auto& cluster1_glo = event.getClustersOnLayer(lastCellLevel + 2).at(clusters[lastCellLevel + 2]);
    const auto& cluster2_glo = event.getClustersOnLayer(lastCellLevel + 1).at(clusters[lastCellLevel + 1]);
    const auto& cluster3_glo = event.getClustersOnLayer(lastCellLevel).at(clusters[lastCellLevel]);

    const auto& cluster3_tf = event.getTrackingFrameInfoOnLayer(lastCellLevel).at(clusters[lastCellLevel]);

    /// FIXME!
    TrackITSExt temporaryTrack{buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_glo, cluster3_tf)};
    for (size_t iC = 0; iC < clusters.size(); ++iC) {
      temporaryTrack.setExternalClusterIndex(iC, clusters[iC], clusters[iC] != constants::its::UnusedIndex);
    }
    intermediateChi2 tmpChi2s;
    bool fitSuccess = fitTrack(event, temporaryTrack, mTrkParams[0].NLayers - 4, -1, -1, tmpChi2s.firstIteration);
    if (!fitSuccess) {
      continue;
    }
    CA_DEBUGGER(fitCounters[nClusters - 4]++);
    temporaryTrack.resetCovariance();
    fitSuccess = fitTrack(event, temporaryTrack, 0, mTrkParams[0].NLayers, 1, tmpChi2s.secondIteration, mTrkParams[0].FitIterationMaxChi2[0]);
    if (!fitSuccess) {
      continue;
    }
    CA_DEBUGGER(backpropagatedCounters[nClusters - 4]++);
    temporaryTrack.getParamOut() = temporaryTrack;
    temporaryTrack.setChi2Out(temporaryTrack.getChi2());
    temporaryTrack.resetCovariance();

    fitSuccess = fitTrack(event, temporaryTrack, mTrkParams[0].NLayers - 1, -1, -1, tmpChi2s.thirdIteration, mTrkParams[0].FitIterationMaxChi2[1]);
    if (!fitSuccess) {
      continue;
    }
    FakeTrackInfo<7> fir{mPrimaryVertexContext, event, road};
    if (fir.nFakeClusters == 0) {
      mDebugRoad.emplace_back(fir, tmpChi2s);
      // LOG(FATAL) << " >>> label: " << fir.mainLabel;
    }
    CA_DEBUGGER(refitCounters[nClusters - 4]++);
    tracks.emplace_back(temporaryTrack);
    chi2s.emplace_back(tmpChi2s);
    CA_DEBUGGER(assert(nClusters == temporaryTrack.getNumberOfClusters()));
  }
  //mTraits->refitTracks(event.getTrackingFrameInfo(), tracks);

  std::sort(tracks.begin(), tracks.end(),
            [](TrackITSExt& track1, TrackITSExt& track2) { return track1.isBetterOutChi2(track2, 1.e6f); });

  for (int iTrack{0}; iTrack < tracks.size(); ++iTrack) {
    auto& track = tracks[iTrack];
    auto& chi2 = chi2s[iTrack];
    CA_DEBUGGER(int nClusters = 0);
    int nShared = 0;
    for (int iLayer{0}; iLayer < mTrkParams[0].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
        continue;
      }
      nShared += int(mPrimaryVertexContext->isClusterUsed(iLayer, track.getClusterIndex(iLayer)));
      CA_DEBUGGER(nClusters++);
    }
    if (nShared > mTrkParams[0].ClusterSharing) {
      continue;
    }
    for (int iLayer{0}; iLayer < mTrkParams[0].NLayers; ++iLayer) {
      if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
        continue;
      }
      mPrimaryVertexContext->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
    }
    mTracks.emplace_back(track);
    mIntermediateTrackChi2.emplace_back(chi2);
  }
}

bool Tracker::fitTrack(const ROframe& event, TrackITSExt& track, int start, int end, int step, std::array<float, 7>& tmpChi2, const float chi2cut)
{
  auto propInstance = o2::base::Propagator::Instance();
  track.setChi2(0);
  for (int iLayer{start}; iLayer != end; iLayer += step) {
    if (track.getClusterIndex(iLayer) == constants::its::UnusedIndex) {
      continue;
    }
    const TrackingFrameInfo& trackingHit = event.getTrackingFrameInfoOnLayer(iLayer).at(track.getClusterIndex(iLayer));

    if (!track.rotate(trackingHit.alphaTrackingFrame)) {
      return false;
    }

    if (!propInstance->propagateToX(track, trackingHit.xTrackingFrame, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mCorrType)) {
      return false;
    }

    auto predChi2{track.getPredictedChi2(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
    if (predChi2 > chi2cut) {
      return false;
    }
    track.setChi2(track.getChi2() + predChi2);
    if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
      return false;
    }
    tmpChi2[iLayer] = predChi2;
  }
  return true;
}

void Tracker::traverseCellsTree(const int currentCellId, const int currentLayerId)
{
  Cell& currentCell{mPrimaryVertexContext->getCells()[currentLayerId][currentCellId]};
  const int currentCellLevel = currentCell.getLevel();

  mPrimaryVertexContext->getRoads().back().addCell(currentLayerId, currentCellId);

  if (currentLayerId > 0 && currentCellLevel > 1) {
    const int cellNeighboursNum{static_cast<int>(
      mPrimaryVertexContext->getCellsNeighbours()[currentLayerId - 1][currentCellId].size())};
    bool isFirstValidNeighbour = true;

    for (int iNeighbourCell{0}; iNeighbourCell < cellNeighboursNum; ++iNeighbourCell) {

      const int neighbourCellId =
        mPrimaryVertexContext->getCellsNeighbours()[currentLayerId - 1][currentCellId][iNeighbourCell];
      const Cell& neighbourCell = mPrimaryVertexContext->getCells()[currentLayerId - 1][neighbourCellId];

      if (currentCellLevel - 1 != neighbourCell.getLevel()) {
        continue;
      }

      if (isFirstValidNeighbour) {
        isFirstValidNeighbour = false;
      } else {
        mPrimaryVertexContext->getRoads().push_back(mPrimaryVertexContext->getRoads().back());
      }

      traverseCellsTree(neighbourCellId, currentLayerId - 1);
    }
  }

  // TODO: crosscheck for short track iterations
  // currentCell.setLevel(0);
}

void Tracker::computeRoadsMClabels(const ROframe& event)
{
  /// Moore's Voting Algorithm
  if (!event.hasMCinformation()) {
    return;
  }

  mPrimaryVertexContext->initialiseRoadLabels();

  int roadsNum{static_cast<int>(mPrimaryVertexContext->getRoads().size())};

  for (int iRoad{0}; iRoad < roadsNum; ++iRoad) {

    Road& currentRoad{mPrimaryVertexContext->getRoads()[iRoad]};
    MCCompLabel maxOccurrencesValue{constants::its::UnusedIndex, constants::its::UnusedIndex,
                                    constants::its::UnusedIndex, false};
    int count{0};
    bool isFakeRoad{false};
    bool isFirstRoadCell{true};

    for (int iCell{0}; iCell < mTrkParams[0].CellsPerRoad(); ++iCell) {
      const int currentCellIndex{currentRoad[iCell]};

      if (currentCellIndex == constants::its::UnusedIndex) {
        if (isFirstRoadCell) {
          continue;
        } else {
          break;
        }
      }

      const Cell& currentCell{mPrimaryVertexContext->getCells()[iCell][currentCellIndex]};

      if (isFirstRoadCell) {

        const int cl0index{mPrimaryVertexContext->getClusters()[iCell][currentCell.getFirstClusterIndex()].clusterId};
        auto& cl0labs{event.getClusterLabels(iCell, cl0index)};
        maxOccurrencesValue = cl0labs;
        count = 1;

        const int cl1index{mPrimaryVertexContext->getClusters()[iCell + 1][currentCell.getSecondClusterIndex()].clusterId};
        const auto& cl1labs{event.getClusterLabels(iCell + 1, cl1index)};

        if (cl1labs == maxOccurrencesValue) {
          ++count;
        } else {
          maxOccurrencesValue = cl1labs;
          count = 1;
          isFakeRoad = true;
        }

        isFirstRoadCell = false;
      }

      const int cl2index{mPrimaryVertexContext->getClusters()[iCell + 2][currentCell.getThirdClusterIndex()].clusterId};
      const auto& cl2labs{event.getClusterLabels(iCell + 2, cl2index)};

      if (cl2labs == maxOccurrencesValue) {
        ++count;
      } else {
        --count;
        isFakeRoad = true;
      }

      if (count == 0) {
        maxOccurrencesValue = cl2labs;
        count = 1;
      }
    }

    mPrimaryVertexContext->setRoadLabel(iRoad, maxOccurrencesValue.getRawValue(), isFakeRoad);
  }
}

void Tracker::computeTracksMClabels(const ROframe& event)
{
  /// Moore's Voting Algorithm
  if (!event.hasMCinformation()) {
    return;
  }

  int tracksNum{static_cast<int>(mTracks.size())};

  for (auto& track : mTracks) {

    MCCompLabel maxOccurrencesValue{constants::its::UnusedIndex, constants::its::UnusedIndex,
                                    constants::its::UnusedIndex, false};
    int count{0};
    bool isFakeTrack{false};

    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index == constants::its::UnusedIndex) {
        continue;
      }
      const MCCompLabel& currentLabel = event.getClusterLabels(iCluster, index);
      if (currentLabel == maxOccurrencesValue) {
        ++count;
      } else {
        if (count != 0) { // only in the first iteration count can be 0 at this point
          isFakeTrack = true;
          --count;
        }
        if (count == 0) {
          maxOccurrencesValue = currentLabel;
          count = 1;
        }
      }
      track.setExternalClusterIndex(iCluster, event.getClusterExternalIndex(iCluster, index));
    }

    if (isFakeTrack) {
      maxOccurrencesValue.setFakeFlag();
    }
    mTrackLabels.emplace_back(maxOccurrencesValue);
  }
}

void Tracker::rectifyClusterIndices(const ROframe& event)
{
  int tracksNum{static_cast<int>(mTracks.size())};
  for (auto& track : mTracks) {
    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index != constants::its::UnusedIndex) {
        track.setExternalClusterIndex(iCluster, event.getClusterExternalIndex(iCluster, index));
      }
    }
  }
}

/// Clusters are given from outside inward (cluster1 is the outermost). The innermost cluster is given in the tracking
/// frame coordinates
/// whereas the others are referred to the global frame. This function is almost a clone of CookSeed, adapted to return
/// a TrackParCov
track::TrackParCov Tracker::buildTrackSeed(const Cluster& cluster1, const Cluster& cluster2,
                                           const Cluster& cluster3, const TrackingFrameInfo& tf3)
{
  const float ca = std::cos(tf3.alphaTrackingFrame), sa = std::sin(tf3.alphaTrackingFrame);
  const float x1 = cluster1.xCoordinate * ca + cluster1.yCoordinate * sa;
  const float y1 = -cluster1.xCoordinate * sa + cluster1.yCoordinate * ca;
  const float z1 = cluster1.zCoordinate;
  const float x2 = cluster2.xCoordinate * ca + cluster2.yCoordinate * sa;
  const float y2 = -cluster2.xCoordinate * sa + cluster2.yCoordinate * ca;
  const float z2 = cluster2.zCoordinate;
  const float x3 = tf3.xTrackingFrame;
  const float y3 = tf3.positionTrackingFrame[0];
  const float z3 = tf3.positionTrackingFrame[1];

  const float crv = math_utils::computeCurvature(x1, y1, x2, y2, x3, y3);
  const float x0 = math_utils::computeCurvatureCentreX(x1, y1, x2, y2, x3, y3);
  const float tgl12 = math_utils::computeTanDipAngle(x1, y1, x2, y2, z1, z2);
  const float tgl23 = math_utils::computeTanDipAngle(x2, y2, x3, y3, z2, z3);

  const float fy = 1. / (cluster2.rCoordinate - cluster3.rCoordinate);
  const float& tz = fy;
  const float cy = (math_utils::computeCurvature(x1, y1, x2, y2 + constants::its::Resolution, x3, y3) - crv) /
                   (constants::its::Resolution * getBz() * o2::constants::math::B2C) *
                   20.f; // FIXME: MS contribution to the cov[14] (*20 added)
  constexpr float s2 = constants::its::Resolution * constants::its::Resolution;

  return track::TrackParCov(tf3.xTrackingFrame, tf3.alphaTrackingFrame,
                            {y3, z3, crv * (x3 - x0), 0.5f * (tgl12 + tgl23),
                             std::abs(getBz()) < o2::constants::math::Almost0 ? o2::constants::math::Almost0
                                                                              : crv / (getBz() * o2::constants::math::B2C)},
                            {s2, 0.f, s2, s2 * fy, 0.f, s2 * fy * fy, 0.f, s2 * tz, 0.f, s2 * tz * tz, s2 * cy, 0.f,
                             s2 * fy * cy, 0.f, s2 * cy * cy});
}

void Tracker::getGlobalConfiguration()
{
  auto& tc = o2::its::TrackerParamConfig::Instance();
  if (tc.useMatCorrTGeo) {
    setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo);
  }
  setUseSmoother(tc.useKalmanSmoother);
}

// void Tracker::smoothTracks(const ROframe& event)
// {
//   for (int iTrack{0}; iTrack < mTracks.size(); ++iTrack) {
//     auto& track = mTracks[iTrack];
//     if (track.getNumberOfClusters() > 4) {
//       FakeTrackInfo<7> fi{mPrimaryVertexContext, event, track, false};
//       if (!fi.isFake) { // correct tracks 5c+
//         Smoother<7> sm{track, 3, event, getBz(), mCorrType};
//         auto& cluster = event.getClusters()[3][track.getClusterIndex(3)];
//         auto phiBin = mPrimaryVertexContext->mIndexTableUtils.getPhiBinIndex(math_utils::getNormalizedPhiCoordinate(cluster.phiCoordinate));
//         auto zBin = mPrimaryVertexContext->mIndexTableUtils.getZBinIndex(3, cluster.zCoordinate);
//         auto tableBin = mPrimaryVertexContext->mIndexTableUtils.getBinIndex(zBin, phiBin);
//         LOG(WARN) << "Cluster z: " << cluster.zCoordinate << " phi: " << cluster.phiCoordinate << " phi bin: "
//                   << phiBin << " z bin: " << zBin << " bin on idxTable: " << tableBin
//                   << " index pointed: " << mPrimaryVertexContext->getIndexTables()[2][tableBin]
//                   << " index pointed next: " << mPrimaryVertexContext->getIndexTables()[2][tableBin+1];
//         int4 rect = mTraits->getBinsRect(cluster, 2, cluster.zCoordinate, cluster.zCoordinate, 0.5f, 0.5f / constants::its2::LayersRCoordinate()[3]);
//         LOG(WARN) << "rect x=" << rect.x << " rect y=" << rect.y;
//         LOG(WARN) << "rect z=" << rect.z << " rect w=" << rect.w;
//         bool here = false;
//         if (rect.x == 0 && rect.y == 0 && rect.z == 0 && rect.w == 0) {
//           continue;
//         }
//         int phiBinsNum{rect.w - rect.y + 1};
//         if (phiBinsNum < 0) {
//           phiBinsNum += 128;
//         }
//         for (int iPhiBin{rect.y}, iPhiCount{0}; iPhiCount < phiBinsNum; iPhiBin = ++iPhiBin == 128 ? 0 : iPhiBin, iPhiCount++) {
//           const int firstBinIndex{mPrimaryVertexContext->mIndexTableUtils.getBinIndex(rect.x, iPhiBin)};
//           const int maxBinIndex{firstBinIndex + rect.z - rect.x + 1};
//           const int firstRowClusterIndex = mPrimaryVertexContext->getIndexTables()[2][firstBinIndex];
//           const int maxRowClusterIndex = mPrimaryVertexContext->getIndexTables()[2][maxBinIndex];
//           LOG(WARN) << "tClusterIndex: " << track.getClusterIndex(3) << " fClusterIndex: " << firstRowClusterIndex << " lClusterIndex: " << maxRowClusterIndex;

// for (int iCluster{firstBinIndex}; iCluster < maxBinIndex; ++iCluster) {
//   if (iCluster == track.getClusterIndex(3)) {
//     continue;
//   }
//   here = true;
//   LOG(WARN) << " iCluster: " << iCluster;
//   sm.testCluster(iCluster, event);
// }
// }
// }

// if (here) {
// LOG(FATAL) << "fatalizing.";
// }
// }

// if (fi.nFakeClusters < 2) { // Correct track
//   for (int iLayerToSmooth{3}; iLayerToSmooth < 4; ++iLayerToSmooth) {
//     if (fi.firstFakeIndex != -1 && fi.firstFakeIndex == iLayerToSmooth) { //
//       auto label = fi.mcLabels[iLayerToSmooth];
//       if (label.isSet() && label != fi.mainLabel) {
//       }
// int layerF = -1;

// for (int i{0}; i < 7; ++i) {
// if (fi.mcLabels[i].isSet() && fi.mcLabels[i] != fi.mainLabel) {
// layerF = i;
// break;
// }
// }
// LOG(WARN) << "Found track with one fake cluster!\nlabels: ";
//for (int il{0}; il < 7; ++il) {
//  std::cout << "(";
//
// for (auto& l : fi.mcLabels /*[il]*/) {
// std::cout << " " << l;
// }
// std::cout << ")\t";
// int correct = event.getFirstClusterIDFromLabel(layerF, fi.mainLabel);
// bool hasCorrect = correct != -1 ? true : false;
// std::cout << "\nFake cluster is on layer: " << layerF << ", id of correct cluster: " << correct << ", id of fake cluster: " << track.getClusterIndex(layerF);
// if (correct != -1) {
// std::cout << ", counter_proof: " << event.getClusterLabels(layerF, correct) << std::endl;
// Smoother<7> sm{track, layerF, event, getBz(), mCorrType};
// if (sm.isValidInit()) {
// float smChi2Initial{sm.getChi2()};
// bool better = sm.testCluster(correct, event);
// float smChi2Tested{-1.f};
// if (better) {
// smChi2Tested = sm.getLastChi2();
// }
// auto trackLength = track.getNumberOfClusters();
// auto startLayer = track.getFirstClusterLayer();
// mDebugger->dumpSmootherChi2(layerF, smChi2Initial, smChi2Tested, startLayer, trackLength);
// } else {
// std::cout << std::endl;
// }
// } else {
// std::cout << std::endl;
// }
// mDebugger->dumpLayerFake(layerF, hasCorrect);
//     }
//   }
// }
//   }
// }

void Tracker::smoothTracks(const ROframe& event)
{
  ////////////////////////////////////////////
  //      Debug finally promoted tracks     //
  ////////////////////////////////////////////
  for (int iTrack{0}; iTrack < mTracks.size(); ++iTrack) {
    auto& track = mTracks[iTrack];
    FakeTrackInfo<7> fi{mPrimaryVertexContext, event, track, false};
    int layerF = -1;
    if (fi.nFakeClusters == 1 && fi.firstFakeIndex == 6 && track.getNClusters() == 7) {
      for (int i{0}; i < 7; ++i) {
        if (fi.mcLabels[i].isSet() && fi.mcLabels[i] != fi.mainLabel) {
          layerF = i;
          break;
        }
      }
      LOG(WARN) << "Found track with one fake cluster!\nlabels: ";
      // for (int il{0}; il < 7; ++il) {
      //   std::cout << "(";

      for (auto& l : fi.mcLabels /*[il]*/) {
        std::cout << " " << l;
      }
      // std::cout << ")\t";
      int correct = event.getFirstClusterIDFromLabel(layerF, fi.mainLabel);
      bool hasCorrect = correct != -1 ? true : false;
      std::cout << "\nFake cluster is on layer: " << layerF << ", id of correct cluster: " << correct << ", id of fake cluster: " << track.getClusterIndex(layerF);
      if (correct != -1) {
        std::cout << ", proof: " << event.getClusterLabels(layerF, correct) << std::endl;
        bool foundRoad{false};
        std::pair<o2::its::FakeTrackInfo<7>, o2::its::intermediateChi2> road;
        for (auto iRoad{0}; iRoad < mDebugRoad.size(); ++iRoad) {
          road = mDebugRoad[iRoad];
          if (road.first.mainLabel == fi.mainLabel) {
            LOG(WARN) << "Found a match, dumping road labels...";
            for (auto& l : road.first.mcLabels /*[il]*/) {
              std::cout << " " << l;
            }
            std::cout << std::endl;
            if (road.first.occurrences.size() == 1 && road.first.occurrences[0].second == 7) {
              foundRoad = true;
              break;
            }
          }
        }
        if (!foundRoad) {
          LOG(WARN) << "Could not find road with 7 correct clusters, browsing cells";
          for (auto iCell{0}; iCell < mPrimaryVertexContext->getCells()[4].size(); ++iCell) {
            auto index_0 = mPrimaryVertexContext->getCells()[4][iCell].getFirstClusterIndex();
            auto index_1 = mPrimaryVertexContext->getCells()[4][iCell].getSecondClusterIndex();
            auto index_2 = mPrimaryVertexContext->getCells()[4][iCell].getThirdClusterIndex();
            // if (index_0 != constants::its::UnusedIndex && index_1 != constants::its::UnusedIndex && index_2 != constants::its::UnusedIndex) {
            index_0 = mPrimaryVertexContext->getClusters()[4][index_0].clusterId;
            index_1 = mPrimaryVertexContext->getClusters()[5][index_1].clusterId;
            index_2 = mPrimaryVertexContext->getClusters()[6][index_2].clusterId;

            o2::MCCompLabel mcLabel_0 = event.getClusterLabels(4, index_0);
            o2::MCCompLabel mcLabel_1 = event.getClusterLabels(5, index_1);
            o2::MCCompLabel mcLabel_2 = event.getClusterLabels(6, index_2);

            LOG(WARN) << "Processing: " << mcLabel_0 << " " << mcLabel_1 << " " << mcLabel_2;
            if (mcLabel_0 == mcLabel_1 && mcLabel_0 == mcLabel_2 && mcLabel_0 == fi.mainLabel) {
              LOG(WARN) << "Cell exists, " << mcLabel_0 << " id: " << iCell << " cluster ids: " << index_0 << " " << index_1 << " " << index_2;
            } else {
              LOG(WARN) << "Cell does not exist!";
            }
            // }
          }
          LOG(FATAL) << "fatalizing...";
        } else {
          if (road.second.secondIteration[6] < mIntermediateTrackChi2[iTrack].secondIteration[6]) {
            LOG(WARN) << "mismatch: " << road.second.secondIteration[6] << " vs: " << mIntermediateTrackChi2[iTrack].secondIteration[6];

            float sumRoadFirst{0.f};
            float sumRoadSecond{0.f};
            float sumRoadThird{0.f};
            float sumWinnerFirst{0.f};
            float sumWinnerSecond{0.f};
            float sumWinnerThird{0.f};
            // 1
            std::cout << "1 iteration, MC road: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << road.second.firstIteration[i] << " | ";
              sumRoadFirst += road.second.firstIteration[i];
            }
            std::cout << " sum: " << sumRoadFirst << std::endl;
            std::cout << "1 iteration, W track: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << mIntermediateTrackChi2[iTrack].firstIteration[i] << " | ";
              sumWinnerFirst += mIntermediateTrackChi2[iTrack].firstIteration[i];
            }
            std::cout << " sum: " << sumWinnerFirst << std::endl;
            // 2
            std::cout << "2 iteration, MC road: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << road.second.secondIteration[i] << " | ";
              sumRoadSecond += road.second.secondIteration[i];
            }
            std::cout << " sum: " << sumRoadSecond << std::endl;
            std::cout << "2 iteration, W track: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << mIntermediateTrackChi2[iTrack].secondIteration[i] << " ";
              sumWinnerSecond += mIntermediateTrackChi2[iTrack].secondIteration[i];
            }
            std::cout << " sum: " << sumWinnerSecond << std::endl;
            // 3
            std::cout << "3 iteration, MC road: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << road.second.thirdIteration[i] << " | ";
              sumRoadThird += road.second.thirdIteration[i];
            }
            std::cout << " sum: " << sumRoadThird << std::endl;
            std::cout << "3 iteration, W track: ";
            for (size_t i{0}; i < 7; ++i) {
              std::cout << mIntermediateTrackChi2[iTrack].thirdIteration[i] << " | ";
              sumWinnerThird += mIntermediateTrackChi2[iTrack].thirdIteration[i];
            }
            std::cout << " sum: " << sumWinnerThird << std::endl;

            LOG(FATAL) << "fatalizing ...";
          }
        }

        Smoother<7> sm{track, layerF, event, getBz(), mCorrType};
        if (sm.isValidInit()) {
          float smChi2Initial{sm.getChi2()};
          bool better = sm.testCluster(correct, event);
          float smChi2Tested{-1.f};
          if (better) {
            smChi2Tested = sm.getLastChi2();
          }
          auto trackLength = track.getNumberOfClusters();
          auto startLayer = track.getFirstClusterLayer();
          mDebugger->dumpSmootherChi2(layerF, smChi2Initial, smChi2Tested, startLayer, trackLength);
        } else {
          std::cout << std::endl;
        }
      } else {
        std::cout << std::endl;
      }
      mDebugger->dumpLayerFake(layerF, hasCorrect);
    }
  }
}
} // namespace its
} // namespace o2
