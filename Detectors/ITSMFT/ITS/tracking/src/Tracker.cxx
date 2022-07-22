// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
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
#include "ITStracking/TrackingConfigParam.h"

#include "ReconstructionDataFormats/Track.h"
#include <cassert>
#include <iostream>
#include <dlfcn.h>
#include <cstdlib>
#include <string>
#include <climits>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace o2
{
namespace its
{

Tracker::Tracker(o2::its::TrackerTraits* traits)
{
  /// Initialise standard configuration with 1 iteration
  mTrkParams.resize(1);
  mTraits = traits;
}

Tracker::~Tracker() = default;

void Tracker::clustersToTracks(std::function<void(std::string s)> logger, std::function<void(std::string s)> error)
{
  double total{0};
  mTraits->UpdateTrackingParameters(mTrkParams);
  for (int iteration = 0; iteration < (int)mTrkParams.size(); ++iteration) {
    total += evaluateTask(&Tracker::initialiseTimeFrame, "Timeframe initialisation", logger, iteration);
    total += evaluateTask(&Tracker::computeTracklets, "Tracklet finding", logger, iteration);
    logger(fmt::format("\t- Number of tracklets: {}", mTimeFrame->getNumberOfTracklets()));
    if (!mTimeFrame->checkMemory(mTrkParams[iteration].MaxMemory)) {
      error("Too much memory used during trackleting, check the detector status and/or the selections.");
      break;
    }
    float trackletsPerCluster = mTimeFrame->getNumberOfClusters() > 0 ? float(mTimeFrame->getNumberOfTracklets()) / mTimeFrame->getNumberOfClusters() : 0.f;
    if (trackletsPerCluster > mTrkParams[iteration].TrackletsPerClusterLimit) {
      error(fmt::format("Too many tracklets per cluster ({}), check the detector status and/or the selections.", trackletsPerCluster));
      break;
    }

    total += evaluateTask(&Tracker::computeCells, "Cell finding", logger, iteration);
    logger(fmt::format("\t- Number of Cells: {}", mTimeFrame->getNumberOfCells()));
    if (!mTimeFrame->checkMemory(mTrkParams[iteration].MaxMemory)) {
      error("Too much memory used during cell finding, check the detector status and/or the selections.");
      break;
    }
    float cellsPerCluster = mTimeFrame->getNumberOfClusters() > 0 ? float(mTimeFrame->getNumberOfCells()) / mTimeFrame->getNumberOfClusters() : 0.f;
    if (cellsPerCluster > mTrkParams[iteration].CellsPerClusterLimit) {
      error(fmt::format("Too many cells per cluster ({}), check the detector status and/or the selections.", cellsPerCluster));
      break;
    }

    total += evaluateTask(&Tracker::findCellsNeighbours, "Neighbour finding", logger, iteration);
    total += evaluateTask(&Tracker::findRoads, "Road finding", logger, iteration);
    logger(fmt::format("\t- Number of Roads: {}", mTimeFrame->getRoads().size()));
    total += evaluateTask(&Tracker::findTracks, "Track finding", logger, iteration);
    total += evaluateTask(&Tracker::extendTracks, "Extending tracks", logger, iteration);
  }

  total += evaluateTask(&Tracker::findShortPrimaries, "Short primaries finding", logger);
  /// TODO: Add desperate tracking, aka the extension of short primaries to recover holes in layer 3

  std::stringstream sstream;
  if (constants::DoTimeBenchmarks) {
    sstream << std::setw(2) << " - "
            << "Timeframe " << mTimeFrameCounter++ << " processing completed in: " << total << "ms using " << mNThreads << " threads.";
  }
  logger(sstream.str());

  if (mTimeFrame->hasMCinformation()) {
    computeTracksMClabels();
  }
  rectifyClusterIndices();
  mNumberOfRuns++;
}

void Tracker::computeTracklets(int& iteration)
{
  mTraits->computeLayerTracklets(iteration);
}

void Tracker::computeCells(int& iteration)
{
  mTraits->computeLayerCells(iteration);
}

void Tracker::findCellsNeighbours(int& iteration)
{
  mTraits->findCellsNeighbours(iteration);
}

void Tracker::findRoads(int& iteration)
{
  mTraits->findRoads(iteration);
}

void Tracker::findTracks(int& iteration)
{
  mTraits->findTracks(iteration);
}

void Tracker::extendTracks(int& iteration)
{
  mTraits->extendTracks(iteration);
}

void Tracker::findShortPrimaries()
{
  mTraits->findShortPrimaries();
}

void Tracker::computeRoadsMClabels()
{
  /// Moore's Voting Algorithm
  if (!mTimeFrame->hasMCinformation()) {
    return;
  }

  mTimeFrame->initialiseRoadLabels();

  int roadsNum{static_cast<int>(mTimeFrame->getRoads().size())};

  for (int iRoad{0}; iRoad < roadsNum; ++iRoad) {

    Road& currentRoad{mTimeFrame->getRoads()[iRoad]};
    std::vector<std::pair<MCCompLabel, size_t>> occurrences;
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

      const Cell& currentCell{mTimeFrame->getCells()[iCell][currentCellIndex]};

      if (isFirstRoadCell) {

        const int cl0index{mTimeFrame->getClusters()[iCell][currentCell.getFirstClusterIndex()].clusterId};
        auto cl0labs{mTimeFrame->getClusterLabels(iCell, cl0index)};
        bool found{false};
        for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
          std::pair<o2::MCCompLabel, size_t>& occurrence = occurrences[iOcc];
          for (auto& label : cl0labs) {
            if (label == occurrence.first) {
              ++occurrence.second;
              found = true;
              // break; // uncomment to stop to the first hit
            }
          }
        }
        if (!found) {
          for (auto& label : cl0labs) {
            occurrences.emplace_back(label, 1);
          }
        }

        const int cl1index{mTimeFrame->getClusters()[iCell + 1][currentCell.getSecondClusterIndex()].clusterId};

        const auto& cl1labs{mTimeFrame->getClusterLabels(iCell + 1, cl1index)};
        found = false;
        for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
          std::pair<o2::MCCompLabel, size_t>& occurrence = occurrences[iOcc];
          for (auto& label : cl1labs) {
            if (label == occurrence.first) {
              ++occurrence.second;
              found = true;
              // break; // uncomment to stop to the first hit
            }
          }
        }
        if (!found) {
          for (auto& label : cl1labs) {
            occurrences.emplace_back(label, 1);
          }
        }

        isFirstRoadCell = false;
      }

      const int cl2index{mTimeFrame->getClusters()[iCell + 2][currentCell.getThirdClusterIndex()].clusterId};
      const auto& cl2labs{mTimeFrame->getClusterLabels(iCell + 2, cl2index)};
      bool found{false};
      for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
        std::pair<o2::MCCompLabel, size_t>& occurrence = occurrences[iOcc];
        for (auto& label : cl2labs) {
          if (label == occurrence.first) {
            ++occurrence.second;
            found = true;
            // break; // uncomment to stop to the first hit
          }
        }
      }
      if (!found) {
        for (auto& label : cl2labs) {
          occurrences.emplace_back(label, 1);
        }
      }
    }

    std::sort(occurrences.begin(), occurrences.end(), [](auto e1, auto e2) {
      return e1.second > e2.second;
    });

    auto maxOccurrencesValue = occurrences[0].first;
    mTimeFrame->setRoadLabel(iRoad, maxOccurrencesValue.getRawValue(), isFakeRoad);
  }
}

void Tracker::computeTracksMClabels()
{
  for (int iROF{0}; iROF < mTimeFrame->getNrof(); ++iROF) {
    for (auto& track : mTimeFrame->getTracks(iROF)) {
      std::vector<std::pair<MCCompLabel, size_t>> occurrences;
      occurrences.clear();

      for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
        const int index = track.getClusterIndex(iCluster);
        if (index == constants::its::UnusedIndex) {
          continue;
        }
        auto labels = mTimeFrame->getClusterLabels(iCluster, index);
        bool found{false};
        for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
          std::pair<o2::MCCompLabel, size_t>& occurrence = occurrences[iOcc];
          for (auto& label : labels) {
            if (label == occurrence.first) {
              ++occurrence.second;
              found = true;
              // break; // uncomment to stop to the first hit
            }
          }
        }
        if (!found) {
          for (auto& label : labels) {
            occurrences.emplace_back(label, 1);
          }
        }
      }
      std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) {
        return e1.second > e2.second;
      });

      auto maxOccurrencesValue = occurrences[0].first;
      uint32_t pattern = track.getPattern();
      // set fake clusters pattern
      for (int ic{TrackITSExt::MaxClusters}; ic--;) {
        auto clid = track.getClusterIndex(ic);
        if (clid != constants::its::UnusedIndex) {
          auto labelsSpan = mTimeFrame->getClusterLabels(ic, clid);
          for (auto& currentLabel : labelsSpan) {
            if (currentLabel == maxOccurrencesValue) {
              pattern |= 0x1 << (16 + ic); // set bit if correct
              break;
            }
          }
        }
      }
      track.setPattern(pattern);
      if (occurrences[0].second < track.getNumberOfClusters()) {
        maxOccurrencesValue.setFakeFlag();
      }
      mTimeFrame->getTracksLabel(iROF).emplace_back(maxOccurrencesValue);
    }
  }
}

void Tracker::rectifyClusterIndices()
{
  for (int iROF{0}; iROF < mTimeFrame->getNrof(); ++iROF) {
    for (auto& track : mTimeFrame->getTracks(iROF)) {
      for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
        const int index = track.getClusterIndex(iCluster);
        if (index != constants::its::UnusedIndex) {
          track.setExternalClusterIndex(iCluster, mTimeFrame->getClusterExternalIndex(iCluster, index));
        }
      }
    }
  }
}

void Tracker::getGlobalConfiguration()
{
  auto& tc = o2::its::TrackerParamConfig::Instance();
  if (tc.useMatCorrTGeo) {
    mTraits->setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo);
  } else if (tc.useFastMaterial) {
    mTraits->setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE);
  } else {
    mTraits->setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT);
  }
  setNThreads(tc.nThreads);
  for (auto& params : mTrkParams) {
    if (params.NLayers == 7) {
      for (int i{0}; i < 7; ++i) {
        params.LayerMisalignment[i] = tc.sysErrZ2[i] > 0 ? std::sqrt(tc.sysErrZ2[i]) : params.LayerMisalignment[i];
      }
    }
    params.PhiBins = tc.LUTbinsPhi > 0 ? tc.LUTbinsPhi : params.PhiBins;
    params.ZBins = tc.LUTbinsZ > 0 ? tc.LUTbinsZ : params.ZBins;
    params.PVres = tc.pvRes > 0 ? tc.pvRes : params.PVres;
    params.NSigmaCut *= tc.nSigmaCut > 0 ? tc.nSigmaCut : 1.f;
    params.CellDeltaTanLambdaSigma *= tc.deltaTanLres > 0 ? tc.deltaTanLres : 1.f;
    params.TrackletMinPt *= tc.minPt > 0 ? tc.minPt : 1.f;
    for (int iD{0}; iD < 3; ++iD) {
      params.Diamond[iD] = tc.diamondPos[iD];
    }
    params.UseDiamond = tc.useDiamond;
    if (tc.maxMemory) {
      params.MaxMemory = tc.maxMemory;
    }
    if (tc.useTrackFollower >= 0) {
      params.UseTrackFollower = tc.useTrackFollower;
    }
    if (tc.cellsPerClusterLimit >= 0) {
      params.CellsPerClusterLimit = tc.cellsPerClusterLimit;
    }
    if (tc.trackletsPerClusterLimit >= 0) {
      params.TrackletsPerClusterLimit = tc.trackletsPerClusterLimit;
    }
    if (tc.findShortTracks >= 0) {
      params.FindShortTracks = tc.findShortTracks;
    }
  }
}

void Tracker::adoptTimeFrame(TimeFrame& tf)
{
  mTimeFrame = &tf;
  mTraits->adoptTimeFrame(&tf);
}

void Tracker::setBz(float bz)
{
  mTraits->setBz(bz);
}

void Tracker::setCorrType(const o2::base::PropagatorImpl<float>::MatCorrType& type)
{
  mTraits->setCorrType(type);
}

bool Tracker::isMatLUT() const
{
  return mTraits->isMatLUT();
}

void Tracker::setNThreads(int n)
{
#ifdef WITH_OPENMP
  mNThreads = n > 0 ? n : 1;
#else
  mNThreads = 1;
#endif
}

} // namespace its
} // namespace o2
