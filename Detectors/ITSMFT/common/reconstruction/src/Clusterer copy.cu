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

/// \file Clusterer.cxx
/// \brief Implementation of the ITS cluster finder
#include <algorithm>
#include <TTree.h>
// #include "Framework/Logger.h"
// #include "ITSMFTBase/GeometryTGeo.h"
#include "ITSMFTReconstruction/Clusterer.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "CommonDataFormat/InteractionRecord.h"
#include <cuda.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif
using namespace o2::itsmft;

namespace o2::itsmft::gpu
{
// Empty CUDA kernel
__global__ void empty_kernel()
{
  printf("Hello world from gpu\n");
}
} // namespace o2::itsmft::gpu

//__________________________________________________
void Clusterer::process(int nThreads, PixelReader& reader, CompClusCont* compClus,
                        PatternCont* patterns, ROFRecCont* vecROFRec, MCTruth* labelsCl)
{
#ifdef _PERFORM_TIMING_
  mTimer.Start(kFALSE);
#endif
  if (nThreads < 1) {
    nThreads = 1;
  }
  auto autoDecode = reader.getDecodeNextAuto();
  int rofcount{0};
  do {
    if (autoDecode) {
      reader.setDecodeNextAuto(false); // internally do not autodecode
      if (!reader.decodeNextTrigger()) {
        break;                         // on the fly decoding was requested, but there were no data left
      }
    }
    if (reader.getInteractionRecord().isDummy()) {
      continue; // No IR info was found
    }
    // pre-fetch all non-empty chips of current ROF
    ChipPixelData* curChipData = nullptr;
    mFiredChipsPtr.clear();

    // NOTE: should be total number of pixels fired in ROF
    size_t nPix = 0;

    // NOTE: cycle until no more chips fired in ROF 

    // NOTE: not really sure how getNextChipData works tbh. should not matter for now
    // NOTE: does it only return chips that fired?
    while ((curChipData = reader.getNextChipData(mChips))) {
      mFiredChipsPtr.push_back(curChipData);
      nPix += curChipData->getData().size();
    }

    auto& rof = vecROFRec->emplace_back(reader.getInteractionRecord(), vecROFRec->size(), compClus->size(), 0); // create new ROF

    uint16_t nFired = mFiredChipsPtr.size();
    if (!nFired) {
      if (autoDecode) {
        continue;
      }
      break; // just 1 ROF was asked to be processed
    }
    // NOTE: probably start (at least) one thread per chip. less threads needed if less chips fired
    if (nFired < nThreads) {
      nThreads = nFired;
    }
#ifndef WITH_OPENMP
    nThreads = 1;
#endif
    uint16_t chipStep = nThreads > 1 ? (nThreads == 2 ? 20 : 10) : nFired;
    int dynGrp = std::min(4, std::max(1, nThreads / 2));
    if (nThreads > mThreads.size()) {
      int oldSz = mThreads.size();
      mThreads.resize(nThreads);
      for (int i = oldSz; i < nThreads; i++) {
        mThreads[i] = std::make_unique<ClustererThread>(this, i);
      }
    }
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic, dynGrp) num_threads(nThreads)
    //>> start of MT region
    for (uint16_t ic = 0; ic < nFired; ic += chipStep) {
      auto ith = omp_get_thread_num();
      if (nThreads > 1) {
        mThreads[ith]->process(ic, std::min(chipStep, uint16_t(nFired - ic)),
                               &mThreads[ith]->compClusters,
                               patterns ? &mThreads[ith]->patterns : nullptr,
                               labelsCl ? reader.getDigitsMCTruth() : nullptr,
                               labelsCl ? &mThreads[ith]->labels : nullptr, rof);
      } else { // put directly to the destination
        mThreads[0]->process(0, nFired, compClus, patterns, labelsCl ? reader.getDigitsMCTruth() : nullptr, labelsCl, rof);
      }
    }
    //<< end of MT region
#else
    mThreads[0]->process(0, nFired, compClus, patterns, labelsCl ? reader.getDigitsMCTruth() : nullptr, labelsCl, rof);
#endif
    // copy data of all threads but the 1st one to final destination
    if (nThreads > 1) {
#ifdef _PERFORM_TIMING_
      mTimerMerge.Start(false);
#endif
      size_t nClTot = 0, nPattTot = 0;
      int chid = 0, thrStatIdx[nThreads];
      for (int ith = 0; ith < nThreads; ith++) {
        std::sort(mThreads[ith]->stats.begin(), mThreads[ith]->stats.end(), [](const ThreadStat& a, const ThreadStat& b) { return a.firstChip < b.firstChip; });
        thrStatIdx[ith] = 0;
        nClTot += mThreads[ith]->compClusters.size();
        nPattTot += mThreads[ith]->patterns.size();
      }
      compClus->reserve(nClTot);
      if (patterns) {
        patterns->reserve(nPattTot);
      }
      while (chid < nFired) {
        for (int ith = 0; ith < nThreads; ith++) {
          if (thrStatIdx[ith] >= mThreads[ith]->stats.size()) {
            continue;
          }
          const auto& stat = mThreads[ith]->stats[thrStatIdx[ith]];
          if (stat.firstChip == chid) {
            thrStatIdx[ith]++;
            chid += stat.nChips; // next chip to look
            const auto clbeg = mThreads[ith]->compClusters.begin() + stat.firstClus;
            auto szold = compClus->size();
            compClus->insert(compClus->end(), clbeg, clbeg + stat.nClus);
            if (patterns) {
              const auto ptbeg = mThreads[ith]->patterns.begin() + stat.firstPatt;
              patterns->insert(patterns->end(), ptbeg, ptbeg + stat.nPatt);
            }
            if (labelsCl) {
              labelsCl->mergeAtBack(mThreads[ith]->labels, stat.firstClus, stat.nClus);
            }
          }
        }
      }
      for (int ith = 0; ith < nThreads; ith++) {
        mThreads[ith]->patterns.clear();
        mThreads[ith]->compClusters.clear();
        mThreads[ith]->labels.clear();
        mThreads[ith]->stats.clear();
      }
#ifdef _PERFORM_TIMING_
      mTimerMerge.Stop();
#endif
    } else {
      mThreads[0]->stats.clear();
    }
    rof.setNEntries(compClus->size() - rof.getFirstEntry()); // update
  } while (autoDecode);
  reader.setDecodeNextAuto(autoDecode);                      // restore setting
#ifdef _PERFORM_TIMING_
  mTimer.Stop();
  LOGP(info, "Time to finish: {}s", mTimer.RealTime());
#endif
  std::vector<Clusterer::BBox> outClusterBBoxes;
  std::vector<std::vector<PixelData>> outClusterPixels;
  mToyProblem.executeToyProblem(outClusterBBoxes, outClusterPixels);
  std::cout << "TOY PROBLEM PRODUCED: " << outClusterBBoxes.size() << " BOUNDING BOXES." << std::endl;
  std::cout << "TOY PROBLEM PRODUCED: " << outClusterPixels.size() << " PIXEL LISTS." << std::endl;
}

//__________________________________________________
void Clusterer::ClustererThread::process(uint16_t chip, uint16_t nChips, CompClusCont* compClusPtr, PatternCont* patternsPtr,
                                         const ConstMCTruth* labelsDigPtr, MCTruth* labelsClPtr, const ROFRecord& rofPtr)
{
  if (stats.empty() || stats.back().firstChip + stats.back().nChips != chip) { // there is a jump, register new block
    stats.emplace_back(ThreadStat{chip, 0, uint32_t(compClusPtr->size()), patternsPtr ? uint32_t(patternsPtr->size()) : 0, 0, 0});
  }
  for (int ic = 0; ic < nChips; ic++) {
    auto* curChipData = parent->mFiredChipsPtr[chip + ic];
    #pragma omp critical
    {
      parent->mToyProblem.addChipAsync(curChipData);
    }
    auto chipID = curChipData->getChipID();
    if (parent->mMaxBCSeparationToMask > 0) { // mask pixels fired from the previous ROF
      const auto& chipInPrevROF = parent->mChipsOld[chipID];
      if (std::abs(rofPtr.getBCData().differenceInBC(chipInPrevROF.getInteractionRecord())) < parent->mMaxBCSeparationToMask) {
        parent->mMaxRowColDiffToMask ? curChipData->maskFiredInSample(parent->mChipsOld[chipID], parent->mMaxRowColDiffToMask) : curChipData->maskFiredInSample(parent->mChipsOld[chipID]);
      }
    }
    auto nclus0 = compClusPtr->size();
    auto validPixID = curChipData->getFirstUnmasked();
    auto npix = curChipData->getData().size();
    if (validPixID < npix) {    // chip data may have all of its pixels masked!
      auto valp = validPixID++;
      if (validPixID == npix) { // special case of a single pixel fired on the chip
        finishChipSingleHitFast(valp, curChipData, compClusPtr, patternsPtr, labelsDigPtr, labelsClPtr);
      } else {
        initChip(curChipData, valp);
        for (; validPixID < npix; validPixID++) {
          if (!curChipData->getData()[validPixID].isMasked()) {
            updateChip(curChipData, validPixID);
          }
        }
        finishChip(curChipData, compClusPtr, patternsPtr, labelsDigPtr, labelsClPtr);
      }
    }
    if (parent->mMaxBCSeparationToMask > 0) { // current chip data will be used in the next ROF to mask overflow pixels
      parent->mChipsOld[chipID].swap(*curChipData);
    }
  }
  auto& currStat = stats.back();
  currStat.nChips += nChips;
  currStat.nClus = compClusPtr->size() - currStat.firstClus;
  currStat.nPatt = patternsPtr ? (patternsPtr->size() - currStat.firstPatt) : 0;
}

//__________________________________________________
void Clusterer::ClustererThread::finishChip(ChipPixelData* curChipData, CompClusCont* compClusPtr,
                                            PatternCont* patternsPtr, const ConstMCTruth* labelsDigPtr, MCTruth* labelsClusPtr)
{
  const auto& pixData = curChipData->getData();
  for (int i1 = 0; i1 < preClusterHeads.size(); ++i1) {
    auto ci = preClusterIndices[i1];
    if (ci < 0) {
      continue;
    }
    BBox bbox(curChipData->getChipID());
    int nlab = 0;
    int next = preClusterHeads[i1];
    pixArrBuff.clear();
    while (next >= 0) {
      const auto& pixEntry = pixels[next];
      const auto pix = pixData[pixEntry.second];
      pixArrBuff.push_back(pix); // needed for cluster topology
      bbox.adjust(pix.getRowDirect(), pix.getCol());
      if (labelsClusPtr) {
        if (parent->mSquashingDepth) { // the MCtruth for this pixel is stored in chip data: due to squashing we lose contiguity
          fetchMCLabels(curChipData->getOrderedPixId(pixEntry.second), labelsDigPtr, nlab);
        } else {                       // the MCtruth for this pixel is at curChipData->startID+pixEntry.second
          fetchMCLabels(pixEntry.second + curChipData->getStartID(), labelsDigPtr, nlab);
        }
      }
      next = pixEntry.first;
    }
    preClusterIndices[i1] = -1;
    for (int i2 = i1 + 1; i2 < preClusterHeads.size(); ++i2) {
      if (preClusterIndices[i2] != ci) {
        continue;
      }
      next = preClusterHeads[i2];
      while (next >= 0) {
        const auto& pixEntry = pixels[next];
        const auto pix = pixData[pixEntry.second]; // PixelData
        pixArrBuff.push_back(pix);                 // needed for cluster topology
        bbox.adjust(pix.getRowDirect(), pix.getCol());
        if (labelsClusPtr) {
          if (parent->mSquashingDepth) { // the MCtruth for this pixel is stored in chip data: due to squashing we lose contiguity
            fetchMCLabels(curChipData->getOrderedPixId(pixEntry.second), labelsDigPtr, nlab);
          } else {                       // the MCtruth for this pixel is at curChipData->startID+pixEntry.second
            fetchMCLabels(pixEntry.second + curChipData->getStartID(), labelsDigPtr, nlab);
          }
        }
        next = pixEntry.first;
      }
      preClusterIndices[i2] = -1;
    }
    if (bbox.isAcceptableSize()) {
      parent->streamCluster(pixArrBuff, &labelsBuff, bbox, parent->mPattIdConverter, compClusPtr, patternsPtr, labelsClusPtr, nlab);
    } else {
      auto warnLeft = MaxHugeClusWarn - parent->mNHugeClus;
      if (warnLeft > 0) {
        LOGP(warn, "Splitting a huge cluster: chipID {}, rows {}:{} cols {}:{}{}", bbox.chipID, bbox.rowMin, bbox.rowMax, bbox.colMin, bbox.colMax,
             warnLeft == 1 ? " (Further warnings will be muted)" : "");
#ifdef WITH_OPENMP
#pragma omp critical
#endif
        {
          parent->mNHugeClus++;
        }
      }
      BBox bboxT(bbox); // truncated box
      std::vector<PixelData> pixbuf;
      do {
        bboxT.rowMin = bbox.rowMin;
        bboxT.colMax = std::min(bbox.colMax, uint16_t(bboxT.colMin + o2::itsmft::ClusterPattern::MaxColSpan - 1));
        do { // Select a subset of pixels fitting the reduced bounding box
          bboxT.rowMax = std::min(bbox.rowMax, uint16_t(bboxT.rowMin + o2::itsmft::ClusterPattern::MaxRowSpan - 1));
          for (const auto& pix : pixArrBuff) {
            if (bboxT.isInside(pix.getRowDirect(), pix.getCol())) {
              pixbuf.push_back(pix);
            }
          }
          if (!pixbuf.empty()) { // Stream a piece of cluster only if the reduced bounding box is not empty
            parent->streamCluster(pixbuf, &labelsBuff, bboxT, parent->mPattIdConverter, compClusPtr, patternsPtr, labelsClusPtr, nlab, true);
            pixbuf.clear();
          }
          bboxT.rowMin = bboxT.rowMax + 1;
        } while (bboxT.rowMin < bbox.rowMax);
        bboxT.colMin = bboxT.colMax + 1;
      } while (bboxT.colMin < bbox.colMax);
    }
  }
}

//__________________________________________________
void Clusterer::ClustererThread::finishChipSingleHitFast(uint32_t hit, ChipPixelData* curChipData, CompClusCont* compClusPtr,
                                                         PatternCont* patternsPtr, const ConstMCTruth* labelsDigPtr, MCTruth* labelsClusPtr)
{
  auto pix = curChipData->getData()[hit];
  uint16_t row = pix.getRowDirect(), col = pix.getCol();

  if (labelsClusPtr) { // MC labels were requested
    int nlab = 0;
    fetchMCLabels(curChipData->getStartID() + hit, labelsDigPtr, nlab);
    auto cnt = compClusPtr->size();
    for (int i = nlab; i--;) {
      labelsClusPtr->addElement(cnt, labelsBuff[i]);
    }
  }

  // add to compact clusters, which must be always filled
  unsigned char patt[ClusterPattern::MaxPatternBytes]{0x1 << (7 - (0 % 8))}; // unrolled 1 hit version of full loop in finishChip
  uint16_t pattID = (parent->mPattIdConverter.size() == 0) ? CompCluster::InvalidPatternID : parent->mPattIdConverter.findGroupID(1, 1, patt);
  if ((pattID == CompCluster::InvalidPatternID || parent->mPattIdConverter.isGroup(pattID)) && patternsPtr) {
    patternsPtr->emplace_back(1); // rowspan
    patternsPtr->emplace_back(1); // colspan
    patternsPtr->insert(patternsPtr->end(), std::begin(patt), std::begin(patt) + 1);
  }
  compClusPtr->emplace_back(row, col, pattID, curChipData->getChipID());
}

//__________________________________________________
Clusterer::Clusterer() : mPattIdConverter()
{
  /*gpu::empty_kernel<<<1, 1>>>();

  // Check for any errors
  cudaError_t cudaerr = cudaDeviceSynchronize();
  if (cudaerr != cudaSuccess)
    printf("Kernel launch failed with error \"%s\".\n", cudaGetErrorString(cudaerr));*/
#ifdef _PERFORM_TIMING_
  mTimer.Stop();
  mTimer.Reset();
  mTimerMerge.Stop();
  mTimerMerge.Reset();
#endif
}

//__________________________________________________
void Clusterer::ClustererThread::initChip(const ChipPixelData* curChipData, uint32_t first)
{
  // init chip with the 1st unmasked pixel (entry "from" in the mChipData)
  prev = column1 + 1;
  curr = column2 + 1;
  resetColumn(curr);

  pixels.clear();
  preClusterHeads.clear();
  preClusterIndices.clear();
  auto pix = curChipData->getData()[first];
  currCol = pix.getCol();
  curr[pix.getRowDirect()] = 0; // can use getRowDirect since the pixel is not masked
  // start the first pre-cluster
  preClusterHeads.push_back(0);
  preClusterIndices.push_back(0);
  pixels.emplace_back(-1, first); // id of current pixel
  noLeftCol = true;               // flag that there is no column on the left to check yet
}

//__________________________________________________
void Clusterer::ClustererThread::updateChip(const ChipPixelData* curChipData, uint32_t ip)
{
  const auto pix = curChipData->getData()[ip];
  uint16_t row = pix.getRowDirect(); // can use getRowDirect since the pixel is not masked
  if (currCol != pix.getCol()) {     // switch the buffers
    swapColumnBuffers();
    resetColumn(curr);
    noLeftCol = false;
    if (pix.getCol() > currCol + 1) {
      // no connection with previous column, this pixel cannot belong to any of the
      // existing preclusters, create a new precluster and flag to check only the row above for next pixels of this column
      currCol = pix.getCol();
      addNewPrecluster(ip, row);
      noLeftCol = true;
      return;
    }
    currCol = pix.getCol();
  }

  Bool_t orphan = true;

  if (noLeftCol) {                              // check only the row above
    if (curr[row - 1] >= 0) {
      expandPreCluster(ip, row, curr[row - 1]); // attach to the precluster of the previous row
      return;
    }
  } else {
#ifdef _ALLOW_DIAGONAL_ALPIDE_CLUSTERS_
    int neighbours[]{curr[row - 1], prev[row], prev[row + 1], prev[row - 1]};
#else
    int neighbours[]{curr[row - 1], prev[row]};
#endif
    for (auto pci : neighbours) {
      if (pci < 0) {
        continue;
      }
      if (orphan) {
        expandPreCluster(ip, row, pci); // attach to the adjascent precluster
        orphan = false;
        continue;
      }
      // reassign precluster index to smallest one
      if (preClusterIndices[pci] < preClusterIndices[curr[row]]) {
        preClusterIndices[curr[row]] = preClusterIndices[pci];
      } else {
        preClusterIndices[pci] = preClusterIndices[curr[row]];
      }
    }
  }
  if (orphan) {
    addNewPrecluster(ip, row); // start new precluster
  }
}

//__________________________________________________
void Clusterer::ClustererThread::fetchMCLabels(int digID, const ConstMCTruth* labelsDig, int& nfilled)
{
  // transfer MC labels to cluster
  if (nfilled >= MaxLabels) {
    return;
  }
  const auto& lbls = labelsDig->getLabels(digID);
  for (int i = lbls.size(); i--;) {
    int ic = nfilled;
    for (; ic--;) { // check if the label is already present
      if (labelsBuff[ic] == lbls[i]) {
        return;     // label is found, do nothing
      }
    }
    labelsBuff[nfilled++] = lbls[i];
    if (nfilled >= MaxLabels) {
      break;
    }
  }
  //
}

//__________________________________________________
void Clusterer::clear()
{
  // reset
#ifdef _PERFORM_TIMING_
  mTimer.Stop();
  mTimer.Reset();
  mTimerMerge.Stop();
  mTimerMerge.Reset();
#endif
}

//__________________________________________________
void Clusterer::print() const
{
  // print settings
  LOGP(info, "Clusterizer squashes overflow pixels separated by {} BC and <= {} in row/col seeking down to {} neighbour ROFs", mMaxBCSeparationToSquash, mMaxRowColDiffToMask, mSquashingDepth);
  LOG(info) << "Clusterizer masks overflow pixels separated by < " << mMaxBCSeparationToMask << " BC and <= "
            << mMaxRowColDiffToMask << " in row/col";

#ifdef _PERFORM_TIMING_
  auto& tmr = const_cast<TStopwatch&>(mTimer); // ugly but this is what root does internally
  auto& tmrm = const_cast<TStopwatch&>(mTimerMerge);
  LOG(info) << "Inclusive clusterization timing (w/o disk IO): Cpu: " << tmr.CpuTime()
            << " Real: " << tmr.RealTime() << " s in " << tmr.Counter() << " slots";
  LOG(info) << "Threads output merging timing                : Cpu: " << tmrm.CpuTime()
            << " Real: " << tmrm.RealTime() << " s in " << tmrm.Counter() << " slots";

#endif
}

//__________________________________________________
void Clusterer::reset()
{
  // reset for new run
  clear();
  mNHugeClus = 0;
}



ToyProblem::ToyProblem(std::unique_ptr<RegionExtractor> regionExtractor, std::unique_ptr<ClusterAlgorithm> clusterAlgorithm)
  : regionExtractor(std::move(regionExtractor)), clusterAlgorithm(std::move(clusterAlgorithm))
{
}

void ToyProblem::addChip(ChipPixelData* chipData)
{
  if (!regionExtractor) {
    throw std::runtime_error("Clustering algorithm is not set");
  }
  std::vector<std::vector<std::vector<int>>> chipRegions = regionExtractor->preprocess(chipData, MAX_DIST_X, MAX_DIST_Y);
  // std::cout << "Added chip " << extractedRegions.size() + 1 << " of " << extractionTasks.size() << std::endl;
  extractedRegions.insert(extractedRegions.end(), chipRegions.begin(), chipRegions.end());
  chipIds.push_back(chipData->getChipID());
}

void ToyProblem::addChipAsync(ChipPixelData* chipData)
{
  extractionTasks.push_back([this, chipData] { addChip(chipData); });
}

void ToyProblem::executeExtractionAsync()
{
  if (extractionTasks.empty())
    return;
  std::cout << "Number of tasks: " << extractionTasks.size() << std::endl;
  Timer timer("executeExtractionAsync");
  for (auto& task : extractionTasks) {
    task();
  }
}

void ToyProblem::executeToyProblem(std::vector<Clusterer::BBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  Timer timer("executeToyProblem");
  std::cout << "starting extraction" << std::endl;
  executeExtractionAsync();
  std::cout << "extraction done, starting clustering" << std::endl;
  performClustering(clusterBBoxes, clusterPixels);
  std::cout << "clustering done, starting postprocessing" << std::endl;
  postProcess();
}

void ToyProblem::performClustering(std::vector<Clusterer::BBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  if (!clusterAlgorithm) {
    throw std::runtime_error("Clustering stragegy is not set");
  }
  Timer timer("performClustering");
  clusterAlgorithm->clusterize(extractedRegions, chipIds, clusterBBoxes, clusterPixels);
}

void ToyProblem::postProcess()
{
  Timer timer("postProcess");
}


std::vector<std::vector<std::vector<int>>> ExpansionRegionExtractor::preprocess(const ChipPixelData* chipData, const int maxdist_x, const int maxdist_y)
{
  std::vector<std::vector<std::vector<int>>> extractedRegions;
  std::vector<PixelData> pixelData = chipData->getData();

  // could be replaced by unordered_set for better runtime complexity
  // requires hash function for PixelData, however
  std::set<PixelData> pixelSet(pixelData.begin(), pixelData.end());

  std::vector<std::vector<int>> fullRegion = convertSparsePixelsToGrid(pixelData);
  
  int numRows = fullRegion.size();
  int numCols = 0;
  if (!fullRegion.empty()) {
    numCols = fullRegion[0].size();
  }

  std::vector<std::vector<PixelData*>> pixelDataPointers(numRows, std::vector<PixelData*>(numCols, nullptr));
  for (const PixelData& pixel : pixelSet) {
    const PixelData* ptr = &pixel;
    pixelDataPointers[pixel.getRow()][pixel.getCol()] = const_cast<PixelData*>(ptr);
  }

  while (!pixelSet.empty()) {
    PixelData nextPixel = *pixelSet.begin();
    int nextRow = nextPixel.getRow();
    int nextCol = nextPixel.getCol();
    Region regionInfo = {nextRow, nextCol, 1, 1};

    while (expandRegion(fullRegion, regionInfo, maxdist_x, maxdist_y, pixelSet, pixelDataPointers)) { }

    std::vector<std::vector<int>> region(regionInfo.height, std::vector<int>(regionInfo.width, 0));

    for (int i = 0; i < regionInfo.height; ++i) {
      for (int j = 0; j < regionInfo.width; ++j) {
        region[i][j] = fullRegion[i + regionInfo.row][j + regionInfo.col];
      }
    }

    for (int row = regionInfo.row; row < regionInfo.row + regionInfo.height; ++row) {
      for (int col = regionInfo.col; col < regionInfo.col + regionInfo.width; ++col) {

        if (row >= fullRegion.size() || col >= fullRegion[row].size()) {
          std::cout << "ERROR: Trying to access out of bounds index in fullRegion" << std::endl;
        }

        if (fullRegion[row][col] != 1) continue;

        fullRegion[row][col] = 0;

        if (row >= pixelDataPointers.size() || col >= pixelDataPointers[row].size()) {
          std::cout << "ERROR: Trying to access out of bounds index in pixelDataPointers" << std::endl;
        }

        PixelData* pixelData = pixelDataPointers[row][col];
        if (pixelData) {
          if (pixelSet.find(*pixelData) == pixelSet.end()) {
            std::cout << "ERROR: Trying to erase an object from pixelSet that doesn't exist" << std::endl;
          }

          pixelSet.erase(*pixelData);
          pixelDataPointers[row][col] = nullptr;
        }
      }
    }
    extractedRegions.push_back(region);
  }

  return extractedRegions;
}

bool ExpansionRegionExtractor::expandRegion(std::vector<std::vector<int>>& fullRegion,
                                            Region& regionInfo,
                                            const int maxdist_x,
                                            const int maxdist_y,
                                            std::set<PixelData>& pixelSet,
                                            std::vector<std::vector<PixelData*>>& pixelDataPointers)
{
  int fullRegionHeight = fullRegion.size();
  int fullRegionWidth = fullRegion[0].size();

  int startRow = std::max(regionInfo.row - maxdist_x, 0);
  int startCol = std::max(regionInfo.col - maxdist_y, 0);
  int endRow = std::min(regionInfo.row + regionInfo.height + maxdist_x, fullRegionHeight);
  int endCol = std::min(regionInfo.col + regionInfo.width + maxdist_y, fullRegionWidth);

  bool regionExpanded = false;

  for (int row = startRow; row < endRow; ++row) {
    for (int col = startCol; col < endCol; ++col) {
      if (row >= regionInfo.row && row < regionInfo.row + regionInfo.height &&
          col >= regionInfo.col && col < regionInfo.col + regionInfo.width) {
        continue;
      }

      if (fullRegion[row][col] == 1) {
        if (col < regionInfo.col) {
          regionInfo.col = col;
          regionInfo.width += (regionInfo.col - col + 1);
        } else if (col >= regionInfo.col + regionInfo.width) {
          regionInfo.width += (col - regionInfo.col - regionInfo.width + 1);
        }

        if (row < regionInfo.row) {
          regionInfo.row = row;
          regionInfo.height += (regionInfo.row - row + 1);
        } else if (row >= regionInfo.row + regionInfo.height) {
          regionInfo.height += (row - regionInfo.row - regionInfo.height + 1);
        }
        regionExpanded = true;
      }
    }
  }
  return regionExpanded;
}

std::vector<std::vector<int>> ExpansionRegionExtractor::convertSparsePixelsToGrid(const std::vector<PixelData> pixels)
{
  uint16_t maxRow = 0, maxCol = 0;
  for (const auto& pixel : pixels) {
    maxRow = std::max(maxRow, pixel.getRow());
    maxCol = std::max(maxCol, pixel.getCol());
  }

  std::vector<std::vector<int>> grid(maxRow + 1, std::vector<int>(maxCol + 1, 0));

  for (const auto& pixel : pixels) {
    grid[pixel.getRow()][pixel.getCol()] = 1;
  }

  return grid;
}


#define CHECK_CUDA_ERROR(err)                                       \
  do {                                                              \
    if (err != cudaSuccess) {                                       \
      fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err)); \
    }                                                               \
  } while (0)

namespace o2::itsmft::gpu
{
using o2::itsmft::Point;
using o2::itsmft::BoundingBox;

__device__ void print(int* data, int* labels, int* parent, int* rank, int nrow, int ncol)
{
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      // Print elements of the Data Matrix
      printf("%d ", data[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Labels Matrix
      printf("%d ", labels[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Parent Matrix
      printf("%d ", parent[i * ncol + j]);
    }
    printf("\t\t");

    for (int j = 0; j < ncol; j++) {
      // Print elements of the Rank Matrix
      printf("%d ", rank[i * ncol + j]);
    }
    printf("\n");
  }
  printf("\n");
}

// recursively find root of cluster
__device__ int find(int x, int* parent)
{
  if (x != parent[x]) {
    parent[x] = find(parent[x], parent);
  }
  return parent[x];
}

// merge two clusters
__device__ int unify(int x, int y, int* parent, int* rank)
{
  int rootX = find(x, parent);
  int rootY = find(y, parent);

  if (rootX == rootY)
    return;

  if (rank[rootX] < rank[rootY]) {
    parent[rootX] = rootY;
    return rootY;
  } else {
    parent[rootY] = rootX;
    if (rank[rootX] == rank[rootY]) {
      rank[rootX]++;
    }
    return rootX;
  }
}

__global__ void ccl_kernel(int N, int* data, int* regionSizes, int* regionHeights, int* regionWidths,
                           int* startIndices, int* labels, int* parent, int* rank, int* numClusters,
                           int* clusterSizes, Point* scratchMemory, int* scratchMemIndex, int* stridedSizes, 
                           int* stridedPositions, int stride, BoundingBox* boxes, BoundingBox* stridedBoxes, 
                           int* chipIds, int* stridedChipIds)
{
  int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadIndex >= N) {
    return;
  }

  for (int regionIdx = threadIndex; regionIdx < N; regionIdx += blockDim.x * gridDim.x) {
    int startIndex = startIndices[regionIdx];
    int width = regionWidths[regionIdx];
    int height = regionHeights[regionIdx];

    int currentLabel = 1;

    // First pass
    for (int r = 0; r < height; ++r) {
      for (int c = 0; c < width; ++c) {
        int idx = startIndex + r * width + c;
        if (data[idx] == 1) {
          int leftIdx = (c == 0 ? -1 : idx - 1);
          int topIdx = (r == 0 ? -1 : idx - width);

          bool connectLeft = leftIdx != -1 && data[leftIdx] == 1;
          bool connectTop = topIdx != -1 && data[topIdx] == 1;

          if (!connectLeft && !connectTop) {
            labels[idx] = currentLabel++;
            parent[idx] = idx;
            rank[idx] = 0;
            numClusters[regionIdx]++;
            clusterSizes[idx] = 1;
            boxes[idx].min_r = r;
            boxes[idx].min_c = c;
            boxes[idx].max_r = r;
            boxes[idx].max_c = c;
          } else if (connectLeft && connectTop) {
            int topRootIdx = find(topIdx, parent);
            int leftRootIdx = find(leftIdx, parent);

            if (topRootIdx == leftRootIdx) {
              labels[idx] = labels[topRootIdx];
              parent[idx] = topRootIdx;
              rank[idx] = 0;

              // in this case, the bounding box already encompasses the newly added pixel
            } else {
              int newRoot = unify(leftRootIdx, topRootIdx, parent, rank);
              labels[idx] = labels[newRoot];
              parent[idx] = newRoot;
              numClusters[regionIdx]--;
              clusterSizes[newRoot] = clusterSizes[leftRootIdx] + clusterSizes[topRootIdx];

              boxes[newRoot].min_r = min(boxes[topRootIdx].min_r, boxes[leftRootIdx].min_r);
              boxes[newRoot].min_c = min(boxes[topRootIdx].min_c, boxes[leftRootIdx].min_c);
              boxes[newRoot].max_r = max(boxes[topRootIdx].max_r, boxes[leftRootIdx].max_r);
              boxes[newRoot].max_c = max(boxes[topRootIdx].max_c, boxes[leftRootIdx].max_c);
            }
          } else {
            int rootIdx;
            if (connectLeft) {
              rootIdx = find(leftIdx, parent);
            } else if (connectTop) {
              rootIdx = find(topIdx, parent);
            }
            labels[idx] = labels[rootIdx];
            parent[idx] = rootIdx;
            rank[idx] = 0;
            clusterSizes[rootIdx]++;

            boxes[rootIdx].min_r = min(boxes[rootIdx].min_r, r);
            boxes[rootIdx].min_c = min(boxes[rootIdx].min_c, c);
            boxes[rootIdx].max_r = max(boxes[rootIdx].max_r, r);
            boxes[rootIdx].max_c = max(boxes[rootIdx].max_c, c);
          }
        }
      }
    }

    int localClusterIndex = 0;

    // Second pass
    for (int r = 0; r < height; ++r) {
      for (int c = 0; c < width; ++c) {
        int idx = startIndex + r * width + c;
        if (data[idx] != 0) {
          int rootIdx = find(idx, parent);
          labels[idx] = labels[rootIdx];

          if (idx == rootIdx) {
            int currentClusterSize = clusterSizes[idx];
            BoundingBox currentClusterBox = boxes[idx];

            // global scratch memory index
            int reservedIndex = atomicAdd(scratchMemIndex, currentClusterSize);

            // pseudo local strided index
            int stridedIndex = stride * threadIndex + localClusterIndex;
            stridedSizes[stridedIndex] = currentClusterSize;
            stridedBoxes[stridedIndex] = currentClusterBox;
            stridedPositions[stridedIndex] = reservedIndex;
            stridedChipIds[stridedIndex] = chipIds[regionIdx];

            int localPixelIndex = 0;
            for (int r2 = 0; r2 < height; ++r2) {
              for (int c2 = 0; c2 < width; ++c2) {
                int idx2 = startIndex + r2 * width + c2;
                if (data[idx2] != 0 && find(idx2, parent) == rootIdx) {
                  scratchMemory[reservedIndex + localPixelIndex].r = r2;
                  scratchMemory[reservedIndex + localPixelIndex].c = c2;
                  localPixelIndex++;
                }
              }
            }
          }
          localClusterIndex++;
        }
      }
    }
  }
}
} // namespace o2::itsmft::gpu

std::vector<int> flatten(const std::vector<std::vector<std::vector<int>>>& data)
{
  std::vector<int> flatData;
  for (const auto& region : data)
    for (const auto& row : region)
      for (const auto& pixel : row)
        flatData.push_back(pixel);
  return flatData;
}

void PlayneClusterAlgorithm::clusterize(const std::vector<std::vector<std::vector<int>>>& data, const std::vector<int>& chipIds, std::vector<Clusterer::BBox>& clusterBBoxes, std::vector<std::vector<PixelData>>& clusterPixels)
{
  std::vector<int> flatData = flatten(data);

  int totalNumPixels = flatData.size();
  int numRegions = data.size();
  int scratchMemLength = 100 * numRegions; // assuming that there will not be more than 100 pixels in clusters per region on average
  int stride = 8;

  std::vector<int> regionSizes;
  std::vector<int> regionHeights;
  std::vector<int> regionWidths;
  std::vector<int> startIndices;
  int* labels = new int[totalNumPixels]();
  int* parent = new int[totalNumPixels]();
  int* rank = new int[totalNumPixels]();
  int* clusterSizes = new int[totalNumPixels]();
  int* numClusters = new int[numRegions]();

  Point* scratchMemory = new Point[scratchMemLength]();
  int* stridedSizes = new int[numRegions * stride]();
  int* stridedPositions = new int[numRegions * stride]();
  int scratchMemIndex = 0;
  BoundingBox* boxes = new BoundingBox[totalNumPixels]();
  BoundingBox* stridedBoxes = new BoundingBox[numRegions * stride]();
  int* stridedChipIds = new int[numRegions * stride]();

  // change that stuff to cudaMemset
  for (int i = 0; i < totalNumPixels; i++) {
    parent[i] = i;
    labels[i] = 0;
    rank[i] = 0;
    clusterSizes[i] = 0;

    boxes[i].min_r = 0;
    boxes[i].max_r = 0;
    boxes[i].min_c = 0;
    boxes[i].max_c = 0;
  }

  for (int i = 0; i < numRegions; i++) {
    numClusters[i] = 0;
  }

  for (int i = 0; i < numRegions * stride; i++) {
    stridedSizes[i] = 0;
    stridedPositions[i] = 0;
    stridedChipIds[i] = 0;
  }

  for (int i = 0; i < scratchMemLength; i++) {
    scratchMemory[i].r = -1;
    scratchMemory[i].c = -1;
  }

  for (int i = 0; i < numRegions * stride; i++) {
    stridedBoxes[i].min_r = 0;
    stridedBoxes[i].max_r = 0;
    stridedBoxes[i].min_c = 0;
    stridedBoxes[i].max_c = 0;
  }

  startIndices.push_back(0);

  for (const auto& region : data) {
    regionSizes.push_back(region.size() * (region.empty() ? 0 : region[0].size()));
    regionHeights.push_back(region.size());
    regionWidths.push_back((region.empty() ? 0 : region[0].size()));
    startIndices.push_back(startIndices.back() + regionSizes.back());
  }

  // Remove the last start index (beyond the end of the data)
  startIndices.pop_back();

  int* deviceData;
  int* deviceChipIds;
  int* deviceSizes;
  int* deviceHeights;
  int* deviceWidths;
  int* deviceStartIndices;
  int* deviceLabels;
  int* deviceParents;
  int* deviceRanks;
  int* deviceClusterSizes;
  int* deviceNumClusters;

  Point* deviceScratchMemory;
  int* deviceScratchMemIndex;
  int* deviceStridedSizes;
  int* deviceStridedPositions;
  int* deviceStridedChipIds;
  BoundingBox* deviceBoxes;
  BoundingBox* deviceStridedBoxes;

  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceData, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceChipIds, numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceSizes, regionSizes.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceHeights, regionHeights.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceWidths, regionWidths.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStartIndices, startIndices.size() * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceLabels, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceParents, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceRanks, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceClusterSizes, totalNumPixels * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceNumClusters, numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceScratchMemory, scratchMemLength * sizeof(Point)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceScratchMemIndex, sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedSizes, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedPositions, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedChipIds, stride * numRegions * sizeof(int)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceBoxes, totalNumPixels * sizeof(BoundingBox)));
  CHECK_CUDA_ERROR(cudaMalloc((void**)&deviceStridedBoxes, stride * numRegions * sizeof(BoundingBox)));

  CHECK_CUDA_ERROR(cudaMemcpy(deviceData, flatData.data(), totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceChipIds, chipIds.data(), numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceSizes, regionSizes.data(), regionSizes.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceWidths, regionWidths.data(), regionWidths.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceHeights, regionHeights.data(), regionHeights.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStartIndices, startIndices.data(), startIndices.size() * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceLabels, labels, totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceParents, parent, totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceRanks, rank, totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceClusterSizes, clusterSizes, totalNumPixels * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceNumClusters, numClusters, numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceScratchMemory, scratchMemory, scratchMemLength * sizeof(Point), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceScratchMemIndex, &scratchMemIndex, sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStridedSizes, stridedSizes, stride * numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStridedPositions, stridedPositions, stride * numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStridedChipIds, stridedChipIds, stride * numRegions * sizeof(int), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceBoxes, boxes, totalNumPixels * sizeof(BoundingBox), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(deviceStridedBoxes, stridedBoxes, stride * numRegions * sizeof(BoundingBox), cudaMemcpyHostToDevice));

  const int threadsPerBlock = 256;
  const int numBlocks = (numRegions + threadsPerBlock - 1) / threadsPerBlock;

  gpu::ccl_kernel<<<numBlocks, threadsPerBlock>>>(numRegions, deviceData, deviceSizes, deviceHeights,
                                                  deviceWidths, deviceStartIndices, deviceLabels, deviceParents, deviceRanks,
                                                  deviceNumClusters, deviceClusterSizes, deviceScratchMemory, deviceScratchMemIndex,
                                                  deviceStridedSizes, deviceStridedPositions, stride, deviceBoxes, deviceStridedBoxes, 
                                                  deviceChipIds, deviceStridedChipIds);

  CHECK_CUDA_ERROR(cudaMemcpy(labels, deviceLabels, totalNumPixels * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(scratchMemory, deviceScratchMemory, scratchMemLength * sizeof(int), cudaMemcpyDeviceToHost)); // copying can be shortened for tiny optimization
  CHECK_CUDA_ERROR(cudaMemcpy(stridedSizes, deviceStridedSizes, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedPositions, deviceStridedPositions, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedChipIds, deviceStridedChipIds, stride * numRegions * sizeof(int), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(stridedBoxes, deviceStridedBoxes, stride * numRegions * sizeof(BoundingBox), cudaMemcpyDeviceToHost));

  for (int clusterIdx = 0; clusterIdx < stride * numRegions; ++clusterIdx) {
    if (stridedSizes[clusterIdx] == 0)
      continue;

    clusterPixels.emplace_back();
    int clusterStartPosition = stridedPositions[clusterIdx];
    int chipId = stridedChipIds[clusterIdx];

    const BoundingBox& sBox = stridedBoxes[clusterIdx];
    Clusterer::BBox bBox(chipId);
    bBox.rowMin = static_cast<uint16_t>(sBox.min_r);
    bBox.colMin = static_cast<uint16_t>(sBox.min_c);
    bBox.rowMax = static_cast<uint16_t>(sBox.max_r);
    bBox.colMax = static_cast<uint16_t>(sBox.max_c);
    clusterBBoxes.push_back(bBox);

    for (int pixelIdx = 0; pixelIdx < stridedSizes[clusterIdx]; ++pixelIdx) {
      clusterPixels.back().push_back(PixelData(scratchMemory[clusterStartPosition + pixelIdx].r, scratchMemory[clusterStartPosition + pixelIdx].c));
    }
  }

  for (int i = 0; i < std::min(static_cast<int>(clusterPixels.size()), 20); i++) {
    for (const PixelData& pixel : clusterPixels[i]) {
      std::cout << "(" << pixel.getRow() << "," << pixel.getCol() << ") ";
    }
    std::cout << std::endl;
  }

  for (int i = 0; i < std::min(static_cast<int>(clusterBBoxes.size()), 20); i++) {
    std::cout << "(" << clusterBBoxes[i].chipID << "," << clusterBBoxes[i].rowMin << "," << clusterBBoxes[i].rowMax << "," << clusterBBoxes[i].colMin << "," << clusterBBoxes[i].colMax << ") "; 
    std::cout << std::endl;
  }

  CHECK_CUDA_ERROR(cudaFree(deviceData));
  CHECK_CUDA_ERROR(cudaFree(deviceSizes));
  CHECK_CUDA_ERROR(cudaFree(deviceWidths));
  CHECK_CUDA_ERROR(cudaFree(deviceHeights));
  CHECK_CUDA_ERROR(cudaFree(deviceStartIndices));
  CHECK_CUDA_ERROR(cudaFree(deviceLabels));
  CHECK_CUDA_ERROR(cudaFree(deviceParents));
  CHECK_CUDA_ERROR(cudaFree(deviceRanks));

  delete[] labels;
  delete[] parent;
  delete[] rank;
}