#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <array>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>

#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITS/TrackITS.h"

#endif

using namespace std;
using namespace o2::gpu;
using namespace o2::itsmft;
using namespace o2::its;

//  using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
//  using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

struct DataFrames {
  void update(int frame, long index)
  {
    if (frame < firstFrame) {
      firstFrame = frame;
      firstIndex = index;
    }
    if (frame > lastFrame) {
      lastFrame = frame;
      lastIndex = index;
    }
    if (frame == firstFrame && index < firstIndex) {
      firstIndex = index;
    }
    if (frame == lastFrame && index > lastIndex) {
      lastIndex = index;
    }
  }

  long firstFrame = 10000;
  int firstIndex = 1;
  long lastFrame = -10000;
  int lastIndex = -1;
};

void CheckTrackerCA(std::string kineFile = "o2sim_Kine.root",
                    std::string clusFile = "o2clus_its.root",
                    std::string tracFile = "o2trac_its.root")
{
  auto f = TFile::Open("CheckTrackerCA.root", "recreate");
  auto nt = new TNtuple("ntt", "tracks ntuple", "rof:evid:isrec:isfake:chi2:mcPt:recPt:nClus:label");

  // Geometry
  o2::base::GeometryManager::loadGeometry();
  auto gman = o2::its::GeometryTGeo::Instance();

  // MC tracks
  auto fKin = TFile::Open(kineFile.data());
  auto mcTree = (TTree*)fKin->Get("o2sim");
  mcTree->SetBranchStatus("*", 0); // disable all branches
  mcTree->SetBranchStatus("MCTrack*", 1);
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree->SetBranchAddress("MCTrack", &mcArr);

  // Clusters
  auto fClu = TFile::Open(clusFile.data());
  auto clusTree = (TTree*)fClu->Get("o2sim");
  std::vector<CompClusterExt>* clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

  // Reconstructed tracks
  auto fRt = TFile::Open(tracFile.data());
  auto recTree = (TTree*)fRt->Get("o2sim");
  std::vector<TrackITS>* recArr = nullptr;
  recTree->SetBranchAddress("ITSTrack", &recArr);
  std::vector<o2::MCCompLabel>* trkLabArr = nullptr;
  recTree->SetBranchAddress("ITSTrackMCTruth", &trkLabArr);

  // Look for events information in frames
  int loadedEventClust{-1};
  int loadedEventTracks = -1;

  auto nEv = mcTree->GetEntriesFast();
  auto nRofClus = clusTree->GetEntriesFast();
  auto nRofRec = recTree->GetEntriesFast();
  vector<DataFrames> clusterFrames(nEv);
  vector<DataFrames> trackFrames(nEv);

  LOG(INFO) << "Find mc events in cluster frames.. ";
  for (auto iRof{0}; iRof < nRofClus; iRof++) { // Cluster frames
    if (!clusTree->GetEvent(iRof))
      continue;
    loadedEventClust = iRof;
    for (size_t iLast{0}; iLast < clusArr->size(); iLast++) { // Find the last MC event within this reconstructed entry
      auto lab = (clusLabArr->getLabels(iLast))[0];
      if (!lab.isValid() || lab.getSourceID() != 0 || lab.getEventID() < 0 || lab.getEventID() >= nEv)
        continue;
      clusterFrames[lab.getEventID()].update(iRof, iLast);
    }
  }

  LOG(INFO) << "Find mc events in track frames.. ";
  for (auto iRof{0}; iRof < nRofRec; iRof++) { // Track frames
    if (!recTree->GetEvent(iRof))
      continue;
    int loadedEventTracks = iRof;
    for (size_t iLast{0}; iLast < recArr->size(); iLast++) { // Find the last MC event within this reconstructed entry
      auto lab = (*trkLabArr)[iLast];
      if (!lab.isValid()) {
        const TrackITS& recTrack = (*recArr)[iLast];
      }
      if (!lab.isValid() || lab.getSourceID() != 0 || lab.getEventID() < 0 || lab.getEventID() >= nEv)
        continue;
      trackFrames[lab.getEventID()].update(iRof, iLast);
    }
  }
}