#if !defined(__CLING__) || defined(__ROOTCLING__)
// o2sim->Draw("ITSTrack.mChi2/9 : ITSTrack.getPt() >> chiptcn(20,0,2,100,0,5)","ITSTrack.getNumberOfClusters()==7 && ITSTrack.getPt()<2 && !ITSTrackMCTruth.isFake()")
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include "DataFormatsITS/TrackITS.h"
#endif

using namespace std;

void inspectTracks(string trackFileName = "o2trac_its.root")
{
  auto f = TFile::Open(trackFileName.data());
  auto o2sim = (TTree *)f->Get("o2sim");
  auto c = new TCanvas("insTracks", "Track inspection");
  c->cd();
  o2sim->Draw("ITSTrack.mChi2/9 : ITSTrack.getPt() >> chiptcn(50,1.5,5,150,0,1000)", "ITSTrack.getNumberOfClusters()==7 && ITSTrack.getPt()>1.5&& ITSTrack.getPt()<5 && !ITSTrackMCTruth.isFake();#it{p}_{T} (GeV/#it{c});#chi^{2", "colz");
  auto p = new TCanvas("insTracksChi2Prof", "Track inspection, #chi^{2} profile");
  p->cd();
  p->SetLogy();
  o2sim->Draw("ITSTrack.mChi2/9 >> chiptcnprof(150,0,1000)", "ITSTrack.getNumberOfClusters()==7 && ITSTrack.getPt()>1.5&& ITSTrack.getPt()<5 && !ITSTrackMCTruth.isFake();#chi^{2");
}
