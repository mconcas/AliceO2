#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "ITSMFTSimulation/Hit.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TH2F.h>
#include <TH3F.h>

#include <vector>
#endif

void plotHits()
{
  TH2* ep = new TH2F("etaph", "hist_etaph;#varphi;#eta", 150, 0., TMath::TwoPi(), 150, -5, 5);
  TH2* zp = new TH2F("zph", "hist_zph;;#varphi;z ", 150, 0., TMath::TwoPi(), 300, -200, 200);
  TH2* zr = new TH2F("zr", "hist_zr;z;r ", 300, -350, 350, 300, 0, 100);
  TH3F* xyz = new TH3F("xyz", "hist_xyz;x;y;z", 300, -100, 100, 300, -100, 100, 300, -350, 350);

  std::vector<TFile*> hitFiles;
  hitFiles.push_back(TFile::Open("o2sim_HitsTF3.root"));
  hitFiles.push_back(TFile::Open("o2sim_HitsFT3.root"));
  hitFiles.push_back(TFile::Open("o2sim_HitsTRK.root"));

  TTree* trkTree = hitFiles[2] ? (TTree*)hitFiles[2]->Get("o2sim") : nullptr;
  TTree* ft3Tree = hitFiles[1] ? (TTree*)hitFiles[1]->Get("o2sim") : nullptr;
  TTree* tf3Tree = hitFiles[0] ? (TTree*)hitFiles[0]->Get("o2sim") : nullptr;

  // TRK
  std::vector<o2::itsmft::Hit>* trkHit = nullptr;
  trkTree->SetBranchAddress("TRKHit", &trkHit);

  // FT3
  std::vector<o2::itsmft::Hit>* ft3Hit = nullptr;
  ft3Tree->SetBranchAddress("FT3Hit", &ft3Hit);

  // TF3
  std::vector<o2::itsmft::Hit>* tf3Hit = nullptr;
  // tf3Tree->SetBranchAddress("TF3Hit", &tf3Hit);

  for (int iev = 0; iev < trkTree->GetEntries(); iev++) {
    trkTree->GetEntry(iev);
    for (const auto& h : *trkHit) {
      TVector3 posvec(h.GetX(), h.GetY(), h.GetZ());
      ep->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Eta());
      zp->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Z());
      zr->Fill(posvec.Z(), TMath::Hypot(posvec.X(), posvec.Y()));
      xyz->Fill(posvec.X(), posvec.Y(), posvec.Z());
    }
    ft3Tree->GetEntry(iev);
    for (const auto& h : *ft3Hit) {
      TVector3 posvec(h.GetX(), h.GetY(), h.GetZ());
      ep->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Eta());
      zp->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Z());
      zr->Fill(posvec.Z(), TMath::Hypot(posvec.X(), posvec.Y()));
      xyz->Fill(posvec.X(), posvec.Y(), posvec.Z());
    }
    // tf3Tree->GetEntry(iev);
    // for (const auto &h : *tf3Hit)
    // {
    //     TVector3 posvec(h.GetX(), h.GetY(), h.GetZ());
    //     ep->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Eta());
    //     zp->Fill(TVector2::Phi_0_2pi(posvec.Phi()), posvec.Z());
    //     zr->Fill (posvec.Z(),TMath::Hypot(posvec.X(), posvec.Y()));
    //     xyz->Fill(posvec.X(), posvec.Y(), posvec.Z());
    // }
  }
  auto* EPcanvas = new TCanvas("EtaPhi", "EP", 1000, 800);
  EPcanvas->cd();
  ep->Draw();
  EPcanvas->SaveAs("EtaPhi.png");
  auto* ZPcanvas = new TCanvas("ZPhi", "ZP", 1000, 800);
  ZPcanvas->cd();
  zp->Draw();
  ZPcanvas->SaveAs("ZPhi.png");
  auto* RZcanvas = new TCanvas("RZ", "RZ", 1000, 800);
  RZcanvas->cd();
  zr->Draw();
  RZcanvas->SaveAs("RZ.png");
  auto* XYZcanvas = new TCanvas("XYZ", "XYZ", 1000, 800);
  XYZcanvas->cd();
  xyz->Draw();
  XYZcanvas->SaveAs("XYZ.png");
}