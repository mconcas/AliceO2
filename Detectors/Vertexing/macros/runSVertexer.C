/// \file runSVertexer.C
/// \brief Simple macro run dummy svertexer, it eventually will support GPU

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DetectorsVertexing/SVertexer.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include <TRandom.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <Math/SVector.h>
#include <array>

#include <dlfcn.h>

#endif

namespace o2
{
namespace vertexing
{

using Vec3D = ROOT::Math::SVector<double, 3>;

template <class FITTER>
float checkResults(o2::utils::TreeStreamRedirector& outs, std::string& treeName, const FITTER& fitter,
                   Vec3D& vgen, TLorentzVector& genPar, const std::vector<double>& dtMass)
{
  int nCand = fitter.getNCandidates();
  std::array<float, 3> p;
  float distMin = 1e9;
  for (int ic = 0; ic < nCand; ic++) {
    const auto& vtx = fitter.getPCACandidate(ic);
    auto df = vgen;
    df -= vtx;

    TLorentzVector moth, prong;
    for (int i = 0; i < fitter.getNProngs(); i++) {
      const auto& trc = fitter.getTrack(i, ic);
      trc.getPxPyPzGlo(p);
      prong.SetVectM({p[0], p[1], p[2]}, dtMass[i]);
      moth += prong;
    }
    auto nIter = fitter.getNIterations(ic);
    auto chi2 = fitter.getChi2AtPCACandidate(ic);
    double dst = TMath::Sqrt(df[0] * df[0] + df[1] * df[1] + df[2] * df[2]);
    distMin = dst < distMin ? dst : distMin;
    //    float genX
    outs << treeName.c_str() << "cand=" << ic << "ncand=" << nCand << "nIter=" << nIter << "chi2=" << chi2
         << "genPart=" << genPar << "recPart=" << moth
         << "genX=" << vgen[0] << "genY=" << vgen[1] << "genZ=" << vgen[2]
         << "dx=" << df[0] << "dy=" << df[1] << "dz=" << df[2] << "dst=" << dst << "\n";
  }
  return distMin;
}

TLorentzVector generate(Vec3D& vtx, std::vector<o2::track::TrackParCov>& vctr, float bz,
                        TGenPhaseSpace& genPHS, double parMass, const std::vector<double>& dtMass, std::vector<int> forceQ)
{
  const float errYZ = 1e-2, errSlp = 1e-3, errQPT = 2e-2;
  std::array<float, 15> covm = {
    errYZ * errYZ,
    0., errYZ * errYZ,
    0, 0., errSlp * errSlp,
    0., 0., 0., errSlp * errSlp,
    0., 0., 0., 0., errQPT * errQPT};
  bool accept = true;
  TLorentzVector parent, d0, d1, d2;
  do {
    accept = true;
    double y = gRandom->Rndm() - 0.5;
    double pt = 0.1 + gRandom->Rndm() * 3;
    double mt = TMath::Sqrt(parMass * parMass + pt * pt);
    double pz = mt * TMath::SinH(y);
    double phi = gRandom->Rndm() * TMath::Pi() * 2;
    double en = mt * TMath::CosH(y);
    double rdec = 10.; // radius of the decay
    vtx[0] = rdec * TMath::Cos(phi);
    vtx[1] = rdec * TMath::Sin(phi);
    vtx[2] = rdec * pz / pt;
    parent.SetPxPyPzE(pt * TMath::Cos(phi), pt * TMath::Sin(phi), pz, en);
    int nd = dtMass.size();
    genPHS.SetDecay(parent, nd, dtMass.data());
    genPHS.Generate();
    vctr.clear();
    float p[4];
    for (int i = 0; i < nd; i++) {
      auto* dt = genPHS.GetDecay(i);
      if (dt->Pt() < 0.05) {
        accept = false;
        break;
      }
      dt->GetXYZT(p);
      float s, c, x;
      std::array<float, 5> params;
      o2::math_utils::sincos(dt->Phi(), s, c);
      o2::math_utils::rotateZInv(vtx[0], vtx[1], x, params[0], s, c);

      params[1] = vtx[2];
      params[2] = 0.; // since alpha = phi
      params[3] = 1. / TMath::Tan(dt->Theta());
      params[4] = (i % 2 ? -1. : 1.) / dt->Pt();
      covm[14] = errQPT * errQPT * params[4] * params[4];
      //
      // randomize
      float r1, r2;
      gRandom->Rannor(r1, r2);
      params[0] += r1 * errYZ;
      params[1] += r2 * errYZ;
      gRandom->Rannor(r1, r2);
      params[2] += r1 * errSlp;
      params[3] += r2 * errSlp;
      params[4] *= gRandom->Gaus(1., errQPT);
      if (forceQ[i] == 0) {
        params[4] = 0.; // impose straight track
      }
      auto& trc = vctr.emplace_back(x, dt->Phi(), params, covm);
      float rad = forceQ[i] == 0 ? 600. : TMath::Abs(1. / trc.getCurvature(bz));
      if (!trc.propagateTo(trc.getX() + (gRandom->Rndm() - 0.5) * rad * 0.05, bz) ||
          !trc.rotate(trc.getAlpha() + (gRandom->Rndm() - 0.5) * 0.2)) {
        printf("Failed to randomize ");
        trc.print();
      }
    }
  } while (!accept);

  return parent;
}

} // namespace vertexing
} // namespace o2

void testCPU()
{
  int NTest = 10000;
  o2::utils::TreeStreamRedirector outStream("dcafitterNTestCPU.root");

  TGenPhaseSpace genPHS;
  double pion = 0.13957;
  double k0 = 0.49761;
  double kch = 0.49368;
  double dch = 1.86965;
  std::vector<double> k0dec = {pion, pion};
  std::vector<double> dchdec = {pion, kch, pion};
  std::vector<o2::track::TrackParCov> vctracks;
  o2::vertexing::Vec3D vtxGen;

  double bz = 5.0;
  // 2 prongs vertices

  LOG(INFO) << "Processing 2-prong Helix - Helix case";
  std::vector<int> forceQ{1, 1};

  o2::vertexing::DCAFitterN<2> ft; // 2 prong fitter
  ft.setBz(bz);
  ft.setPropagateToPCA(true);  // After finding the vertex, propagate tracks to the DCA. This is default anyway
  ft.setMaxR(200);             // do not consider V0 seeds with 2D circles crossing above this R. This is default anyway
  ft.setMaxDZIni(4);           // do not consider V0 seeds with tracks Z-distance exceeding this. This is default anyway
  ft.setMinParamChange(1e-3);  // stop iterations if max correction is below this value. This is default anyway
  ft.setMinRelChi2Change(0.9); // stop iterations if chi2 improves by less that this factor

  std::string treeName2A = "pr2a", treeName2W = "pr2w";
  TStopwatch swA, swW;
  int nfoundA = 0, nfoundW = 0;
  double meanDA = 0, meanDW = 0;
  swA.Stop();
  swW.Stop();
  for (int iev = 0; iev < NTest; iev++) {
    auto genParent = o2::vertexing::generate(vtxGen, vctracks, bz, genPHS, k0, k0dec, forceQ);

    ft.setUseAbsDCA(true);
    swA.Start(false);
    int ncA = ft.process(vctracks[0], vctracks[1]); // HERE WE FIT THE VERTICES
    swA.Stop();
    LOG(DEBUG) << "fit abs.dist " << iev << " NC: " << ncA << " Chi2: " << (ncA ? ft.getChi2AtPCACandidate(0) : -1);
    if (ncA) {
      auto minD = checkResults(outStream, treeName2A, ft, vtxGen, genParent, k0dec);
      meanDA += minD;
      nfoundA++;
    }

    ft.setUseAbsDCA(false);
    swW.Start(false);
    int ncW = ft.process(vctracks[0], vctracks[1]); // HERE WE FIT THE VERTICES
    swW.Stop();
    LOG(DEBUG) << "fit wgh.dist " << iev << " NC: " << ncW << " Chi2: " << (ncW ? ft.getChi2AtPCACandidate(0) : -1);
    if (ncW) {
      auto minD = checkResults(outStream, treeName2W, ft, vtxGen, genParent, k0dec);
      meanDW += minD;
      nfoundW++;
    }
  }
  ft.print();
  meanDA /= nfoundA ? nfoundA : 1;
  meanDW /= nfoundW ? nfoundW : 1;
  LOG(INFO) << "Processed " << NTest << " 2-prong vertices Helix : Helix";
  LOG(INFO) << "2-prongs with abs.dist minization: eff= " << float(nfoundA) / NTest
            << " mean.dist to truth: " << meanDA << " CPU time: " << swA.CpuTime();
  LOG(INFO) << "2-prongs with wgh.dist minization: eff= " << float(nfoundW) / NTest
            << " mean.dist to truth: " << meanDW << " CPU time: " << swW.CpuTime();
}

void runSVertexer()
{
  testCPU();

  const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
  o2::base::GeometryManager::loadGeometry();
  o2::base::Propagator::initFieldFromGRP(grp);
  o2::vertexing::SVertexerCUDA sV{};
  sV.init();
  // sV.process();

  return;
}
