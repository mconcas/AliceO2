// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITStracking/Smoother.h"

namespace o2
{
namespace its
{

constexpr std::array<double, 3> getInverseSymm2D(const std::array<double, 3>& mat)
{
  const double det = mat[0] * mat[2] - mat[1] * mat[1];
  return std::array<double, 3>{mat[2] / det, -mat[1] / det, mat[0] / det};
}

// Smoother
template <unsigned int D>
Smoother<D>::Smoother(TrackITSExt& track, int smoothingLayer, const ROframe& event, float bZ, o2::base::PropagatorF::MatCorrType corr) : mLayerToSmooth{smoothingLayer},
                                                                                                                                         mBz(bZ),
                                                                                                                                         mCorr(corr)
{
  auto propInstance = o2::base::Propagator::Instance();
  const TrackingFrameInfo& smTf = event.getTrackingFrameInfoOnLayer(mLayerToSmooth).at(track.getClusterIndex(mLayerToSmooth));

  mInwardsTrack = track;
  mOutwardsTrack = {mInwardsTrack.getParamOut(),
                    static_cast<short>(mInwardsTrack.getNumberOfClusters()), -999, static_cast<std::uint32_t>(event.getROFrameId()),
                    mInwardsTrack.getParamOut(), mInwardsTrack.getClusterIndexes()};

  mInwardsTrack.resetCovariance();
  mInwardsTrack.setChi2(0);
  mOutwardsTrack.resetCovariance();
  mOutwardsTrack.setChi2(0);

  bool statusIw{false};
  bool statusOw{false};

  //////////////////////
  // Outward propagation
  for (size_t iLayer{0}; iLayer < mLayerToSmooth; ++iLayer) {
    if (mInwardsTrack.getClusterIndex(iLayer) == constants::its::UnusedIndex) { // Shorter tracks
      continue;
    }
    const TrackingFrameInfo& tF = event.getTrackingFrameInfoOnLayer(iLayer).at(mInwardsTrack.getClusterIndex(iLayer));
    statusIw = mInwardsTrack.rotate(tF.alphaTrackingFrame);
    statusIw &= propInstance->propagateToX(mInwardsTrack,
                                           tF.xTrackingFrame,
                                           mBz,
                                           o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                           o2::base::PropagatorImpl<float>::MAX_STEP,
                                           mCorr);
    mInwardsTrack.setChi2(mInwardsTrack.getChi2() + mInwardsTrack.getPredictedChi2(tF.positionTrackingFrame, tF.covarianceTrackingFrame));
    statusIw &= mInwardsTrack.o2::track::TrackParCov::update(tF.positionTrackingFrame, tF.covarianceTrackingFrame);
    LOG(INFO) << "Outwards loop on inwards track, layer: " << iLayer << " x: " << mInwardsTrack.getX();
  }

  // Prediction on the previously outwards-propagated track is done on a copy, as the process seems to be not reversible
  auto inwardsClone = mInwardsTrack;
  statusIw = inwardsClone.rotate(smTf.alphaTrackingFrame);
  statusIw &= propInstance->propagateToX(inwardsClone,
                                         smTf.xTrackingFrame,
                                         mBz,
                                         o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                         o2::base::PropagatorImpl<float>::MAX_STEP,
                                         mCorr);
  /////////////////////
  // Inward propagation
  for (size_t iLayer{D - 1}; iLayer > mLayerToSmooth; --iLayer) {
    if (mOutwardsTrack.getClusterIndex(iLayer) == constants::its::UnusedIndex) { // Shorter tracks
      continue;
    }
    const TrackingFrameInfo& tF = event.getTrackingFrameInfoOnLayer(iLayer).at(mOutwardsTrack.getClusterIndex(iLayer));
    statusOw = mOutwardsTrack.rotate(tF.alphaTrackingFrame);
    statusOw &= propInstance->propagateToX(mOutwardsTrack,
                                           tF.xTrackingFrame,
                                           mBz,
                                           o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                           o2::base::PropagatorImpl<float>::MAX_STEP,
                                           mCorr);
    mOutwardsTrack.setChi2(mOutwardsTrack.getChi2() + mOutwardsTrack.getPredictedChi2(tF.positionTrackingFrame, tF.covarianceTrackingFrame));
    statusOw &= mOutwardsTrack.o2::track::TrackParCov::update(tF.positionTrackingFrame, tF.covarianceTrackingFrame);
    LOG(INFO) << "Inwards loop on outwards track, layer: " << iLayer << " x: " << mOutwardsTrack.getX();
  }
  LOG(INFO) << "TF info, X: " << smTf.xTrackingFrame;
  // Prediction on the previously inwards-propagated track is done on a copy, as the process seems to be not revesible
  auto outwardsClone = mOutwardsTrack;
  statusOw = outwardsClone.rotate(smTf.alphaTrackingFrame);
  statusOw &= propInstance->propagateToX(outwardsClone,
                                         smTf.xTrackingFrame,
                                         mBz,
                                         o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                         o2::base::PropagatorImpl<float>::MAX_STEP,
                                         mCorr);
  // Compute weighted local chi2
  mBestChi2 = getSmoothedPredictedChi2(outwardsClone, inwardsClone, smTf.positionTrackingFrame, smTf.covarianceTrackingFrame);
}

template <unsigned int D>
float Smoother<D>::getSmoothedPredictedChi2(const o2::track::TrackParCov& outwTrack, // outwards track: from innermost cluster to outermost
                                            const o2::track::TrackParCov& inwTrack,  // inwards track: from outermost cluster to innermost
                                            const std::array<float, 2>& cls,
                                            const std::array<float, 3>& clCov)
{
  // Tracks need to be already propagated, compute only chi2
  // Symmetric covariances assumed

  if (outwTrack.getX() != inwTrack.getX()) {
    LOG(FATAL) << "Tracks need to be propagated to the same point! inwTrack.X=" << inwTrack.getX() << " outwTrack.X=" << outwTrack.getX();
  }

  std::array<double, 2> pp1 = {static_cast<double>(outwTrack.getY()), static_cast<double>(outwTrack.getZ())}; // P1: predicted Y,Z points
  std::array<double, 2> pp2 = {static_cast<double>(inwTrack.getY()), static_cast<double>(inwTrack.getZ())};   // P2: predicted Y,Z points

  std::array<double, 3> c1 = {static_cast<double>(outwTrack.getSigmaY2()),
                              static_cast<double>(outwTrack.getSigmaZY()),
                              static_cast<double>(outwTrack.getSigmaZ2())}; // Cov. track 1

  std::array<double, 3> c2 = {static_cast<double>(inwTrack.getSigmaY2()),
                              static_cast<double>(inwTrack.getSigmaZY()),
                              static_cast<double>(inwTrack.getSigmaZ2())}; // Cov. track 2

  std::array<double, 3> w1 = getInverseSymm2D(c1); // weight matrices
  std::array<double, 3> w2 = getInverseSymm2D(c2);

  std::array<double, 3> w1w2 = {w1[0] + w2[0], w1[1] + w2[1], w1[2] + w2[2]}; // (W1 + W2)
  std::array<double, 3> C = getInverseSymm2D(w1w2);                           // C = (W1+W2)^-1

  std::array<double, 2> w1pp1 = {w1[0] * pp1[0] + w1[1] * pp1[1], w1[1] * pp1[0] + w1[2] * pp1[1]}; // W1 * P1
  std::array<double, 2> w2pp2 = {w2[0] * pp2[0] + w2[1] * pp2[1], w2[1] * pp2[0] + w2[2] * pp2[1]}; // W2 * P2

  float Y = static_cast<float>(C[0] * (w1pp1[0] + w2pp2[0]) + C[1] * (w1pp1[1] + w2pp2[1])); // Pp: weighted normalized combination of the predictions:
  float Z = static_cast<float>(C[1] * (w1pp1[0] + w2pp2[0]) + C[2] * (w1pp1[1] + w2pp2[1])); // Pp = [(W1 * P1) + (W2 * P2)] / (W1 + W2)

  std::array<double, 2> delta = {Y - cls[0], Z - cls[1]};                                                                                         // Δ = Pp - X, X: space point of cluster (Y,Z)
  std::array<double, 3> CCp = {C[0] + static_cast<double>(clCov[0]), C[1] + static_cast<double>(clCov[1]), C[2] + static_cast<double>(clCov[2])}; // Transformation of cluster covmat: CCp = C * Cov
  std::array<double, 3> Wp = getInverseSymm2D(CCp);                                                                                               // Get weight matrix: Wp = CCp^-1

  float chi2 = static_cast<float>(delta[0] * (Wp[0] * delta[0] + Wp[1] * delta[1]) + delta[1] * (Wp[1] * delta[0] + Wp[2] * delta[1])); // chi2 = tΔ * (Wp * Δ)
#ifdef CA_DEBUG
  LOG(INFO) << "Propagated t1_y: " << pp1[0] << " t1_z: " << pp1[1];
  LOG(INFO) << "Propagated t2_y: " << pp2[0] << " t2_z: " << pp2[1];
  LOG(INFO) << "Smoothed prediction Y: " << Y << " Z: " << Z;
  LOG(INFO) << "cov t1: 0: " << c1[0] << " 1: " << c1[1] << " 2: " << c1[2];
  LOG(INFO) << "cov t2: 0: " << c2[0] << " 1: " << c2[1] << " 2: " << c2[2];
  LOG(INFO) << "cov Pr: 0: " << C[0] << " 1: " << C[1] << " 2: " << C[2];
  LOG(INFO) << "chi2: " << chi2;
  LOG(INFO) << "";
#endif
  return chi2;
}

template <unsigned int D>
bool Smoother<D>::testCluster(const int clusterId, const ROframe& event)
{
  auto propInstance = o2::base::Propagator::Instance();
  const TrackingFrameInfo& testTf = event.getTrackingFrameInfoOnLayer(mLayerToSmooth).at(mOutwardsTrack.getClusterIndex(mLayerToSmooth));

  bool statusIw{false};
  bool statusOw{false};

  // Prediction on the previously outwards-propagated track is done on a copy, as the process seems to be not reversible
  auto inwardsClone = mInwardsTrack;
  statusIw = inwardsClone.rotate(testTf.alphaTrackingFrame);
  statusIw &= propInstance->propagateToX(inwardsClone,
                                         testTf.xTrackingFrame,
                                         mBz,
                                         o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                         o2::base::PropagatorImpl<float>::MAX_STEP,
                                         mCorr);

  // Prediction on the previously inwards-propagated track is done on a copy, as the process seems to be not revesible
  auto outwardsClone = mOutwardsTrack;
  statusOw = outwardsClone.rotate(testTf.alphaTrackingFrame);
  statusOw &= propInstance->propagateToX(outwardsClone,
                                         testTf.xTrackingFrame,
                                         mBz,
                                         o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                         o2::base::PropagatorImpl<float>::MAX_STEP,
                                         mCorr);
  // Compute weighted local chi2
  float testedLocalChi2 = getSmoothedPredictedChi2(outwardsClone, inwardsClone, testTf.positionTrackingFrame, testTf.covarianceTrackingFrame);

  return testedLocalChi2 < mBestChi2;
}

template class Smoother<7>;

// bool Tracker::smoothTrack(TrackITSExt& track, const int testedClusterIndex, const int level, const ROframe& event)
// {
//   // This method applies on inward fitted tracks
//   // It is the decision-taking function based on a Kalman smoothing approach.
//   // If the substitution of a given cluster in a track provides a better local "smoothed chi2" the cluster is substituted.
//   // Selection is performed comparing the local chi2 obtained with propagation at the layer of the track models in two concurrent directions.
//   // Returns true if tested cluster provides better smoothed chi2; track definition is updated accordingly and refitted.
//   // Returns false otherwise; track is left untouched.
//   bool status;

//   // Need to copy the input track in case we change its clusters.
//   o2::its::TrackITSExt inwardsTrack = track;
//   o2::its::TrackITSExt outwardsTrack = {inwardsTrack.getParamOut(),
//                                         static_cast<short>(inwardsTrack.getNumberOfClusters()), -999, static_cast<std::uint32_t>(event.getROFrameId()),
//                                         inwardsTrack.getParamOut(), inwardsTrack.getClusterIndexes()}; // inwards track: from outermost cluster to innermost; // outwards track: from innermost cluster to outermost

//   outwardsTrack.resetCovariance();
//   outwardsTrack.setChi2(0);
//   inwardsTrack.resetCovariance();
//   inwardsTrack.setChi2(0);

//   // Seeding of tracks is done using their first two clusters respectively.
//   // They are propagated up to the cluster before the moothing level
//   // kalmanPropagateOutwardsTrack(event, outwardsTrack, 0, level);
//   // kalmanPropagateOutwardsTrack(event, inwardsTrack, inwardsTrack.getNumberOfClusters() - 1, inwardsTrack.getNumberOfClusters() - level);

//   std::array<TrackingFrameInfo, 2> choices; // put into array to reduce code repetition
//   choices[0] = event.getTrackingFrameInfoOnLayer(level).at(track.getClusterIndex(level));
//   choices[1] = event.getTrackingFrameInfoOnLayer(level).at(testedClusterIndex);

//   float bestChi2{o2::constants::math::VeryBig};
//   bool changed{false};
//   for (auto tf : choices) {
//     // perform first two steps of the fits but do not update the track
//     auto outwardsTrackCopy = outwardsTrack;
//     auto inwardsTrackCopy = inwardsTrack;
//     status = inwardsTrackCopy.rotate(tf.alphaTrackingFrame);
//     if (!status) {
//       LOG(INFO) << "Failed rotation inwards track in smoother";
//       return status;
//     }
//     status &= inwardsTrackCopy.propagateTo(tf.xTrackingFrame, getBz());
//     if (!status) {
//       LOG(INFO) << "Failed propagation inwards track in smoother";
//       return status;
//     }
//     status &= outwardsTrackCopy.rotate(tf.alphaTrackingFrame);
//     if (!status) {
//       LOG(INFO) << "Failed rotation outwards track in smoother";
//       return status;
//     }
//     status &= outwardsTrackCopy.propagateTo(tf.xTrackingFrame, getBz());
//     if (!status) {
//       LOG(INFO) << "Failed propagation outwards track in smoother";
//       return status;
//     }
//     float localChi2 = getSmoothedPredictedChi2(outwardsTrackCopy, inwardsTrackCopy, tf.positionTrackingFrame, tf.covarianceTrackingFrame);

//     if (localChi2 < bestChi2) {
//       if (bestChi2 != o2::constants::math::VeryBig) {
//         // this is only if tested cluster gave better local chi2 -> update track cluster list
//         outwardsTrack.setExternalClusterIndex(level, testedClusterIndex);
//         changed = true;
//       }
//       bestChi2 = localChi2;
//     }
//   }

//   // If track info changed: continue fitting process with new cluster
//   if (changed) {
//     float density = 2.33f;                       // Density of Si [g/cm^3]
//     float distance;                              // Default thickness
//     float xx0 = ((level > 2) ? 0.008f : 0.003f); // Thickness in units of X0
//     distance = xx0 * 9.37f;                      // Thickness in cm

//     // Finalize fitting step here
//     outwardsTrack.setChi2(outwardsTrack.getChi2() +
//                           outwardsTrack.getPredictedChi2(choices[1].positionTrackingFrame, choices[1].covarianceTrackingFrame));
//     status &= outwardsTrack.o2::track::TrackParCov::update(choices[1].positionTrackingFrame, choices[1].covarianceTrackingFrame);
//     // if (mMatLayerCylSet) {
//     //   const auto startingCluster = mPrimaryVertexContext->getClusters()[level][outwardsTrack.getClusterIndex(level)];
//     //   const auto endingCluster = mPrimaryVertexContext->getClusters()[level + 1][outwardsTrack.getClusterIndex(level + 1)]; // first/second in the direction of the kalman propagation

//     //   auto matbud = mMatLayerCylSet->getMatBudget(startingCluster.xCoordinate, startingCluster.yCoordinate, startingCluster.zCoordinate,
//     //                                               endingCluster.xCoordinate, endingCluster.yCoordinate, endingCluster.zCoordinate);
//     //   xx0 = matbud.meanX2X0;
//     //   density = matbud.meanRho;
//     //   distance = matbud.length;
//     // }

//     // if (!outwardsTrack.correctForMaterial(xx0, -distance * density, !mMatLayerCylSet)) { // (first < last) ? -1. : 1.)
//     //   return false;
//     // }
//     // Finalize kalman propagation and set new track in place of the old one.
//     // kalmanPropagateOutwardsTrack(event, outwardsTrack, level, outwardsTrack.getNumberOfClusters());
//     track = outwardsTrack;

//     return true;
//   }

//   return false;
// }
} // namespace its
} // namespace o2