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

/// \file Detector.h
/// \brief Definition of the Detector class

#ifndef ALICEO2_FCT_DETECTOR_H_
#define ALICEO2_FCT_DETECTOR_H_

#include <vector>                             // for vector
#include "DetectorsBase/GeometryManager.h"    // for getSensID
#include "DetectorsBase/Detector.h"           // for Detector
#include "DetectorsCommonDataFormats/DetID.h" // for Detector
#include "ITSMFTSimulation/Hit.h"             // for Hit
#include "Rtypes.h"                           // for Int_t, Double_t, Float_t, Bool_t, etc
#include "TArrayD.h"                          // for TArrayD
#include "TGeoManager.h"                      // for gGeoManager, TGeoManager (ptr only)
#include "TLorentzVector.h"                   // for TLorentzVector
#include "TVector3.h"                         // for TVector3
#include "FCTBase/FCTBaseParam.h"

class FairVolume;
class TGeoVolume;

class TParticle;

class TString;

namespace o2
{
namespace fct
{
class GeometryTGeo;
}
} // namespace o2
namespace o2
{
namespace fct
{
class FCTLayer;
}
} // namespace o2

namespace o2
{
namespace fct
{
class FCTLayer;

class Detector : public o2::base::DetImpl<Detector>
{
 public:
  /// Name : Detector Name
  /// Active: kTRUE for active detectors (ProcessHits() will be called)
  ///         kFALSE for inactive detectors
  Detector(Bool_t active);

  /// Default constructor
  Detector();

  /// Default destructor
  ~Detector() override;

  /// Initialization of the detector is done here
  void InitializeO2Detector() override;

  /// This method is called for each step during simulation (see FairMCApplication::Stepping())
  Bool_t ProcessHits(FairVolume* v = nullptr) override;

  /// Registers the produced collections in FAIRRootManager
  void Register() override;

  /// Gets the produced collections
  std::vector<o2::itsmft::Hit>* getHits(Int_t iColl) const
  {
    if (iColl == 0) {
      return mHits;
    }
    return nullptr;
  }

  /// Has to be called after each event to reset the containers
  void Reset() override;

  /// Base class to create the detector geometry
  void ConstructGeometry() override;

  /// This method is an example of how to add your own point of type Hit to the clones array
  o2::itsmft::Hit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                          const TVector3& startMom, double startE, double endTime, double eLoss,
                          unsigned char startStatus, unsigned char endStatus);

  Int_t chipVolUID(Int_t id) const { return o2::base::GeometryManager::getSensID(o2::detectors::DetID::FCT, id); }

  void EndOfEvent() override;

  void FinishPrimary() override { ; }
  virtual void finishRun() { ; }
  void BeginPrimary() override { ; }
  void PostTrack() override { ; }
  void PreTrack() override { ; }

  /// Returns the number of layers
  Int_t getNumberOfLayers() const { return mNumberOfLayers; }

  void buildBasicFCT(const FCTBaseParam& param);
  void buildFCTV1();
  void buildFCTFromFile(std::string);

  GeometryTGeo* mGeometryTGeo; //! access to geometry details

  void exportLayout();

 private:
  std::vector<Int_t> mLayerID;
  std::vector<Int_t> mConverterLayerId;
  std::vector<TString> mLayerName;
  std::vector<TString> mConverterLayerName;
  Int_t mNumberOfLayers;
  Int_t mNumberOfConverterLayers;
  Bool_t mOnlyChargedParticles;

  /// this is transient data about track passing the sensor
  struct TrackData {               // this is transient
    bool mHitStarted;              //! hit creation started
    unsigned char mTrkStatusStart; //! track status flag
    TLorentzVector mPositionStart; //! position at entrance
    TLorentzVector mMomentumStart; //! momentum
    double mEnergyLoss;            //! energy loss
  } mTrackData;                    //!

  /// Container for hit data
  std::vector<o2::itsmft::Hit>* mHits;

  /// Create the detector materials
  virtual void createMaterials();

  /// Create the detector geometry
  void createGeometry(const FCTBaseParam& param);

  /// Define the sensitive volumes of the geometry
  void defineSensitiveVolumes();

  Detector(const Detector&);

  Detector& operator=(const Detector&);

  std::vector<FCTLayer> mLayers;
  std::vector<FCTLayer> mConverterLayers;
  bool mIsPipeActivated = true; //! If Alice 3 pipe is present append inner disks to vacuum volume to avoid overlaps

  template <typename Det>
  friend class o2::base::DetImpl;
  ClassDefOverride(Detector, 1);
};

} // namespace fct
} // namespace o2

#ifdef USESHM
namespace o2
{
namespace base
{
template <>
struct UseShm<o2::fct::Detector> {
  static constexpr bool value = true;
};
} // namespace base
} // namespace o2
#endif

#endif
