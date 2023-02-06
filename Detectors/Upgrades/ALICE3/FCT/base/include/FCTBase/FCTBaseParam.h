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

#ifndef ALICEO2_FCT_BASEPARAM_H_
#define ALICEO2_FCT_BASEPARAM_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2
{
namespace fct
{

// **
// ** Parameters for FCT base configuration
// **

enum FCTGeometry {
  Default = 0,
  Telescope = 1
};

struct FCTBaseParam : public o2::conf::ConfigurableParamHelper<FCTBaseParam> {
  // Geometry Builder parameters

  Int_t geoModel = FCTGeometry::Default;

  // FCTGeometry::Telescope parameters
  Int_t nLayers = 10;
  Float_t z0 = -16.0;      // First layer z position
  Float_t zLength = 263.0; // Distance between first and last layers
  Float_t etaIn = 4.5;
  Float_t etaOut = 1.5;
  Float_t Layerx2X0 = 0.01;
  
  Bool_t OnlyChargedParticles = true; // When true, only charged partciles are saved. When false, all particles are saved.

  Int_t specialSetup = 0; // 0 when no special setup is active
  // 1: Windowed beam pipe on 90 degrees
  // 2: Windowed beam pipe on 45 degrees
  // 3: Vacuum Vessel + VacV Wall Optimistic
  // 4: Vacuum Vessel + VacV Wall Pesimistic
  // 5: Vacuum Vessel + VacV Wall Optimistic + Window 45 degrees
  // 5: Vacuum Vessel + VacV Wall Pesimistic + Window 45 degrees
  Bool_t drawFCT = false; // If true, draw fct on a canvas and save to .pdf
  Float_t winEtaMin = -3.5;
  Float_t winEtaMax = -5;

  // FCTGeometry::External file
  std::string configFile = ""; // Overrides geoModel parameter when provided

  O2ParamDef(FCTBaseParam, "FCTBase");
};

} // end namespace fct
} // end namespace o2

#endif // ALICEO2_FCT_BASEPARAM_H_
