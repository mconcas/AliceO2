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

/// \file Ray.cxx
/// \brief Implementation of ray between start-end points for material budget estimate

#include "DetectorsBase/Ray.h"
#include "GPUCommonMath.h"

using namespace o2::base;
using namespace o2::gpu;

//______________________________________________________
GPUd() int Ray::crossLayer(const MatLayerCyl& lr)
{
  // Calculate parameters t of intersection with cyl.layer
  // Calculated as solution of equation for ray crossing with circles of r (rmin and rmax)
  // t^2*mDistXY2 +- sqrt( mXDxPlusYDy^2 - mDistXY2*(mR02 - r^2) )
  // Region of valid t is 0:1.
  // Straigh line may have 2 crossings with cyl. layer
  float detMax = mXDxPlusYDy2 - mDistXY2 * (mR02 - lr.getRMax2());
  if (detMax < 0) {
    return 0; // does not reach outer R, hence inner also
  }
  float detMaxRed = CAMath::Sqrt(detMax) * mDistXY2i;
  float tCross0Max = mXDxPlusYDyRed + detMaxRed; // largest possible t
  float tmpsq = CAMath::Sqrt(detMax);
  float tmpRMax2 = lr.getRMax2();
// #ifdef __CUDACC__
//   if (!(blockIdx.x * blockDim.x + threadIdx.x))
//     printf("crossLayer debug: mXDxPlusYDy2 = %x mDistXY2 = %x mR02 = %x lr.getRMax2() = %x detMax = %x mDistXY2i = %x sqDetMax = %x detMaxRed = %x \n",
//            __float_as_uint(mXDxPlusYDy2),
//            __float_as_uint(mDistXY2),
//            __float_as_uint(mR02),
//            __float_as_uint(tmpRMax2),
//            __float_as_uint(detMax),
//            __float_as_uint(mDistXY2i),
//            __float_as_uint(tmpsq),
//            __float_as_uint(detMaxRed));
// #else
//   printf("crossLayer debug: mXDxPlusYDy2 = %x mDistXY2 = %x mR02 = %x lr.getRMax2() = %x detMax = %x mDistXY2i = %x sqDetMax = %x detMaxRed = %x \n",
//          reinterpret_cast<unsigned int&>(mXDxPlusYDy2),
//          reinterpret_cast<unsigned int&>(mDistXY2),
//          reinterpret_cast<unsigned int&>(mR02),
//          reinterpret_cast<unsigned int&>(tmpRMax2),
//          reinterpret_cast<unsigned int&>(detMax),
//          reinterpret_cast<unsigned int&>(mDistXY2i),
//          reinterpret_cast<unsigned int&>(tmpsq),
//          reinterpret_cast<unsigned int&>(detMaxRed));
// #endif
  if (tCross0Max < 0) { // max t is outside of the limiting point -> other t's also
    return 0;
  }

  float tCross0Min = mXDxPlusYDyRed - detMaxRed; // smallest possible t
  if (tCross0Min > 1.f) {                        // min t is outside of the limiting point -> other t's also
    return 0;
  }
  float detMin = mXDxPlusYDy2 - mDistXY2 * (mR02 - lr.getRMin2());
  if (detMin < 0) { // does not reach inner R -> just 1 tangential crossing
    mCrossParams1[0] = tCross0Min > 0.f ? tCross0Min : 0.f;
    mCrossParams2[0] = tCross0Max < 1.f ? tCross0Max : 1.f;
    return validateZRange(mCrossParams1[0], mCrossParams2[0], lr);
  }
  int nCross = 0;
  float detMinRed = CAMath::Sqrt(detMin) * mDistXY2i;
  float tCross1Max = mXDxPlusYDyRed + detMinRed;
  float tCross1Min = mXDxPlusYDyRed - detMinRed;

  if (tCross1Max < 1.f) {
    mCrossParams1[0] = tCross0Max < 1.f ? tCross0Max : 1.f;
    mCrossParams2[0] = tCross1Max > 0.f ? tCross1Max : 0.f;
    if (validateZRange(mCrossParams1[nCross], mCrossParams2[nCross], lr)) {
      nCross++;
    }
  }

  if (tCross1Min > -0.f) {
    mCrossParams1[nCross] = tCross1Min < 1.f ? tCross1Min : 1.f;
    mCrossParams2[nCross] = tCross0Min > 0.f ? tCross0Min : 0.f;
    if (validateZRange(mCrossParams1[nCross], mCrossParams2[nCross], lr)) {
      nCross++;
    }
  }
  return nCross;
}
