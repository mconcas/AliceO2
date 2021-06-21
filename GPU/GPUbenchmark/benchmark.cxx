// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file benchmark.cxx
/// \author: mconcas@cern.ch

#include "Shared/Kernels.h"

int main()
{
  // o2::benchmark::GPUbenchmark<char> bm_char{};
  // bm_char.run();
  o2::benchmark::GPUbenchmark<size_t> bm_size_t{};
  bm_size_t.run();
  // o2::benchmark::GPUbenchmark<int> bm_int{};
  // bm_int.run();

  return 0;
}