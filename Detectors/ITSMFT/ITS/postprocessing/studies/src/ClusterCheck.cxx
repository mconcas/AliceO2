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

#include "ITSStudies/ClusterCheck.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/Task.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsParameters/GRPECSObject.h"

namespace o2
{
namespace its
{
namespace study
{
using namespace o2::framework;
using namespace o2::framework;
using namespace o2::globaltracking;

class ClusterUsageStudy : public Task
{
 public:
  ClusterUsageStudy(std::shared_ptr<DataRequest> dr,
                    std::shared_ptr<o2::base::GRPGeomRequest> gr,
                    bool useMC = false) : mDataRequest(dr),
                                          mGGCCDBRequest(gr),
                                          mUseMC(useMC) {}
  ~ClusterUsageStudy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void process(o2::globaltracking::RecoContainer&);

 private:
  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC{false};
};

void ClusterUsageStudy::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
}

void ClusterUsageStudy::run(ProcessingContext& pc)
{
}

void ClusterUsageStudy::endOfStream(EndOfStreamContext& ec)
{
}

void ClusterUsageStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
}

DataProcessorSpec getClusterUsageStudy(bool useMC)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  // ask for ITS clusters
  // ask for ITS tracks
  // ask for cluster usage reference
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              true,                              // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              true,                              // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);
  return DataProcessorSpec{
    "its-study-cluster-usage",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<ClusterUsageStudy>(dataRequest, ggRequest, useMC)},
    Options{}};
}
} // namespace study
} // namespace its
} // namespace o2