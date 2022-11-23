#include "ITSMFTSimulation/AlpideSimResponse.h"
#include <TFile.h>
#include <TSystem.h>
#include <cstdio>
#include <cstddef>
#include <fstream>
#include <iostream>

void alpideResponse(const string path = "./")
{

  o2::itsmft::AlpideSimResponse resp0, resp1;
  std::string responseFile = "AlpideResponseData.root";

  resp0.initData(0, path.data());
  resp1.initData(1, path.data());

  auto file = TFile::Open(responseFile.data(), "recreate");
  file->WriteObjectAny(&resp0, "o2::itsmft::AlpideSimResponse", "response0");
  file->WriteObjectAny(&resp1, "o2::itsmft::AlpideSimResponse", "response1");
  file->Close();
}