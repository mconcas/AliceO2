#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#endif

void CheckSquashedCluster(const string noSquashFileName = "o2clus_its.root", const string squashedFileName = "o2clus_its_squashed.root")
{
  auto fileNoSquashed = TFile::Open(noSquashFileName.data(), "read");
  auto fileSquashed = TFile::Open(squashedFileName.data(), "read");

  auto treeNoSq = (TTree*)fileNoSquashed->Get("o2sim");
  auto treeSq = (TTree*)fileSquashed->Get("o2sim");

  
}