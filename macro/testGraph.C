#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "ITStracking/Algorithms.h"
#include <vector>

void testGraph() {
  std::vector<float> data = {1.,2.,3.,4.,5.,6.,7.,8.,9.};
  // std::function<bool(int&,int&) fun;
  bool c = true;
  o2::its::Graph<float> grafo{data, [c](const float& a, const float& b){ return a < b && c ;}};
  
}
#endif

// [](int& v1, int& v2, Args&&) { return v1 < v2; }