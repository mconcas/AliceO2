#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "ITStracking/ROframe.h"
#include "ITStracking/IOUtils.h"
#endif

using namespace o2::gpu;
using o2::its::MemoryParameters;
using o2::its::TrackingParameters;

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

void CheckTrackerCA()
{
}