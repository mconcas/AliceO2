#! /bin/bash -x

# # Cleanup
# shopt -s extglob
# rm -v !("CheckTrackerCA.C"|"run_checks.sh"|"CreateDictionaries.C"|"run_trac_ca_its.C"|"matbud.root"|"DisplayTrack.C"|"CheckTracks.C"|"CheckClusters.C")

# # Simulation
# # o2-sim-serial -m PIPE ITS -g boxgen --configKeyValues 'BoxGun.pdg=2212 ; BoxGun.eta[0]=0 ; BoxGun.eta[1]=0; BoxGun.number=2; BoxGun.prange[0]=0.5; BoxGun.prange[1]=0.5' -n 1
# o2-sim -m PIPE ITS -g pythia8 -n 1000

# # Digitization
# o2-sim-digitizer-workflow

# # Reconstruction, create clusters
# o2-its-reco-workflow --trackerCA

# # Dictionary production from clusters
# root -q CreateDictionaries.C++

# # Reconstruction, assign pattern id to cluster 
# o2-its-reco-workflow --trackerCA

# --tracking-mode async --trackerCA --configKeyValues 'ITSVertexerParam.phiCut=0.5; ITSVertexerParam.zCut=0.1'
[[ -f errorDumps.txt ]] && rm -v errorDumps.txt
root -q run_trac_ca_its.C++

# # Check clusters macro
# root -q CheckClusters.C++

# # Check read covariances diff
# git diff errorDumps.txt errorDumpsClustersCheck.txt

# # Processing
# root -q CheckTrackerCA.C++
root CheckTracks.C

