# #! /bin/bash

# Cleanup
shopt -s extglob
rm -v !("CheckTrackerCA.C"|"run.sh"|"CreateDictionaries.C"|"run_trac_ca_its.C"|"matbud.root"|"DisplayTrack.C"|"CheckTracks.C")

# Simulation
o2-sim-serial -m PIPE ITS -g boxgen --configKeyValues 'BoxGun.pdg=2212 ; BoxGun.eta[0]=0 ; BoxGun.eta[1]=0; BoxGun.number=2; BoxGun.prange[0]=0.5; BoxGun.prange[1]=0.5' -n 1

# Digitization
o2-sim-digitizer-workflow

# Reconstruction
o2-its-reco-workflow --trackerCA

# Dictionary production
root -q CreateDictionaries.C++

# --tracking-mode async --trackerCA --configKeyValues 'ITSVertexerParam.phiCut=0.5; ITSVertexerParam.zCut=0.1'
root -q run_trac_ca_its.C++

# Processing
root -q CheckTrackerCA.C++