#!/bin/bash                                                                                                                                                                                                                           
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW CMSSW_8_0_20`
cd CMSSW_8_0_20/src/
eval `scram runtime -sh`
echo "CMSSW: "$CMSSW_BASE

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
let "sample=${1}+1"
xrdcp root://cmseos.fnal.gov//store/user/tmitchel/toCondor/CMSSW_8_0_20.tar.gz .
tar -xf CMSSW_8_0_20.tar.gz
rm CMSSW_8_0_20.tar.gz
cd CMSSW_8_0_20/src/Analysis/VLQAna
mkdir test
cd test
xrdcp root://cmseos.fnal.gov//store/user/tmitchel/toCondor/Zmumu/SAMPLE_${sample}.py .
cmsRun SAMPLE_${sample}.py
xrdcp Evt*.root root://cmseos.fnal.gov/OUTPUT/SAMPLE
rm SAMPLE_${sample}.py
rm -rf CMSSW_8_0_20
ls
echo "DONE!"
