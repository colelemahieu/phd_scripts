#!/bin/bash

dirR2="HiForward_Hiforest_PbPb_Run2018_ak2/231111_150349"
dirR4="HiForward_Hiforest_PbPb_Run2018_ak4/231104_201920"
dirR6="HiForward_Hiforest_PbPb_Run2018_ak6/231110_223836"


if [ ${1} -lt 10 ] 
then
	path="/eos/user/c/clemahie/HIForward/${dirR4}/000${1}"
else
	path="/eos/user/c/clemahie/HIForward/${dirR4}/00${1}"
fi


outPath="/eos/cms/store/group/phys_heavyions/clemahie/miniAOD_r4_${1}.root"

if test -d "${path}"; then
    root -b -l -q "/afs/cern.ch/user/c/clemahie/DijetAnalysis/miniAOD_thesis_skim.C(\"${path}/HiForestMiniAOD_*.root\", \"${outPath}\")" 

fi





