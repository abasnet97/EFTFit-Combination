echo "Setting up your local area to run the TopEFT EFTFit module"

echo "Installing CMSSW"
export SCRAM_ARCH=slc7_amd64_gcc700
scram project CMSSW CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
scram b -j8

echo "Installing combine"
git clone git@github.com:cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/
git checkout v8.2.0
cd -
scram b -j8

echo "Installing the EFTFit-Combination"
cd $CMSSW_BASE/src/
git clone git@github.com:abasnet97/EFTFit-Combination.git
scram b -j8
cd -

echo "Installing the  CombineHarvester"
cd $CMSSW_BASE/src/
## Old commands for custom fork
#git clone git@github.com:TopEFT/CombineHarvester.git --branch crab_random
git clone git@github.com:cms-analysis/CombineHarvester.git
cd - # This one is just so we can use `cd -` again later
cd $CMSSW_BASE/src/CombineHarvester/CombineTools/python/combine/
cp $CMSSW_BASE/src/EFTFit-Combination/Fitter/test/crab_random.patch .
# Apply patch to CombineToolBase.py
git apply crab_random.patch
scram b -j8
cd -


echo "All done! Please see the READMEs for more details on how to run EFTFit"
