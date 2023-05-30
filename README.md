# EFTFit
This repository holds the custom files needed to run a EFT fit topcoffea datacards. This is primarily designed for use for combination work.

## New fancy install script
To quickly install this repo, simply run:<br>
`wget -O - https://raw.githubusercontent.com/TopEFT/EFTFit/master/install.sh | sh`<br>
NOTE: This will install the TopEFT custom CombineHarvester fork. If you need to use `-s -1` as implemented in combine, you'll need to install the main CombineHarvester repo.

### Setting up
 
  In order to run combine, you will need to get the appropriate CMSSW release and to clone several repositories.

#### Set up the CMSSW release
Install CMSSW_10_2_13
```
export SCRAM_ARCH=slc7_amd64_gcc700
scram project CMSSW CMSSW_10_2_13
cd CMSSW_10_2_13/src
scram b -j8
```

#### Get the Combine repository
Currently working with tag `v8.2.0`:

```
git clone git@github.com:cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/
git checkout v8.2.0
cd -
scram b -j8
```

#### Get the EFTFit repository
```
cd $CMSSW_BASE/src/
git clone git@github.com:abasnet97/EFTFit-Combination.git
scram b -j8
```

#### Get the CombineHarvester repository
This package is designed to be used with the CombineHarvester fork. This might cause errors when compiling, but you can safely ignore them.

```
git clone git@github.com:cms-analysis/CombineHarvester.git
cd CombineHarvester
git checkout 128e41eb
scram b -j8
```


### Fitting

Now we can actually run combine to perform the fits.

#### Running the fits
- Make sure you have done a `cmsenv` inside of `CMSSW_10_2_13/src/` (wherever you have it installed)
- Enter `CMSSW_10_2_13/src/EFTFit/Fitter/test`
- Copy all the relevant .txt and .root files for the two analyses. Make sure to copy `selectedWCs.txt` and `EFTParam_v8.npy` files as well. As a temporary workaround to make it easy to fetch these pertinent files, following commands can be executed: 
  ```
  xrdcp root://eosuser.cern.ch//eos/user/a/abasnet/EFT/combination_top21003_top22006/EFTParam_v8.npy . 
  xrdcp root://eosuser.cern.ch//eos/user/a/abasnet/EFT/combination_top21003_top22006/top22006_top21003_noSys_datacards .
  ```
- Run `combineCards.py` to merge them all into one txt file.
  ```
  cd top22006_top21003_noSys_datacards
  combineCards.py top21003=datacard_recoeft_run2_noSys.txt top22006=combinedcards.txt > top21003_top22006_combinedcard.txt
  cd -
  ```
- NOTE: combine uses a lot of recursive function calls to create the workspace. When running with systematics, this can cause a segmentation fault. You must run `ulimit -s unlimited` once per session to avoid this.
- Run the following command to generate the workspace file:
    ```
    text2workspace.py top22006_top21003_noSys_datacards/top21003_top22006_combinedcard.txt -o wps.root -P EFTFit-Combination.Fitter.AlternatePhysicsModel:alternatePhysicsModel --PO selectedWCs=top22006_top21003_noSys_datacards/selectedWCs.txt --PO fits=EFTParam_v8.npy --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints=cms
    ``` 
- Run combine with our EFTFit tools
  - Example:
    ```
    python -i ../scripts/EFTFitter.py
    fitter.batch1DScanEFT(basename='.081921.njet.ptbl.Float', batch='condor', workspace='wps.root', other=['-t', '-1'])
    ```
  - Once all jobs are finished, run the following (again inside `python -i ../scripts/EFTFitter.py`) to collect them in the `EFTFit/Fitter/fit_files` folder: 
    ```
    fitter.batchRetrieve1DScansEFT(basename='.081921.njet.ptbl.Float', batch='condor')
    ````

#### Plot making

To make simple 1D plots, use:
```
python -i ../scripts/EFTPlotter.py
plotter.BatchLLPlot1DEFT(basename_lst=['.081121.njet.16wc.Float'])
```
To make comparison plots (e.g. `njets` vs. `njets+ptbl`):
```
python -i ../scripts/EFTPlotter.py
plotter.BestScanPlot(basename_float_lst='.081721.njet.Float', basename_freeze_lst='.081821.njet.ptbl.Float', filename='_float_njet_ptbl', titles=['N_{jet} prof.', 'N_{jet}+p_{T}(b+l) prof.'], printFOM=True)
```

# Making impact plots
Impact plots must be done in three stages:
### Initial fit
Run 
```python
fitter.ImpactInitialFit(workspace='ptz-lj0pt_fullR2_anatest17_noAutostats_withSys.root', wcs=[])
```
to produce the initial fits. A blank `wcs` will run over all WCs.
### Nuisance fit
Run 
```python
fitter.ImpactNuisance(workspace='ptz-lj0pt_fullR2_anatest17_noAutostats_withSys.root', wcs=[])
```
to fit each NP. A blank `wcs` will run over all WCs.
### Produce plots
Run 
```python
fitter.ImpactCollect(workspace='ptz-lj0pt_fullR2_anatest17_noAutostats_withSys.root', wcs=[])
```
to collect all jobs and create the final pdf plots. A blank `wcs` will run over all WCs.
