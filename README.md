# TIDAL

```
git clone --recurse-submodules git@github.com:Ksavva1021/TIDAL.git
```


## Setup Instructions

`Note: SVFit is a separate package so, install and load in a clean terminal`


### Install the environment:

```
micromamba env create -f environment.yml
```

### Setting up SVFIT - Careful:

```
source setup_svfit.sh
```

To check that SVFIT was set up properly. Run interactively in python3:

```
from TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper import FastMTT
x = FastMTT()
```

#### Using SVFIT:

Activate the environment and load ROOT+SVFIT: 

```
micromamba activate TIDAL
source load_package.sh # Select Option 2
```

For the time being SVFIT is only set up for the $\tau_h \tau_h$ channel. There are extra command-line arguments inside the `scripts\run_svfit.py` script but, you only need `--source_dir` and `--use_condor` 

Example Command - This will submit jobs to the batch (again reminder only for $\tau_h \tau_h$ channel at the moment, skips the rest): 
```
python3 scripts/run_svfit.py --source_dir /vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/ --use_condor
```

Output file is called `svfit.parquet` and contains `run, lumi, event, svfit mass, svfit error`. These are saved in the directory where the merged.parquet file is found e.g. 

```/vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/tt/DYto2L_M-50_madgraphMLM/nominal/svfit.parquet```


### Setting up Draw Utilities (Separate from SVFIT, clean terminal):

Setting up Draw utilities (Multidraw):

```
source setup_multidraw.sh
pip install -e .
```

#### Using Draw:

Activate the environment and load ROOT: 

```
micromamba activate TIDAL
source load_package.sh # Select Option 1
```
