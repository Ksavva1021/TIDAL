# TIDAL

# Setup Instructions (Please follow in order to avoid issues :) )

Setting up Draw utilities (Multidraw):

`source setup_multidraw.sh`


Setting up SVFIT - Careful:

`source setup_svfit.sh`

To check that SVFIT was set up properly. Run interactively in python3:

```
from TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper import FastMTT
x = FastMTT()
If this goes through, you are fine
```

Install the environment:

`micromamba env create -f environment.yml`

# Every-time usage instructions (Important):

Activate the environment: `micromamba activate TIDAL`

Load ROOT for Draw: `source load_package.sh` $\rightarrow$ Select 1

Load ROOT + SVFIT for SVFIT: `source load_package.sh` $\rightarrow$ Select 2

# Running SVFIT:
For the time being SVFIT is only set up for the $\tau_h \tau_h$ channel. There are extra command-line arguments inside the `scripts\run_svfit.py` script but, you only need `--source_dir` and `--use_condor` 

Example Command - This will submit jobs to the batch (again reminder only for $\tau_h \tau_h$ channel at the moment, skips the rest): 
```
python3 scripts/run_svfit.py --source_dir /vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/ --use_condor
```

Output file is called `svfit.parquet` and contains `run, lumi, event, svfit mass, svfit error`. These are saved in the directory where the merged.parquet file is found e.g. 

```/vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/tt/DYto2L_M-50_madgraphMLM/nominal/svfit.parquet```








