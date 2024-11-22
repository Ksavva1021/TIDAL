from glob import glob
import os
import json
import yaml
import argparse
import subprocess

### COMMANDS TO RUN:
# ```
# python Draw/scripts/makeDatacards_cpdecay.py --config=Draw/scripts/cpdecay_medium.yaml --vsjet=5 --IPcut=1.25 --Ecut=0.0 --name=IP_1p25_E_0p0
# -> SPECIFY vsjet time/number and name
# ```

# ```
# python Draw/scripts/hadd_cp_datacards.py -i Draw/Run3_Nov20_medium/IP_1p25_E_0p0/*/*/*/*.root -o Draw/Run3_Nov20_medium/IP_1p25_E_0p0/added_Nov20_medium.root
# -> SPECIFY OUTPUT FOLDER (from NAME + VSJET)
# ```


def get_args():
    parser = argparse.ArgumentParser(description="Run combine workflow for CP analysis")
    parser.add_argument('--vsjet', type=str, help="vsjet cut medium, tight, vtight", required=True)
    parser.add_argument('--IPcut', type=str, help="IP cut", required=True)
    parser.add_argument('--Ecut', type=str, help="Ecut", required=True)
    parser.add_argument('--name', type=str, help="name of the trial", required=True)
    return parser.parse_args()

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
        print("\n")
        print("*"*140)
        print(f"\033[1;32m\nCommand '{command}' ran successfully!\n\033[0m")
        print("*"*140)
        print("\n")
    except subprocess.CalledProcessError as e:
        print("\n")
        print("*"*140)
        print(f"\033[1;31m\nCommand '{command}' failed with exit code: {e.returncode}\n\033[0m")
        print("*"*140)
        print("\n")
        raise

def run_scan(vsjet, name, IPcut, Ecut):
    config_path = f"Draw/scripts/cpdecay_{vsjet}.yaml"
    datacard_path = f"Draw/Run3_Nov20_{vsjet}/{name}/*/*/*/*.root"
    stacked_path = f"Draw/Run3_Nov20_{vsjet}/{name}/added_Nov20_{vsjet}.root"
    if vsjet == "medium":
        n_vsjet = 5
    elif vsjet == "tight":
        n_vsjet = 6
    elif vsjet == "vtight":
        n_vsjet = 7
    # CREATE DATACARDS
    run_command(f"python Draw/scripts/makeDatacards_cpdecay.py --config={config_path} --vsjet={n_vsjet} --IPcut={IPcut} --Ecut={Ecut} --name={name}")
    # ADD HISTOGRAMS ACROSS CHANNELS
    run_command(f"python Draw/scripts/hadd_cp_datacards.py -i {datacard_path} -o {stacked_path}")


def main():
    args = get_args()
    try:
        run_scan(args.vsjet, args.name, args.IPcut, args.Ecut)
    except subprocess.CalledProcessError:
        raise RuntimeError("Production failed!")

if __name__ == "__main__":
    main()