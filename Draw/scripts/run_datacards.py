from glob import glob
import os
import json
import yaml
import argparse
import subprocess



def get_args():
    parser = argparse.ArgumentParser(description="Run datacard creation for CP analysis")
    parser.add_argument('--vsjet', type=str, help="vsjet cut medium, tight, vtight", required=True)
    parser.add_argument('--IPcut', type=str, help="IP cut", required=True)
    parser.add_argument('--Ecut', type=str, help="Ecut", required=True)
    parser.add_argument('--name', type=str, help="name of the trial", required=True)
    parser.add_argument('--batch', action='store_true', help="Run in batch mode (datacards only)")
    parser.add_argument('--stack', action='store_true', help="Run datacard stacking (if datacards done in batch)")
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

def run_scan(vsjet, name, IPcut, Ecut, batch, stack):
    config_path = f"Draw/scripts/cpdecay_{vsjet}.yaml"
    with open(config_path, 'r') as file:
        cfg = yaml.safe_load(file)
    # path to use for stacking input
    datacard_path = f"{cfg['output_path']}/{vsjet}/{name}/*/*/*/*.root"
    stacked_path = f"{cfg['output_path']}/{vsjet}/{name}/added_histo_{vsjet}.root"
    # convert vsjet string to number
    if vsjet == "medium":
        n_vsjet = 5
    elif vsjet == "tight":
        n_vsjet = 6
    elif vsjet == "vtight":
        n_vsjet = 7
    # IF RUNNING ON BATCH ONLY CREATE DATACARDS
    if batch:
        run_command(f"python Draw/scripts/makeDatacards_cpdecay.py --config={config_path} --vsjet={n_vsjet} --IPcut={IPcut} --Ecut={Ecut} --name={vsjet}/{name} --batch")
    # IF RAN ON BATCH THEN STACK AS SEPARATE STEP
    elif stack:
        run_command(f"python Draw/scripts/hadd_cp_datacards.py -i {datacard_path} -o {stacked_path}")
    # IF RUNNING LOCAL THEN CREATE DATACARDS AND STACK
    else:
        run_command(f"python Draw/scripts/makeDatacards_cpdecay.py --config={config_path} --vsjet={n_vsjet} --IPcut={IPcut} --Ecut={Ecut} --name={vsjet}/{name}")
        run_command(f"python Draw/scripts/hadd_cp_datacards.py -i {datacard_path} -o {stacked_path}")


def main():
    args = get_args()
    try:
        run_scan(args.vsjet, args.name, args.IPcut, args.Ecut, args.batch, args.stack)
    except subprocess.CalledProcessError:
        raise RuntimeError("Production failed!")

if __name__ == "__main__":
    main()