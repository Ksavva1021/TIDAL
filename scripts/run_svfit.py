import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper import ClassicSVfit, TauTauHistogramAdapter, MeasuredTauLepton
import ROOT as root
from coffea.nanoevents.methods import vector
import awkward as ak
import numpy as np
import os
from alive_progress import alive_bar
import argparse


def calculate_p4(pt: ak.Array, eta: ak.Array, phi: ak.Array, mass: ak.Array) -> ak.Array:
    return ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "mass": mass
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    )


def is_valid_p4(p4: ak.Array) -> bool:
    is_nan_or_inf = np.isnan(p4.pt) | np.isnan(p4.eta) | np.isnan(p4.phi) | np.isnan(p4.mass) | np.isinf(p4.pt) | np.isinf(p4.eta) | np.isinf(p4.phi) | np.isinf(p4.mass)
    return not is_nan_or_inf


def process_event(lep1, lep2, mode, kappa, metx, mety, covMatrix,
                  lep1_MeasuredTauLepton, lep2_MeasuredTauLepton,
                  do_massConstraint, mass_constraint, verbosity):
    if is_valid_p4(lep1) and is_valid_p4(lep2):
        svFitAlgo = ClassicSVfit(verbosity)
        svFitAlgo.setHistogramAdapter(TauTauHistogramAdapter())
        svFitAlgo.addLogM_fixed(True, kappa)
        if do_massConstraint:
            svFitAlgo.setDiTauMassConstraint(mass_constraint)

        measuredTauLeptons = [lep1_MeasuredTauLepton, lep2_MeasuredTauLepton]

        svFitAlgo.integrate(measuredTauLeptons, metx, mety, covMatrix)

        isValidSolution = svFitAlgo.isValidSolution()

        if isValidSolution:
            return (svFitAlgo.getHistogramAdapter().getMass(), svFitAlgo.getHistogramAdapter().getMassErr())
        else:
            return (-1, -1)
    else:
        return (-1, -1)


def run_svfit(input_file: str, channel_: str):
    columns = [
        "met_pt", "met_phi", "met_covXX", "met_covXY", "met_covYY",
        "pt_1", "eta_1", "phi_1", "mass_1", "decayMode_1",
        "pt_2", "eta_2", "phi_2", "mass_2", "decayMode_2",
        "event", "run", "lumi"
    ]
    events = ak.from_parquet(input_file, columns=columns)
    nevents = len(events)

    channel = ak.Array([channel_] * nevents)
    dummy = ak.Array([0] * nevents)

    mode = ak.where(
        channel == "tt",
        ak.ones_like(dummy),
        ak.where(
            channel == "et",
            ak.ones_like(dummy) * 2,
            ak.where(
                channel == "mt",
                ak.ones_like(dummy) * 3,
                ak.zeros_like(dummy)
            )
        )
    )

    kappa = ak.where(
        channel == "em",
        ak.ones_like(dummy) * 3,
        ak.where(
            channel == "mt",
            ak.ones_like(dummy) * 4,
            ak.where(
                channel == "et",
                ak.ones_like(dummy) * 4,
                ak.where(
                    channel == "tt",
                    ak.ones_like(dummy) * 5,
                    ak.ones_like(dummy) * 3
                )
            )
        )
    )

    metx = events.met_pt * np.cos(events.met_phi)
    mety = events.met_pt * np.sin(events.met_phi)

    covXX = events.met_covXX
    covXY = events.met_covXY
    covYY = events.met_covYY

    covMET = np.zeros((len(covXX), 2, 2))
    covMET[:, 0, 0] = covXX
    covMET[:, 1, 0] = covXY
    covMET[:, 0, 1] = covXY
    covMET[:, 1, 1] = covYY

    rootMETMatrices = []

    for cov in covMET:
        matrix = root.TMatrixD(2, 2)
        matrix[0][0] = cov[0][0]
        matrix[1][0] = cov[1][0]
        matrix[0][1] = cov[0][1]
        matrix[1][1] = cov[1][1]
        rootMETMatrices.append(matrix)


    lep_1 = calculate_p4(events.pt_1, events.eta_1, events.phi_1, events.mass_1)
    lep_2 = calculate_p4(events.pt_2, events.eta_2, events.phi_2, events.mass_2)

    leps1_MeasuredTauLepton = np.array(
        [MeasuredTauLepton(*args) for args in zip(mode, lep_1.pt, lep_1.eta, lep_1.phi, lep_1.mass, events.decayMode_1)]
    )

    leps2_MeasuredTauLepton = np.array(
        [MeasuredTauLepton(*args) for args in zip(mode, lep_2.pt, lep_2.eta, lep_2.phi, lep_2.mass, events.decayMode_2)]
    )

    del events

    verbosity = 0
    do_massConstraint = False
    mass_constraint = -1

    svfitMass_values = [0] * nevents
    svfitMass_err_values = [0] * nevents

    with alive_bar(nevents) as bar:
        for index in range(nevents):
            try:
                mass, mass_err = process_event(
                    lep_1[index], lep_2[index], mode[index], kappa[index],
                    metx[index], mety[index], rootMETMatrices[index],
                    leps1_MeasuredTauLepton[index], leps2_MeasuredTauLepton[index],
                    do_massConstraint, mass_constraint, verbosity
                )
                svfitMass_values[index] = mass
                svfitMass_err_values[index] = mass_err
            except Exception as exc:
                print(f"Event {index} generated an exception: {exc}")
                svfitMass_values[index] = -1
                svfitMass_err_values[index] = -1

            bar()

    columns_to_store = ["event", "run", "lumi"]
    events = ak.from_parquet(input_file, columns=columns_to_store)
    events = ak.with_field(events, svfitMass_values, "svfitMass")
    events = ak.with_field(events, svfitMass_err_values, "svfitMass_err")

    output_file = input_file.replace("merged", "svfit")
    ak.to_parquet(events, output_file)


def run_processes(base_dir, use_condor=False):
    for channel in os.listdir(base_dir):
        channel_dir = os.path.join(base_dir, channel)
        if channel != "tt" :
            continue
        if os.path.isdir(channel_dir):
            for process in os.listdir(channel_dir):
                process_dir = os.path.join(channel_dir, process)
                if os.path.isdir(process_dir):
                    variation_dir = os.path.join(process_dir, "nominal")
                    if os.path.isdir(variation_dir):
                        parquet_file = os.path.join(variation_dir, 'merged.parquet')
                        if os.path.exists(parquet_file):
                            run_process(parquet_file, channel, use_condor)


def run_process(parquet_file, channel, use_condor):
    if not use_condor:
        run_svfit(parquet_file, channel)
    else:
        shell_script_path = create_shell_script(parquet_file, channel)
        submission_file_path = create_condor_submission_file(shell_script_path)
        submit_condor_job(submission_file_path)


def create_shell_script(parquet_file, channel):
    script_content = f"""
#!/bin/bash
LD_LIBRARY_PATH=/usr/local/lib
CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/

python3 scripts/run_svfit.py --input_file {parquet_file} --channel {channel} --running_dir
"""
    shell_script_path = os.path.join(os.path.dirname(parquet_file), 'svfit', 'svfit.sh')
    if not os.path.exists(os.path.dirname(shell_script_path)):
        os.makedirs(os.path.dirname(shell_script_path))
    with open(shell_script_path, 'w') as f:
        f.write(script_content)
    print(f"Shell script created at {shell_script_path}")
    os.system(f"chmod +x {shell_script_path}")
    return shell_script_path


def create_condor_submission_file(shell_script_path):
    directory = os.path.dirname(shell_script_path)
    submission_file_content = f"""
executable = {shell_script_path}
output = {directory}/svfit.out
error = {directory}/svfit.err
log = {directory}/svfit.log
request_memory = 8G
getenv = True
+MaxRuntime = 35000
queue
"""
    submission_file_path = os.path.join(directory, 'condor_submit.sub')
    with open(submission_file_path, 'w') as f:
        f.write(submission_file_content)

    return submission_file_path


def submit_condor_job(submission_file_path):
    os.system(f"condor_submit {submission_file_path}")


def main():
    parser = argparse.ArgumentParser(description="Run SVFit")
    parser.add_argument("--source_dir", type=str, help="Source directory with parquet files e.g. /vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/")
    parser.add_argument("--use_condor", action="store_true", help="Use condor submission")
    parser.add_argument("--input_file", type=str, help="Input file (used by condor submission)")
    parser.add_argument("--channel", type=str, help="Channel (used by condor submission)")
    parser.add_argument("--running_dir", action="store_true", help="Running directory (used by condor submission)")
    args = parser.parse_args()

    if not args.running_dir:
        run_processes(args.source_dir, args.use_condor)
    else:

        run_svfit(args.input_file, args.channel)


if __name__ == "__main__":
    main()
