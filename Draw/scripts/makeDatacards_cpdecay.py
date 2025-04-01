import argparse
import yaml
import numpy
import subprocess
import os

def create_bins(variable: str) -> str:
    if '(' in variable:
        variable_name = variable.split('(')[0]
        number_of_bins = variable.split('(')[1].split(',')[0]
        lower_bound = variable.split('(')[1].split(',')[1].strip()
        upper_bound = variable.split('(')[1].split(',')[2].split(')')[0].strip()

        # Function to handle pi expressions
        def evaluate_bound(bound: str) -> float:
            if "pi" in bound:
                bound = bound.replace("pi", "*numpy.pi")
                if bound.startswith('*'):
                    bound = bound[1:]  # Remove leading '*' for valid syntax
                if bound.startswith('-*'):
                    bound = '-' + bound[2:]  # Fix '-*' syntax
            return eval(bound)  # Safely evaluate

        # Evaluate bounds
        lower_bound = str(evaluate_bound(lower_bound))
        upper_bound = str(evaluate_bound(upper_bound))

        # Create the binning array
        binning = numpy.linspace(float(lower_bound), float(upper_bound), int(number_of_bins) + 1)
        binning = ','.join([str(i) for i in binning])
        variable = f"{variable_name}[{binning}]"

    return variable

def create_condor_submit_file(logs_path: str, variable_name: str, cat_name: str, submit_file: str, script_path: str):
    condor_template = f"""
executable = {script_path}
output = {logs_path}/condor_{variable_name}_{cat_name}.out
error = {logs_path}/condor_{variable_name}_{cat_name}.err
log = {logs_path}/condor_{variable_name}_{cat_name}.log
request_memory = 8000
getenv = True
+MaxRuntime = 10800
queue
"""

    with open(submit_file, "w") as f:
        f.write(condor_template)
    os.system(f"chmod +x {submit_file}")

def create_shell_script(input_folder, output_folder, parameter_file, channel, era, method, category, variable, script_path, blind=False, LO_DY=False, do_ss=False):
    if args.TES_variation != 'nominal' and 'm_vis' in variable:
        variable = variable.split('[')[0] + f"_{args.TES_variation}[" + variable.split('[')[1]
        print("Changing variable to TES variation")
    shell_script = f"""
#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
python3 Draw/scripts/HiggsTauTauPlot.py \\
--input_folder {input_folder} \\
--output_folder {output_folder} \\
--parameter_file {parameter_file} \\
--channel {channel} \\
--era {era} \\
--method {method} \\
--category {category} \\
--var {variable}"""

    if blind: shell_script+=' \\\n--blind'

    if LO_DY: shell_script+=' \\\n--LO_DY'

    if do_ss: shell_script+=' \\\n--do_ss'

    with open(script_path, "w") as script_file:
        script_file.write(shell_script)
    os.system(f"chmod +x {script_path}")

parser = argparse.ArgumentParser()

parser.add_argument('--config', type=str, help='Configuration file')
parser.add_argument('--batch', action='store_true', help='Run in batch mode')
parser.add_argument('--LO_DY', action='store_true', help='Use LO instead of NLO DY')
parser.add_argument('--do_ss', action='store_true', help='Same sign')

args = parser.parse_args()

config_file = args.config

# read config file which is a yaml file
with open(config_file) as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

input_folder = config['input_folder']
output_path = config['output_path']
channels = config['channels']
eras = config['eras']
parameter_path = config['parameter_path']
schemes = config['schemes']
variables = config['variables']

available_channels = ['ee','mm', 'et', 'mt', 'tt']
for channel in channels:
    if channel not in available_channels:
        raise ValueError(f"Channel {channel} is not a valid channel. Please choose from {available_channels}")

available_eras = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
for era in eras:
    if era not in available_eras:
        raise ValueError(f"Era {era} is not a valid era. Please choose from {available_eras}")

available_schemes = ['cpdecay', 'control']
for scheme in schemes:
    if scheme not in available_schemes:
        raise ValueError(f"Scheme {scheme} is not a valid scheme. Please choose from {available_schemes}")


for era in eras:
    for channel in channels:
        for scheme in schemes:
            settings = variables[scheme][channel]
            output_folder = f"{output_path}/{era}/{scheme}/{channel}"
            subprocess.run(["mkdir", "-p", output_folder])
            for setting in settings:
                method = setting[0]
                category = setting[1]
                variable = setting[2]

                blind = len(setting)>=4 and setting[3]

                variable = config['variables'][scheme]["definitions"][variable]
                variable = create_bins(variable)
                variable_name = variable.split('[')[0]

                parameter_file = f"{parameter_path}/{era}/params.yaml"

                if args.batch:
                    logs = f"{output_folder}/logs"
                    subprocess.run(["mkdir", "-p", logs])
                    script_path = os.path.join(logs, f"{variable_name}_{category}.sh")
                    create_shell_script(input_folder, output_folder, parameter_file, channel, era, method, category, variable, script_path, blind=blind, LO_DY=args.LO_DY, do_ss=args.do_ss)

                    submit_file = os.path.join(logs, f"submit_{variable_name}_{category}.sub")
                    create_condor_submit_file(logs, variable_name, category, submit_file, script_path)

                    subprocess.run(["condor_submit", submit_file])
                else:
                    # Directly run the process
                    process = [
                        "python3", "Draw/scripts/HiggsTauTauPlot.py",
                        "--input_folder", input_folder,
                        "--output_folder", output_folder,
                        "--parameter_file", parameter_file,
                        "--channel", channel,
                        "--era", era,
                        "--method", method,
                        "--category", category,
                        "--var", variable,
                    #    "--do_ss"
                    ]
                    if blind: process.append("--blind")
                    if args.LO_DY: process.append("--LO_DY")
                    if args.do_ss: process.append("--do_ss")
                    print(" ".join(process))
                    subprocess.run(process)
