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


def create_condor_submit_file(logs_path: str, variable_name: str, submit_file: str, script_path: str):
    condor_template = f"""
executable = {script_path}
output = {logs_path}/condor_{variable_name}.out
error = {logs_path}/condor_{variable_name}.err
log = {logs_path}/condor_{variable_name}.log
request_memory = 8000
getenv = True
+MaxRuntime = 2000
queue
"""

    with open(submit_file, "w") as f:
        f.write(condor_template)
    os.system(f"chmod +x {submit_file}")

def create_shell_script(input_folder, output_folder, parameter_file, channel, era, method, category, variable, additional_selection, additional_weight, datacard_name, script_path, blind=False, auto_rebin=False):
    shell_script = f"""
#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
python3 Draw/scripts/HiggsTauTauPlot.py \\
--input_folder {input_folder} \\
--output_folder {output_folder} \\
--channel {channel} \\
--era {era} \\
--parameter_file {parameter_file} \\
--method {method} \\
--category {category} \\
--sel '{additional_selection}' \\
--add_weight '{additional_weight}' \\
--var {variable} \\
--datacard_name {datacard_name}"""

    if blind: shell_script+=' \\\n--blind'
    if auto_rebin: shell_script+=' \\\n--auto_rebin'

    with open(script_path, "w") as script_file:
        script_file.write(shell_script)
    os.system(f"chmod +x {script_path}")

parser = argparse.ArgumentParser()

parser.add_argument('--config', type=str, help='Configuration file')
parser.add_argument('--batch', action='store_true', help='Run in batch mode')

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

available_channels = ['ee','mm']
for channel in channels:
    if channel not in available_channels:
        raise ValueError(f"Channel {channel} is not a valid channel. Please choose from {available_channels}")

available_eras = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
for era in eras:
    if era not in available_eras:
        raise ValueError(f"Era {era} is not a valid era. Please choose from {available_eras}")

available_schemes = ['ip_calculation']
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
                auto_rebin = True

                variable = config['variables'][scheme]["definitions"][variable]
                variable = create_bins(variable)
                variable_name = variable.split('[')[0]

                parameter_file = f"{parameter_path}/{era}/params.yaml"

                if scheme == 'ip_calculation':
                    default_variable_name = variable_name
                    additional_selections = config['variables'][scheme]["additional_selections"][channel]
                    for selection_name, selection_string in additional_selections.items():
                        variable_name = f"{default_variable_name}_{selection_name}"
                        additional_selection = f"({selection_string})"
                        additional_weight = "(1)"

                        if args.batch:
                            logs = f"{output_folder}/logs"
                            subprocess.run(["mkdir", "-p", logs])
                            script_path = os.path.join(logs, f"{variable_name}_{category}.sh")
                            create_shell_script(input_folder, output_folder, parameter_file, channel, era, method, category, variable, additional_selection, additional_weight, variable_name, script_path, blind=blind, auto_rebin=auto_rebin)

                            submit_file = os.path.join(logs, f"submit_{variable_name}_{category}.sub")
                            create_condor_submit_file(logs, variable_name, submit_file, script_path)

                            subprocess.run(["condor_submit", submit_file])
                        else:
                            # Directly run the process
                            process = [
                                "python3", "Draw/scripts/HiggsTauTauPlot.py",
                                "--input_folder", input_folder,
                                "--output_folder", output_folder,
                                "--channel", channel,
                                "--era", era,
                                "--parameter_file", parameter_file,
                                "--method", method,
                                "--category", category,
                                "--sel", f'{additional_selection}',
                                "--add_weight", f'{additional_weight}',
                                "--var", variable,
                                "--datacard_name", variable_name,
                            #    "--do_ss"
                            ]
                            if blind: process.append("--blind")
                            print(" ".join(process))
                            subprocess.run(process)

                elif scheme == 'ip_control':
                    additional_selection = ""
                    additional_weight = "(1)"
                    if args.batch:
                        logs = f"{output_folder}/logs"
                        subprocess.run(["mkdir", "-p", logs])
                        script_path = os.path.join(logs, f"{variable_name}_{category}.sh")
                        create_shell_script(input_folder, output_folder, parameter_file, channel, era, method, category, variable, additional_selection, additional_weight, variable_name, script_path, blind=blind)

                        submit_file = os.path.join(logs, f"submit_{variable_name}_{category}.sub")
                        create_condor_submit_file(logs, variable_name, submit_file, script_path)
                        subprocess.run(["condor_submit", submit_file])

                    else:
                        # Directly run the process
                        process = [
                            "python3", "Draw/scripts/HiggsTauTauPlot.py",
                            "--input_folder", input_folder,
                            "--output_folder", output_folder,
                            "--channel", channel,
                            "--era", era,
                            "--parameter_file", parameter_file,
                            "--method", method,
                            "--category", category,
                            "--sel", f'{additional_selection}',
                            "--add_weight", f'{additional_weight}',
                            "--var", variable,
                            "--datacard_name", variable_name,
                        #    "--do_ss"
                        ]
                        if blind: process.append("--blind")
                        print(" ".join(process))
                        subprocess.run(process)
