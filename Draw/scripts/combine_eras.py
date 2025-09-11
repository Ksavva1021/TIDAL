import os
import shutil
import shlex
from Draw.python.PlotHistograms import HTT_Histogram

def find_sh_file(root_file):
    # Find the shell file with parameters used to make the plot
    job_dir = os.path.join(os.path.dirname(root_file), 'logs')
    log_file_name = root_file.split('/')[-1].split('datacard_')[1].split(f'_{channel}_Run3')[0] + ".sh"
    return os.path.join(job_dir, log_file_name)

def combine_histos(root_files, output_file):
    # Combine multiple ROOT files into one using hadd
    command = f"hadd -f {output_file} " + " ".join(root_files)
    print("\n>> COMBINING HISTOGRAMS:", command)
    os.system(command)

def plot_combined(cmb_file, tree_name, channel, era, variable, blind=False):
    # Plot the combined root files
    print(f"\n>> PLOTTING FILE: {cmb_file}")
    Histo_Plotter = HTT_Histogram(
        cmb_file,
        tree_name,
        channel,
        era,
        variable,
        blind=blind,
        log_y=False,
        is2Dunrolled=False,
    )
    Histo_Plotter.plot_1D_histo()

def combine_eras(eras_to_combine, channel, scheme, directory, f_name):
    # example args: [Run3_2022, Run3_2022EE], et, control, /path/to/dir, datacard_m_vis_cp_inclusive_tt_Run3_2022.root
    print(f"\n>> ERAS TO COMBINE: {eras_to_combine}")

    # Get the files
    files_to_combine = []
    old_args = {}
    for era in eras_to_combine:
        src_file = os.path.join(directory, era, scheme, channel, f_name.replace("$ERA", era))
        files_to_combine.append(src_file)
        print(f"\n>> FILE: {src_file}")
        log_file = find_sh_file(src_file)
        with open(log_file, 'r') as f:
            tokens = shlex.split(f.read())
        # Build a dictionary with all arguments
        args = {
            tokens[i].split("\n--")[1]: tokens[i + 1]
            for i in range(len(tokens) - 1)
            if tokens[i].startswith("\n--")
        }
        print("-> Arguments found")
        # Check that all parameters match in the arguments!
        if old_args:
            if args['channel'] == old_args['channel'] and args['category'] == old_args['category'] and args['var'] == old_args['var']:
                print("-> channels, category and variable match :)")
            else:
                raise ValueError("Different parameters found in the files to combine!")
        else:
            old_args = args

    # Combine histograms
    if ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine) and ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2022")
        cmb_file = os.path.join(cmb_directory, "datacard_m_vis_cp_inclusive_tt_full2223.root")
        era = "earlyrun3"
    elif ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2023")
        cmb_file = os.path.join(cmb_directory, "datacard_m_vis_cp_inclusive_tt_full22.root")
        era = "full22"
    elif ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_earlyRun3")
        cmb_file = os.path.join(cmb_directory, "datacard_m_vis_cp_inclusive_tt_full23.root")
        era = "full23"
    else:
        raise ValueError("Combination of eras not supported!")
    # make output directory
    os.makedirs(cmb_directory, exist_ok=True)

    # Check for extra name in node
    if 'nodename' in args.keys():
        tree_name = args['channel'] + '_' + args['category'] + args['nodename']
    else:
        tree_name = args['channel'] + '_' + args['category']

    # combine the histograms
    combine_histos(files_to_combine, cmb_file)

    # plot the combined histograms
    plot_combined(cmb_file, tree_name, channel, era, args['var'])




if __name__ == "__main__":

    directory = "/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/Plots/Control/Sep10_KlitosUnblinding"


    eras_to_combine = ['Run3_2023', 'Run3_2023BPix'] #, 'Run3_2022EE']
    # eras_to_combine = ['Run3_2022', 'Run3_2022EE'] #, 'Run3_2023', 'Run3_2023BPix']
    eras_to_combine = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']

    variable = 'm_vis'
    category = 'cp_inclusive'
    channel = "tt"

    f_name = f"datacard_{variable}_{category}_{channel}_$ERA.root"

    scheme = 'control'



    combine_eras(eras_to_combine, channel, scheme, directory, f_name)

