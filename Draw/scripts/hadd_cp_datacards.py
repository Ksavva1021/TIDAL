import ROOT
import argparse
import os
from Draw.python.PlotHistograms import HTT_Histogram

def hadd_root_files(input_files, output_file, dir_combinations, exp_num=None):
    """
    Combine ROOT files and merge specific directories by adding their histograms.
    
    Args:
        input_files (list of str): List of input ROOT files to combine.
        output_file (str): Name of the output ROOT file.
        dir_combinations (dict): Dictionary where keys are new directory names and values are lists
                                 of directories to combine.
                                 E.g., {'higgs_pirho': ['tt_higgs_pirho', 'tt_higgs_rhopi']}
    """
    # Create an empty output ROOT file
    output = ROOT.TFile(output_file, "RECREATE")
    
    # First, create a map of (hist_name, dir_name, file_name) to histogram
    histogram_dir_file_map = {}

    for file_name in input_files:

        input_root = ROOT.TFile.Open(file_name, "READ")
        
        # Loop over all the keys in the input file
        for key in input_root.GetListOfKeys():

            obj = key.ReadObj()

            # Check if the object is a directory
            if isinstance(obj, ROOT.TDirectoryFile):
                dir_name = obj.GetName()

                # Loop over all histograms in the directory
                input_root.cd(dir_name)
                for hist_key in ROOT.gDirectory.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    hist = hist_key.ReadObj()
                        
                    # Check if it's a histogram
                    if isinstance(hist, ROOT.TH1):
                        out_key = (hist_name,dir_name,file_name)
                        histogram_dir_file_map[out_key] = hist.Clone()
                        histogram_dir_file_map[out_key].SetDirectory(0)
       
    # now we convert this map into a histogram_dir_map with format (hist_name, dir_name) : hist
    # we check that a histogram exists for all files, if it is a systematic histogram identifed by ending in 'Up', 'Up_2D', 'Down', or 'Down_2D' then it might not exist in all files
    # in this case we will add the nominal histogram from that file instead

    # first create empty histogram_dir_map by looping over the keys in histogram_dir_file_map
    # we also get a set of the nominal histograms names which are these ones not ending in 'Up', 'Up_2D', 'Down', or 'Down_2D'
    histogram_dir_map = {}
    nominal_hist_names = set()

    for (hist_name, dir_name, file_name) in histogram_dir_file_map.keys():
        out_key = (hist_name, dir_name)
        if out_key not in histogram_dir_map:
            histogram_dir_map[out_key] = None
        if not (hist_name.endswith('Up') or hist_name.endswith('Down') or hist_name.endswith('Up_2D') or hist_name.endswith('Down_2D')):
            nominal_hist_names.add(hist_name)

    for key in histogram_dir_map.keys():
        hist_name, dir_name = key
        hist = None
        hist_count = 0
        for file_name in input_files:
            file_key = (hist_name, dir_name, file_name)
        
            if file_key in histogram_dir_file_map:
                if hist is None:
                    hist = histogram_dir_file_map[file_key].Clone()
                    hist.SetDirectory(0)
                    hist_count += 1
                else:
                    hist.Add(histogram_dir_file_map[file_key])
                    hist_count += 1
            else:
                # check if the histogram is a systematic histogram
                if hist_name.endswith('Up') or hist_name.endswith('Down') or hist_name.endswith('Up_2D') or hist_name.endswith('Down_2D'):
                    # get the nominal histogram name, loop over the nominal_hist_names and find the one that matches
                    nominal_hist_name = None
                    for nh in nominal_hist_names:
                        if hist_name.startswith(nh.replace('_2D','')) and (hist_name.endswith('_2D') == nh.endswith('_2D')):
                            nominal_hist_name = nh
                            break

                    if nominal_hist_name is not None:
                        # add the nominal histogram from this file instead

                        nominal_file_key = (nominal_hist_name, dir_name, file_name)
                        if nominal_file_key in histogram_dir_file_map:
                            if hist is None:
                                hist = histogram_dir_file_map[nominal_file_key].Clone()
                                hist.SetDirectory(0)
                                hist_count += 1
                            else:
                                hist.Add(histogram_dir_file_map[nominal_file_key])
                                hist_count += 1
                    
        if hist is None:
            print(f"Warning: Histogram {hist_name} in directory {dir_name} not found in any input files. Not adding to output.")
        elif exp_num is None or (hist_count == exp_num):
            histogram_dir_map[key] = hist
        elif hist_count != exp_num:
            print(f"Warning: Histogram {hist_name} in directory {dir_name} found in {hist_count} out of expected {exp_num} input files. Not adding to output.")


    # now we take care of combining directories together if needed based on dir_combinations
    if dir_combinations:
        combined_histogram_dir_map = {}
        for (hist_name, dir_name), hist in histogram_dir_map.items():
            if hist is None: continue
            combined = False
            for combined_dir, dirs_to_combine in dir_combinations.items():
                if dir_name in dirs_to_combine:
                    combined = True
                    out_key = (hist_name, combined_dir)
                    if out_key not in combined_histogram_dir_map:
                        combined_histogram_dir_map[out_key] = hist.Clone()
                        combined_histogram_dir_map[out_key].SetDirectory(0)
                    else:
                        combined_histogram_dir_map[out_key].Add(hist)
            if not combined:
                out_key = (hist_name, dir_name)
                if out_key not in combined_histogram_dir_map:
                    combined_histogram_dir_map[out_key] = hist.Clone()
                    combined_histogram_dir_map[out_key].SetDirectory(0)
                else:
                    combined_histogram_dir_map[out_key].Add(hist)
        histogram_dir_map = combined_histogram_dir_map

    dir_names = []

    # now we write the histograms to the output file
    for (hist_name, dir_name), hist in histogram_dir_map.items():
        if hist is None: continue
        # Create the directory in the output file if it doesn't exist
        if not output.GetDirectory(dir_name):
            output.mkdir(dir_name)
            dir_names.append(dir_name)
        output.cd(dir_name)
        hist.SetName(hist_name)
        hist.Write(hist_name)

    # Close the output file
    output.Close()

    for dir_name in dir_names:

        blind = True
        if 'mva_fake' in dir_name or 'mva_tau' in dir_name or 'aiso' in dir_name or '_ss' in dir_name:
            blind = False

        var_name = "Bin number"
        if 'BDT_score' in dir_name or 'mva_fake' in dir_name or 'mva_tau' in dir_name or 'mva_higgs' in dir_name:
            var_name = "BDT score"
        elif 'aiso' in dir_name or dir_name == 'tt_higgs_pipi_ss':
            var_name = r"$\phi_{CP}$"
    
        # make a plot of the combined histograms
        Histo_Plotter = HTT_Histogram(
            output_file,
            dir_name,
            'tt',
            '...',
            var_name,
            blind=blind,
            log_y=False,
            is2Dunrolled=False,
            save_name=output_file.replace('.root', f'_{dir_name}')
        )
        Histo_Plotter.plot_1D_histo()



if __name__ == "__main__":
    # Argument parser for command-line arguments
    parser = argparse.ArgumentParser()
    
    # Add input and output files arguments
    parser.add_argument('-i', '--input', nargs='+', required=True, help="List of input ROOT files")
    parser.add_argument('-o', '--output', required=True, help="Name of the output ROOT file")
    parser.add_argument('-e', '--expected', type=int, default=None, help="Expected number of input files each histogram should appear in i.e the number of eras")
    
    # Parse arguments
    args = parser.parse_args()

    # List of input ROOT files from command-line
    input_files = args.input

    # Output ROOT file from command-line
    output_file = args.output

    exp_num = args.expected

    # Define which directories to combine (customize this based on your case)
    dir_combinations = {
        'tt_higgs_pirho': ['tt_higgs_pirho', 'tt_higgs_rhopi'],
        'tt_higgs_rhoa1': ['tt_higgs_a1rho', 'tt_higgs_rhoa1'],
        'tt_higgs_pia1': ['tt_higgs_a1pi', 'tt_higgs_pia1'],
        'tt_higgs_pia11pr': ['tt_higgs_pia11pr', 'tt_higgs_a11prpi'],
        'tt_higgs_a11pra1': ['tt_higgs_a11pra1', 'tt_higgs_a1a11pr'],

        'tt_fake_pirho': ['tt_fake_pirho', 'tt_fake_rhopi'],
        'tt_fake_rhoa1': ['tt_fake_a1rho', 'tt_fake_rhoa1'],
        'tt_fake_pia1': ['tt_fake_a1pi', 'tt_fake_pia1'],
        'tt_fake_pia11pr': ['tt_fake_pia11pr', 'tt_fake_a11prpi'],
        'tt_fake_a11pra1': ['tt_fake_a11pra1', 'tt_fake_a1a11pr'],       

        'tt_tau_pirho': ['tt_tau_pirho', 'tt_tau_rhopi'],
        'tt_tau_rhoa1': ['tt_tau_a1rho', 'tt_tau_rhoa1'],
        'tt_tau_pia1': ['tt_tau_a1pi', 'tt_tau_pia1'],
        'tt_tau_pia11pr': ['tt_tau_pia11pr', 'tt_tau_a11prpi'],
        'tt_tau_a11pra1': ['tt_tau_a11pra1', 'tt_tau_a1a11pr'],
    }

    # add aiso directories to the combinations
    dir_combinations_extra = {}
    for dir_name in dir_combinations.keys():
        dir_combinations_extra[dir_name + '_aiso'] = [d + '_aiso' for d in dir_combinations[dir_name]]
        dir_combinations_extra[dir_name + '_aco_aiso'] = [d + '_aco_aiso' for d in dir_combinations[dir_name]]
        dir_combinations_extra[dir_name + '_BDT_score_aiso'] = [d + '_BDT_score_aiso' for d in dir_combinations[dir_name]]
    dir_combinations.update(dir_combinations_extra)

    # Call the hadd function with the provided arguments
    hadd_root_files(input_files, output_file, dir_combinations, exp_num=exp_num)
    # should be called added_histo_Run2Bins.root
