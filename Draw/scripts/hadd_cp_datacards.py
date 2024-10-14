import ROOT
import argparse
import os

def hadd_root_files(input_files, output_file, dir_combinations):
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
    
    # Dictionary to store the combined histograms
    combined_dirs_histograms = {}
   
    # Dictionary to store histograms for directories not in dir_combinations
    independent_dirs_histograms = {}

    # Open each input ROOT file
    for file_name in input_files:
        input_root = ROOT.TFile.Open(file_name, "READ")
        
        # Loop over all the keys in the input file
        for key in input_root.GetListOfKeys():
            obj = key.ReadObj()
            
            # Check if the object is a directory
            if isinstance(obj, ROOT.TDirectoryFile):
                dir_name = obj.GetName()
                
                # Find if the directory matches any in dir_combinations
                combined = False
                for combined_dir, dirs_to_combine in dir_combinations.items():
                    if dir_name in dirs_to_combine:
                        combined = True
                        # If this combined directory doesn't exist yet, create it
                        if combined_dir not in combined_dirs_histograms:
                            combined_dirs_histograms[combined_dir] = {}
                        
                        # Loop over all histograms in the directory
                        input_root.cd(dir_name)
                        for hist_key in ROOT.gDirectory.GetListOfKeys():
                            hist_name = hist_key.GetName()
                            hist = hist_key.ReadObj()
                            
                            # Check if it's a histogram
                            if isinstance(hist, ROOT.TH1):
                                # If it's the first time we encounter this histogram, clone it
                                if hist_name not in combined_dirs_histograms[combined_dir]:
                                    combined_dirs_histograms[combined_dir][hist_name] = hist.Clone()
                                    combined_dirs_histograms[combined_dir][hist_name].SetDirectory(0)
                                else:
                                    # Otherwise, add the histogram to the existing one
                                    combined_dirs_histograms[combined_dir][hist_name].Add(hist)
   
                # If directory is not in dir_combinations, handle it separately
                if not combined:
                    if dir_name not in independent_dirs_histograms:
                        independent_dirs_histograms[dir_name] = {}   

                    # Loop over all histograms in the directory
                    input_root.cd(dir_name)
                    for hist_key in ROOT.gDirectory.GetListOfKeys():
                        hist_name = hist_key.GetName()
                        hist = hist_key.ReadObj()
                        
                        # Check if it's a histogram
                        if isinstance(hist, ROOT.TH1):
                            # If it's the first time we encounter this histogram, clone it
                            if hist_name not in independent_dirs_histograms[dir_name]:
                                independent_dirs_histograms[dir_name][hist_name] = hist.Clone()
                                independent_dirs_histograms[dir_name][hist_name].SetDirectory(0)
                            else:
                                # Otherwise, add the histogram to the existing one
                                independent_dirs_histograms[dir_name][hist_name].Add(hist)


    all_dirs_histograms = combined_dirs_histograms | independent_dirs_histograms

    for dir_name, histograms in all_dirs_histograms.items():
        # Create the directory in the output file
        output.mkdir(dir_name)
        output.cd(dir_name)

        # Write all histograms in this directory
        for hist_name, hist in histograms.items():
            hist.Write(hist_name)

    # Close the output file
    output.Close()

if __name__ == "__main__":
    # Argument parser for command-line arguments
    parser = argparse.ArgumentParser()
    
    # Add input and output files arguments
    parser.add_argument('-i', '--input', nargs='+', required=True, help="List of input ROOT files")
    parser.add_argument('-o', '--output', required=True, help="Name of the output ROOT file")
    
    # Parse arguments
    args = parser.parse_args()

    # List of input ROOT files from command-line
    input_files = args.input

    # Output ROOT file from command-line
    output_file = args.output

    # Define which directories to combine (customize this based on your case)
    dir_combinations = {
        'tt_higgs_pirho': ['tt_higgs_pirho', 'tt_higgs_rhopi'],
        'tt_higgs_rhoa1': ['tt_higgs_a1rho', 'tt_higgs_rhoa1'],
        'tt_higgs_pia1': ['tt_higgs_a1pi', 'tt_higgs_pia1'],
        'tt_higgs_pia11pr': ['tt_higgs_pia11pr', 'tt_higgs_a11prpi'],
        'tt_higgs_a11pra1': ['tt_higgs_a11pra1', 'tt_higgs_a1a11pr'],
    }

    # Call the hadd function with the provided arguments
    hadd_root_files(input_files, output_file, dir_combinations)

