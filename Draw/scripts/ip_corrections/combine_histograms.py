from array import array
import ROOT
import glob
import re
import logging
from prettytable import PrettyTable
from collections import defaultdict
import os


# Setup logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def FindRebinning(hist, BinThreshold=100, BinUncertFraction=0.5):
    # Get initial binning
    binning = [hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)]

    def rebin_histogram(hist, binning, new_name_suffix):
        """Helper function to rebin histograms and assign unique names."""
        bin_array = array('d', binning)  # Convert to C-style array
        rebinned_hist = hist.Rebin(len(bin_array) - 1, f"{hist.GetName()}_{new_name_suffix}", bin_array)
        return rebinned_hist

    # Remove outer bins with zero content and error
    still_zero = True
    for i in range(1, hist.GetNbinsX()):
        if not (hist.GetBinContent(i) == 0 and hist.GetBinError(i) == 0):
            still_zero = False
        if still_zero:
            binning.remove(min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i))))
        else:
            break
    hist = rebin_histogram(hist, binning, "outer_left")

    still_zero = True
    for i in reversed(range(2, hist.GetNbinsX() + 1)):
        if not (hist.GetBinContent(i) == 0 and hist.GetBinError(i) == 0):
            still_zero = False
        if still_zero:
            binning.remove(min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i + 1))))
        else:
            break
    hist = rebin_histogram(hist, binning, "outer_right")

    # Rebin left to right based on uncertainty fraction and bin threshold
    finished = False
    k = 0
    while not finished and k < 1000:
        k += 1
        for i in range(1, hist.GetNbinsX()):
            if hist.GetBinContent(i) > 0:
                uncert_frac = hist.GetBinError(i) / hist.GetBinContent(i)
            else:
                uncert_frac = BinUncertFraction + 1
            if uncert_frac > BinUncertFraction and hist.GetBinContent(i) < BinThreshold:
                binning.remove(min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i + 1))))
                hist = rebin_histogram(hist, binning, f"left_{k}")
                break
            elif i + 1 == hist.GetNbinsX():
                finished = True

    # Rebin right to left based on uncertainty fraction and bin threshold
    finished = False
    k = 0
    while not finished and k < 1000:
        k += 1
        for i in reversed(range(2, hist.GetNbinsX() + 1)):
            if hist.GetBinContent(i) > 0:
                uncert_frac = hist.GetBinError(i) / hist.GetBinContent(i)
            else:
                uncert_frac = BinUncertFraction + 1
            if uncert_frac > BinUncertFraction and hist.GetBinContent(i) < BinThreshold:
                binning.remove(min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i))))
                hist = rebin_histogram(hist, binning, f"right_{k}")
                break
            elif i == 2:
                finished = True

    return binning


def RebinHist(hist, binning, new_name_suffix="rebinned"):
    # Convert binning to a C-style array
    new_binning = array("f", map(float, binning))

    # Create a new histogram with a unique name
    rebinned_name = f"{hist.GetName()}_{new_name_suffix}"
    hout = ROOT.TH1D(rebinned_name, hist.GetTitle(), len(new_binning) - 1, new_binning)

    # Fill the rebinned histogram with content and errors
    total_entries = 0  # Track total entries
    for i in range(1, hout.GetNbinsX() + 1):
        for j in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinCenter(j) >= hout.GetBinLowEdge(i) and hist.GetBinCenter(j) < hout.GetBinLowEdge(i + 1):
                new_content = hout.GetBinContent(i) + hist.GetBinContent(j)
                new_error = (hout.GetBinError(i) ** 2 + hist.GetBinError(j) ** 2) ** 0.5
                hout.SetBinContent(i, new_content)
                hout.SetBinError(i, new_error)
                total_entries += hist.GetBinContent(j)  # Accumulate entries

    # Add underflow and overflow contributions
    hout.SetBinContent(1, hout.GetBinContent(1) + hist.GetBinContent(0))
    hout.SetBinError(1, (hout.GetBinError(1) ** 2 + hist.GetBinError(0) ** 2) ** 0.5)
    total_entries += hist.GetBinContent(0)  # Add underflow entries

    hout.SetBinContent(hout.GetNbinsX(), hout.GetBinContent(hout.GetNbinsX()) + hist.GetBinContent(hist.GetNbinsX() + 1))
    hout.SetBinError(hout.GetNbinsX(), (hout.GetBinError(hout.GetNbinsX()) ** 2 + hist.GetBinError(hist.GetNbinsX() + 1) ** 2) ** 0.5)
    total_entries += hist.GetBinContent(hist.GetNbinsX() + 1)  # Add overflow entries

    # Update total entries in the rebinned histogram
    hout.SetEntries(hist.GetEntries())  # Preserve the original number of entries

    return hout


# Function to log histogram yields
def log_yield(histogram):
    return histogram.Integral() if histogram else None


# Function to retrieve histograms from a single file
def get_histograms(file_name, process_names, channel):
    histograms = {}
    file = ROOT.TFile.Open(file_name)
    if not file or file.IsZombie():
        logger.warning(f"Could not open file {file_name}")
        return histograms

    for process_name in process_names:
        hist_path = f"{channel}/{process_name}"
        hist = file.Get(hist_path)
        if hist:
            hist.SetDirectory(0)  # Detach histogram from file so we can close the file
            histograms[process_name] = hist
        else:
            logger.warning(f"Histogram '{hist_path}' not found in {file_name}")

    file.Close()
    return histograms


# Function to subtract histograms from data_obs
def subtract_histograms(data_histogram, process_histograms):
    result = data_histogram.Clone()
    for hist in process_histograms:
        result.Add(hist, -1)  # Subtract each histogram
    return result


def get_eta_index(variable):
    # Look for a digit at the end of the variable name
    match = re.search(r'_(\d+)', variable)
    return match.group(1) if match else None


# Function to process files based on variable, eta binning, and channel
def process_files(directory, variables, eta_boundaries, channel, year, merge_mapping):

    input_directory = f"{directory}/{channel}"
    output_directory = f"{directory}/histograms"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define the process names needed
    all_processes = ["ZL", "data_obs", "ZTT", "ZJ", "TTT", "TTJ", "VVT", "VVJ", "W", "QCD"]

    # Determine channel path
    channel_path = f"{channel}_inclusive"

    # Generate eta ranges by pairing adjacent values in eta_boundaries
    eta_ranges = [(eta_boundaries[i], eta_boundaries[i + 1]) for i in range(len(eta_boundaries) - 1)]

    # Dictionary to store intermediate histograms for each variable and eta binning
    histograms_by_eta = defaultdict(lambda: {"ZL": None, "data_obs": None})

    # Table to store initial and final yields
    initial_table = PrettyTable()
    final_table = PrettyTable()

    # Set up table headers
    initial_table.field_names = ["Variable", "Eta Binning"] + all_processes
    final_table.field_names = ["Combined Variable", "Eta Binning", "ZL", "data_obs (after subtraction)"]

    for variable in variables:
        # Determine eta index from the variable name
        eta_index = get_eta_index(variable)
        if not eta_index:
            logger.warning(f"Could not determine eta index for variable '{variable}', skipping.")
            continue

        for eta_min, eta_max in eta_ranges:
            # Format eta binning for the file pattern
            eta_binning = f"{str(eta_min).replace('.', 'p')}to{str(eta_max).replace('.', 'p')}"

            # Pattern for the specific file
            pattern = f"{input_directory}/datacard_{variable}_eta_{eta_index}_{eta_binning}_inclusive_{channel}_{year}.root"
            files = glob.glob(pattern)

            if not files:
                logger.warning(f"No file found for pattern: {pattern}")
                continue

            file_name = files[0]  # Since there's only one file per pattern

            # Retrieve histograms from the file
            histograms = get_histograms(file_name, all_processes, channel_path)

            # Initialize row for initial table
            row_yields = {process: log_yield(histograms.get(process)) for process in all_processes}
            initial_table.add_row([variable, eta_binning] + [row_yields[process] for process in all_processes])

            # Get the ZL and data_obs histograms, subtracting other processes from data_obs
            zl_hist = histograms.get("ZL")
            data_obs_hist = histograms.get("data_obs")
            if zl_hist and data_obs_hist:
                # Subtract other processes from data_obs
                other_process_histograms = [histograms[proc] for proc in all_processes if proc not in ["ZL", "data_obs"] and proc in histograms]
                data_obs_final = subtract_histograms(data_obs_hist, other_process_histograms)

                # Store ZL and data_obs final histograms for the variable and eta binning
                key = (variable, eta_binning)
                histograms_by_eta[key]["ZL"] = zl_hist.Clone()
                histograms_by_eta[key]["data_obs"] = data_obs_final.Clone()
            else:
                logger.warning(f"Missing ZL or data_obs histogram for variable '{variable}' and eta binning '{eta_binning}'.")

    # Combine histograms based on the user-defined merge mapping
    combined_histograms = {}

    for combined_var, var_list in merge_mapping.items():
        for eta_min, eta_max in eta_ranges:
            eta_binning = f"{str(eta_min).replace('.', 'p')}to{str(eta_max).replace('.', 'p')}"

            zl_combined = None
            data_obs_combined = None

            for var in var_list:
                key = (var, eta_binning)
                if key in histograms_by_eta:
                    zl_hist = histograms_by_eta[key]["ZL"]
                    data_obs_hist = histograms_by_eta[key]["data_obs"]

                    if zl_hist:
                        if zl_combined is None:
                            zl_combined = zl_hist.Clone()
                        else:
                            zl_combined.Add(zl_hist)

                    if data_obs_hist:
                        if data_obs_combined is None:
                            data_obs_combined = data_obs_hist.Clone()
                        else:
                            data_obs_combined.Add(data_obs_hist)

            if zl_combined and data_obs_combined:
                combined_key = f"{combined_var}_{eta_binning}"
                combined_histograms[combined_key] = {"ZL": zl_combined, "data_obs": data_obs_combined}

                # Append to the final yields table
                final_table.add_row([combined_var, eta_binning, zl_combined.Integral(), data_obs_combined.Integral()])

    # Open a single ROOT file for saving all final histograms, including the year in the name
    output_file_name = f"{output_directory}/combined_histograms_{channel}_{year}.root"
    output_file = ROOT.TFile(output_file_name, "RECREATE")

    # Save each final histogram to the output file
    for key, hist_data in combined_histograms.items():
        hist_data["ZL"].Write(f"{key}_MC")
        hist_data["data_obs"].Write(f"{key}_data")

    output_file.Close()
    logger.info(f"All combined histograms saved in '{output_file_name}'.")

    # Print all tables
    print("\nInitial Yields Table")
    print(initial_table)
    print("\nFinal Yields Table")
    print(final_table)

    # Perform rebinning of the histograms
    outfile_rebin = ROOT.TFile(f"{output_directory}/combined_histograms_{channel}_{year}_rebin.root", "RECREATE")

    # binning array first eta bin
    binning_dict = {}

    for key in merge_mapping.keys():
        first_eta_bin = str(eta_boundaries[0]).replace(".","p") + "to" + str(eta_boundaries[1]).replace(".", "p")
        binning_dict[key] = FindRebinning(combined_histograms[f"{key}_{first_eta_bin}"]["ZL"].Clone())
        #print(key,len(binning_dict[key]))

    # look at MC to decide the new binning and apply it to data_obs
    for key, hist_data in combined_histograms.items():
        var = "_".join(key.split("_")[:-1])
        binning = binning_dict[var]
        hist_data_ZL = RebinHist(hist_data["ZL"].Clone(), binning)
        hist_data_data_obs = RebinHist(hist_data["data_obs"].Clone(), binning)
        #print(var, len(binning), hist_data_ZL.GetNbinsX(), hist_data_data_obs.GetNbinsX())
        hist_data_ZL.Write(f"{key}_MC")
        hist_data_data_obs.Write(f"{key}_data")
        del hist_data_ZL
        del hist_data_data_obs

    outfile_rebin.Close()

# Define variables, eta boundaries, channel, and year
variables = [
    "ip_x_1", "ip_x_2",
    "ip_y_1", "ip_y_2",
    "ip_z_1", "ip_z_2",
    "log_err00_1", "log_err00_2",
    "log_err11_1", "log_err11_2",
    "log_err22_1", "log_err22_2",
    "ip_LengthSig_1", "ip_LengthSig_2",
    "ip_LengthTotal_1", "ip_LengthTotal_2",
    "ip_ErrorTotal_1", "ip_ErrorTotal_2",
#    "ip_alt_LengthSig_1", "ip_alt_LengthSig_2",
#    "ip_x_1_Err", "ip_x_2_Err",
#    "ip_y_1_Err", "ip_y_2_Err",
#    "ip_z_1_Err", "ip_z_2_Err",
#    "ip_cov10_1", "ip_cov10_2",
#    "ip_cov20_1", "ip_cov20_2",
#    "ip_cov21_1", "ip_cov21_2",
#    "ip_x_1_Err_ratio", "ip_x_2_Err_ratio",
#    "ip_y_1_Err_ratio", "ip_y_2_Err_ratio",
#    "ip_z_1_Err_ratio", "ip_z_2_Err_ratio",
#    "ip_x_1_log_Err", "ip_x_2_log_Err",
#    "ip_y_1_log_Err", "ip_y_2_log_Err",
#    "ip_z_1_log_Err", "ip_z_2_log_Err",
#    "ip_cov00_1", "ip_cov00_2",
#    "ip_cov11_1", "ip_cov11_2",
#    "ip_cov22_1", "ip_cov22_2",
#    "ip_cov10_1", "ip_cov10_2",
#    "ip_cov20_1", "ip_cov20_2",
#    "ip_cov21_1", "ip_cov21_2",
]

channel = "mm"

if channel == "mm":
    eta_boundaries = [0.0, 0.9, 1.2, 2.1, 2.4]
elif channel == "ee":
    eta_boundaries = [0.0, 1.0, 1.48, 1.65, 2.1]
else:
    logger.error(f"Channel '{channel}' not supported.")
    exit(1)

year = "Run3_2022"  # Specify the year
directory = f"/vols/cms/ks1021/TIDAL/Draw/plots/IP_corrections/before_2/{year}/ip_calculation/"

# Define the merge mapping
merge_mapping = {
    "ip_x": ["ip_x_1", "ip_x_2"],
    "ip_y": ["ip_y_1", "ip_y_2"],
    "ip_z": ["ip_z_1", "ip_z_2"],
    "log_err00": ["log_err00_1", "log_err00_2"],
    "log_err11": ["log_err11_1", "log_err11_2"],
    "log_err22": ["log_err22_1", "log_err22_2"],
    "ip_LengthSig": ["ip_LengthSig_1", "ip_LengthSig_2"],
    "ip_LengthTotal": ["ip_LengthTotal_1", "ip_LengthTotal_2"],
    "ip_ErrorTotal": ["ip_ErrorTotal_1", "ip_ErrorTotal_2"],
    # "PCA_ip_cov00": ["PCA_ip_cov00_1", "PCA_ip_cov00_2"],
    # "PCA_ip_cov11": ["PCA_ip_cov11_1", "PCA_ip_cov11_2"],
    # "PCA_ip_cov22": ["PCA_ip_cov22_1", "PCA_ip_cov22_2"],
    # "PCA_ip_cov10": ["PCA_ip_cov10_1", "PCA_ip_cov10_2"],
    # "PCA_ip_cov20": ["PCA_ip_cov20_1", "PCA_ip_cov20_2"],
    # "PCA_ip_cov21": ["PCA_ip_cov21_1", "PCA_ip_cov21_2"],
    #"ip_alt_LengthSig": ["ip_alt_LengthSig_1", "ip_alt_LengthSig_2"],
#    "ip_x_Err": ["ip_x_1_Err", "ip_x_2_Err"],
#    "ip_y_Err": ["ip_y_1_Err", "ip_y_2_Err"],
#    "ip_z_Err": ["ip_z_1_Err", "ip_z_2_Err"],
    #"ip_cov10": ["ip_cov10_1", "ip_cov10_2"],
    #"ip_cov20": ["ip_cov20_1", "ip_cov20_2"],
    #"ip_cov21": ["ip_cov21_1", "ip_cov21_2"],
#    "ip_x_Err_ratio": ["ip_x_1_Err_ratio", "ip_x_2_Err_ratio"],
#    "ip_y_Err_ratio": ["ip_y_1_Err_ratio", "ip_y_2_Err_ratio"],
#    "ip_z_Err_ratio": ["ip_z_1_Err_ratio", "ip_z_2_Err_ratio"],
    #"ip_x_log_Err": ["ip_x_1_log_Err", "ip_x_2_log_Err"],
    #"ip_y_log_Err": ["ip_y_1_log_Err", "ip_y_2_log_Err"],
    #"ip_z_log_Err": ["ip_z_1_log_Err", "ip_z_2_log_Err"],
#    "ip_cov00": ["ip_cov00_1", "ip_cov00_2"],
#    "ip_cov11": ["ip_cov11_1", "ip_cov11_2"],
#    "ip_cov22": ["ip_cov22_1", "ip_cov22_2"],
#    "ip_cov10": ["ip_cov10_1", "ip_cov10_2"],
#    "ip_cov20": ["ip_cov20_1", "ip_cov20_2"],
#    "ip_cov21": ["ip_cov21_1", "ip_cov21_2"],
}

# Run the processing
process_files(directory, variables, eta_boundaries, channel, year, merge_mapping)
