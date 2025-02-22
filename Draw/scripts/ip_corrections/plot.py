import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import re
from matplotlib.gridspec import GridSpec

hep.style.use("CMS")

channel = "mm"
year = "Run3_2022"
input_dir = f"/vols/cms/ks1021/TIDAL/Draw/plots/IP_corrections/before_2/{year}/ip_calculation/histograms"

output_dir = os.path.join(input_dir, channel, year)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load the ROOT file
file_path = f"{input_dir}/combined_histograms_{channel}_{year}_rebin.root"  # Replace with your actual file path
root_file = uproot.open(file_path)

variable_mapping = {
    "ip_x": "Impact parameter x [cm]",
    "ip_y": "Impact parameter y [cm]",
    "ip_z": "Impact parameter z [cm]",
    "log_err00": "Log error 00",
    "log_err11": "Log error 11",
    "log_err22": "Log error 22",
    "ip_LengthSig": "Impact parameter significance",
    "ip_LengthTotal": "Impact parameter total length [cm]",
    "ip_ErrorTotal": "Impact parameter total error [cm]",
    # "PCA_ip_cov00": "PCA covariance 00",
    # "PCA_ip_cov11": "PCA covariance 11",
    # "PCA_ip_cov22": "PCA covariance 22",
    # "PCA_ip_cov10": "PCA covariance 10",
    # "PCA_ip_cov20": "PCA covariance 20",
    # "PCA_ip_cov21": "PCA covariance 21",
#    "ip_alt_LengthSig": "Impact parameter significance (alternative)",
#    "ip_x_Err": "Impact parameter x error [cm]",
#    "ip_y_Err": "Impact parameter y error [cm]",
#    "ip_z_Err": "Impact parameter z error [cm]",
#    "ip_x_Err_ratio": "Impact parameter x error ratio",
#    "ip_y_Err_ratio": "Impact parameter y error ratio",
#    "ip_z_Err_ratio": "Impact parameter z error ratio",
#    "ip_x_log_Err": "Impact parameter x error log scale",
#    "ip_y_log_Err": "Impact parameter y error log scale",
#    "ip_z_log_Err": "Impact parameter z error log scale",
#    "ip_cov00": "Covariance 00",
#    "ip_cov11": "Covariance 11",
#    "ip_cov22": "Covariance 22",
#    "ip_cov10": "Covariance 10",
#    "ip_cov20": "Covariance 20",
#    "ip_cov21": "Covariance 21",
}

sorted_keys = sorted(variable_mapping.keys(), key=len, reverse=True)
pattern = re.compile(r"(" + "|".join(map(re.escape, sorted_keys)) + r")_.*")

# Function to clean ROOT names by removing ';' and any following numbers
def clean_name(name):
    return re.sub(r";\d+", "", name)

def get_eta_range(name):
    eta_range = re.search(r"(\d+p\d+to\d+p\d+)", name).group(1)
    eta = eta_range.replace("p", ".")
    eta_low, eta_high = eta.split("to")
    return eta_low, eta_high

# Extract all histogram names, checking if each item is a histogram
hist_names = [name for name, obj in root_file.items() if obj.classname.startswith("TH1")]
cleaned_names = {clean_name(name): name for name in hist_names}  # Map cleaned names to original

# Loop over all _data and _MC pairs and plot each one
for clean_hist_name in cleaned_names:
    if "_data" in clean_hist_name:
        base_name = clean_hist_name.replace("_data", "")
        mc_name = base_name + "_MC"

        # Check if both _data and _MC versions exist
        if mc_name in cleaned_names:
            data_hist_name = cleaned_names[clean_hist_name]
            mc_hist_name = cleaned_names[mc_name]

            # Load the histograms
            data_hist = root_file[data_hist_name]
            mc_hist = root_file[mc_hist_name]

            # Convert histogram data to arrays for plotting
            data_vals, data_edges = data_hist.to_numpy()
            mc_vals, mc_edges = mc_hist.to_numpy()

            # Calculate bin centers for plotting data points
            data_centers = (data_edges[:-1] + data_edges[1:]) / 2

            # Calculate ratio and handle division by zero
            ratio_vals = np.divide(data_vals, mc_vals, out=np.zeros_like(data_vals), where=mc_vals != 0)
            ratio_errors = np.divide(np.sqrt(data_vals), mc_vals, out=np.zeros_like(data_vals), where=mc_vals != 0)

            # Set up figure with two subplots (main plot and ratio plot)
            fig = plt.figure(figsize=(12, 8))
            gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)

            # Main histogram plot
            ax_main = fig.add_subplot(gs[0])
            ax_main.errorbar(data_centers, data_vals, fmt='o', color="black", label="Data", markersize=3, capsize=2)
            hep.histplot(mc_vals, mc_edges, color="royalblue", alpha=0.8, histtype="fill", linewidth=1.5, label="ZLL", ax=ax_main)

            # Check if "Err" is in the base_name and set y-axis to log scale if true
            if "log" in base_name or "cov" in base_name:
                ax_main.set_yscale("log")

            if year == "Run3_2022":
                lumi = 7.98
            elif year == "Run3_2022EE":
                lumi = 26.67
            elif year == "Run3_2023":
                lumi = 17.79
            elif year == "Run3_2023BPIX":
                lumi = 9.45
            hep.cms.label(ax=ax_main, label="", data=False, lumi=lumi, com=13.6, fontsize=14)

            # Adding labels and customizations for main plot
            ax_main.legend(fontsize=13, loc="upper right")
            ax_main.set_ylabel("Events", fontsize=13)
            ax_main.tick_params(axis='both', labelsize=13)
            ax_main.tick_params(labelbottom=False)  # Hide x-axis labels for the main plot

            # Add text to indicate the eta range
            eta_low, eta_high = get_eta_range(base_name)
            eta_text = rf"${eta_low} < |\eta| < {eta_high}$"
            ax_main.text(0.80, 0.75, eta_text, transform=ax_main.transAxes, fontsize=14, bbox=dict(facecolor='white', alpha=0.8))

            # Ratio plot
            ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
            ax_ratio.errorbar(data_centers, ratio_vals, yerr=ratio_errors, fmt='o', color="black", markersize=3, capsize=2)
            ax_ratio.axhline(1, color="gray", linestyle="--", linewidth=1)  # Reference line at ratio = 1

            # Ratio plot labels and customizations
            match = pattern.match(base_name)
            if match:
                variable = variable_mapping[match.group(1)]
            else:
                variable = base_name
            ax_ratio.set_xlabel(variable, fontsize=13)
            ax_ratio.set_ylabel("Data / MC", fontsize=12)
            ax_ratio.set_ylim(0.5, 1.5)  # Set y-limits to focus on ratio deviations
            ax_ratio.tick_params(axis='both', labelsize=12)

            # Save plot with a unique name based on the base histogram name
            output_path = os.path.join(output_dir, f"{base_name}_comparison_with_ratio.png")
            plt.savefig(output_path, dpi=300)
            plt.close(fig)  # Close the figure to free memory
