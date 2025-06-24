import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import re

hep.style.use("CMS")

channel = "mm"
year = "Run3_2023BPix"
input_dir = f"/vols/cms/ks1021/TIDAL/Draw/plots/IP_corrections/{year}/ip_calculation/histograms"

output_dir_png = os.path.join(input_dir, channel, year, "pngs")
output_dir_pdf = os.path.join(input_dir, channel, year, "pdfs")
os.makedirs(output_dir_png, exist_ok=True)
os.makedirs(output_dir_pdf, exist_ok=True)

# Load ROOT file
file_path = f"{input_dir}/combined_histograms_{channel}_{year}_rebin.root"
root_file = uproot.open(file_path)

# Variable names
variable_mapping = {
    "ip_x": r"$\mathrm{IP}_x$ [cm]",
    "ip_y": r"$\mathrm{IP}_y$ [cm]",
    "ip_z": r"$\mathrm{IP}_z$ [cm]",
}

sorted_keys = sorted(variable_mapping.keys(), key=len, reverse=True)
pattern = re.compile(r"(" + "|".join(map(re.escape, sorted_keys)) + r")_.*")

# Utilities
def clean_name(name):
    return re.sub(r";\d+", "", name)

def get_eta_range(name):
    eta_range = re.search(r"(\d+p\d+to\d+p\d+)", name).group(1)
    eta = eta_range.replace("p", ".")
    eta_low, eta_high = eta.split("to")
    return eta_low, eta_high

# Extract histograms
hist_names = [name for name, obj in root_file.items() if obj.classname.startswith("TH1")]
cleaned_names = {clean_name(name): name for name in hist_names}

# Loop through data/MC pairs
for clean_hist_name in cleaned_names:
    if "_data" in clean_hist_name:
        base_name = clean_hist_name.replace("_data", "")
        mc_name = base_name + "_MC"

        if mc_name in cleaned_names:
            data_hist_name = cleaned_names[clean_hist_name]
            mc_hist_name = cleaned_names[mc_name]

            # Load histograms
            data_hist = root_file[data_hist_name]
            mc_hist = root_file[mc_hist_name]

            data_vals, data_edges = data_hist.to_numpy()
            data_vals = np.where(data_vals < 0, 0, data_vals)
            mc_vals, mc_edges = mc_hist.to_numpy()
            mc_vals = np.where(mc_vals < 0, 0, mc_vals)

            data_centers = (data_edges[:-1] + data_edges[1:]) / 2
            bin_widths = data_edges[1:] - data_edges[:-1]

            # MC uncertainty — use sqrt(N) as placeholder
            mc_errors = np.sqrt(mc_vals)

            # Ratio and errors
            ratio_vals = np.divide(data_vals, mc_vals, out=np.zeros_like(data_vals), where=mc_vals != 0)
            ratio_errors_data = np.divide(np.sqrt(data_vals), mc_vals, out=np.zeros_like(data_vals), where=mc_vals != 0)
            ratio_errors_data = np.where(np.isnan(ratio_errors_data), 0, ratio_errors_data)

            ratio_errors_MC_down = np.divide(mc_errors, mc_vals, out=np.zeros_like(mc_errors), where=mc_vals != 0)
            ratio_errors_MC_up = ratio_errors_MC_down.copy()  # symmetric

            # Setup figure
            fig, (ax, ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True, figsize=(9, 8))

            # MAIN PLOT

            # Light fill under MC
            ax.fill_between(mc_edges[:-1], mc_vals, step="post", color='royalblue', alpha=0.3, linewidth=0)

            # MC step line
            ax.step(mc_edges[:-1], mc_vals, where='post', color='royalblue', linewidth=1.5, label='ZLL')

            # MC uncertainty band with hatch
            ax.fill_between(mc_edges,
                            np.append(mc_vals - mc_errors, 0),
                            np.append(mc_vals + mc_errors, 0),
                            step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0,
                            label="Background Uncertainty")

            # Data points with x-errors
            ax.errorbar(data_centers, data_vals, yerr=np.sqrt(data_vals), xerr=bin_widths/2,
                        fmt='o', color='black', markersize=3, linewidth=0.6, label='Observation')

            # Log y-scale if needed
            if "log" in base_name or "cov" in base_name:
                ax.set_yscale("log")
                ax.set_ylim(0.1, 10 * np.max(mc_vals))
            else:
                ax.set_ylim(0, 2.0 * np.max(mc_vals))

            # CMS label
            if year == "Run3_2022":
                lumi = 7.98
            elif year == "Run3_2022EE":
                lumi = 26.67
            elif year == "Run3_2023":
                lumi = 17.79
            elif year == "Run3_2023BPix":
                lumi = 9.45
            hep.cms.label(ax=ax, label="Preliminary", data=True, lumi=lumi, com=13.6, fontsize=16)

            # eta range label
            eta_low, eta_high = get_eta_range(base_name)
            eta_text = rf"${eta_low} < |\eta| < {eta_high}$"
            ax.text(0.035, 0.86, eta_text, fontsize=16, fontweight="bold", transform=ax.transAxes)

            # Legend
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], loc='upper right', frameon=1, framealpha=1, bbox_to_anchor=(0.98, 0.98))

            # Y label
            ax.set_ylabel("Events", fontsize=14)
            ax.tick_params(axis='both', labelsize=13)
            ax.tick_params(labelbottom=False)

            # RATIO PLOT

            # Ref line at 1
            ax_ratio.axhline(1, color='black', linestyle=':', linewidth=1)

            # Ratio points
            ax_ratio.errorbar(data_centers, ratio_vals, yerr=ratio_errors_data, xerr=bin_widths/2,
                              fmt='o', color='black', markersize=3, linewidth=0.6)

            # MC uncertainty band in ratio plot (NO fill below MC)
            ax_ratio.fill_between(mc_edges,
                                  1 - np.append(ratio_errors_MC_down, 0),
                                  1 + np.append(ratio_errors_MC_up, 0),
                                  step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0)

            # X and Y labels
            match = pattern.match(base_name)
            if match:
                variable = variable_mapping[match.group(1)]
            else:
                variable = base_name
            ax_ratio.set_xlabel(variable, fontsize=14)
            ax_ratio.set_ylabel("Obs/Exp", fontsize=12)
            ax_ratio.set_ylim(0.5, 1.5)
            ax_ratio.tick_params(axis='both', labelsize=12)

            fig.subplots_adjust(hspace=0.05)

            # Save PNG and PDF
            output_name = f"{base_name}_comparison_with_ratio"
            save_path_png = os.path.join(output_dir_png, f"{output_name}.png")
            save_path_pdf = os.path.join(output_dir_pdf, f"{output_name}.pdf")

            fig.savefig(save_path_png, bbox_inches='tight', dpi=300)
            print(f"Saved PNG to {save_path_png}")

            fig.savefig(save_path_pdf, bbox_inches='tight')
            print(f"Saved PDF to {save_path_pdf}")

            plt.close(fig)
