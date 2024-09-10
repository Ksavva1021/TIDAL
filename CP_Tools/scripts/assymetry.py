import uproot3 as uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D
import pandas as pd
import sys


def AcoPlot(bins, r, x, weight_sm, weight_ps, title, savefig):
    # Create a new figure and axis
    fig, (ax, ax_ratio) = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.02})
    hep.cms.label(rlabel="", ax=ax,fontsize=13)
    # Histogram for Scalar weights
    hist1, bin_edges, _ = ax.hist(x, bins=bins, range=r, color='b', weights=weight_sm,
                                  label='CP-even', histtype='step', lw=1)
    # Histogram for Pseudoscalar weights
    hist2, bin_edges, _ = ax.hist(x, bins=bins, range=r, color='r', weights=weight_ps,
                                  label='CP-odd', histtype='step', lw=1)

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(hist2, hist1, out=np.zeros_like(hist2), where=hist1!=0)
    ax_ratio.plot((bin_edges[:-1] + bin_edges[1:]) / 2, ratio, color='red', marker='_', linestyle='none')
    ax_ratio.axhline(y=1, color='blue', linestyle='-', linewidth=1)

    # Calculate asymmetry between Scalar and Pseudoscalar
    asymmetry_overall = (1 / bins) * (np.sum(abs(hist1 - hist2) / (hist1 + hist2)))
    # Add text annotation to the plot
    ax.text(0.94, 1.02, "A={:.3f}".format(asymmetry_overall), ha='center', va='center', transform=ax.transAxes, fontsize=12)

    # Set titles and labels
    ax.set_title(title,fontsize=12)
    ax.set_ylabel('a.u.', fontsize=12)
    ax_ratio.set_xlabel('$\Delta\phi$', fontsize=13)
    ax_ratio.set_ylabel('Ratio', fontsize=12)

    ax.set_xlim([0, 2 * np.pi])  # Set the x-axis range from 0 to 2π
    ax_ratio.set_xlim([0, 2 * np.pi])  # Set the x-axis range from 0 to 2π

    # Hide x-tick labels on ax to avoid overlap
    ax.tick_params(axis='x', which='both', labelbottom=False)
    ax.tick_params(axis='both', which='both', labelsize=11)
    ax_ratio.tick_params(axis='both', which='both', labelsize=11)
    ax.minorticks_on()
    ax_ratio.minorticks_on()

    # Add legend to the main plot
    legend_entries = [Line2D([0], [0], color='b', lw=1, label='CP-even'),
                     Line2D([0], [0], color='r', lw=1, label='CP-odd')
    ]

    ax.legend(loc=(0.29,0.1),handles=legend_entries, ncol=2, frameon=False, fontsize=11)

    # Save the figure
    fig.savefig(savefig)
    plt.close(fig)


# Open the ROOT file
Run = "Run3"
DM = sys.argv[1]
deeptau_wp = sys.argv[2]
deeptau_wp = int(deeptau_wp)

if Run == "Run2":
    # for Run 3 read from a pickle file
    df = pd.read_pickle("/vols/cms/ks1021/offline/Htt_CPinDecay/Regression/samples/pickle/ggH_DM_1010.pkl")
else:
    file = uproot.open("/vols/cms/ks1021/TIDAL/CP_Tools/samples/merged_cp.root")
    tree = file["ntuple"]

    df = tree.pandas.df()
    df = df[(df['os'] == 1)]
    #df = df[(df['idDeepTau2018v2p5VSjet_1'] >= deeptau_wp)]
    #df = df[(df['idDeepTau2018v2p5VSjet_2'] >= deeptau_wp)]
    #df = df[(df['ip_LengthSig_1'] > 1.5) & (df['ip_LengthSig_2'] > 1.5)]
    #df = df[(df['decayMode_1'] == 10) & (df['decayMode_2'] == 10)]
    df.reset_index(drop=True, inplace=True)

if Run == "Run2":
    variable1 = "pv_angle"
    data1 = df[variable1].to_numpy()
    weight_sm = df["wt_cp_sm"].to_numpy()
    weight_ps = df["wt_cp_ps"].to_numpy()
else:
    variable1 = "aco_angle_" + DM
    variable2 = "wt_cp_sm"
    variable3 = "wt_cp_ps"

    # Extract the data
    data1 = df[variable1].to_numpy()
    weight_sm = df[variable2].to_numpy()
    weight_ps = df[variable3].to_numpy()

bins = 20
r = (0, 2 * np.pi)

#if DM == "0_0":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $\pi$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")
#elif DM == "0_1":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $\rho$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")
#elif DM == "0_10":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $a_{1}$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")
#elif DM == "1_1":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\rho$ - $\rho$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")
#elif DM == "1_10":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\rho$ - $a_{1}$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")
#elif DM == "10_10":
#    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$a_{1}$ - $a_{1}$",f"CP_Tools/plots_sel_IP/WP_{deeptau_wp}/Aco_CP_{Run}_{DM}.png")

if DM == "0_0":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $\pi$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")
elif DM == "0_1":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $\rho$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")
elif DM == "0_10":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\pi$ - $a_{1}$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")
elif DM == "1_1":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\rho$ - $\rho$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")
elif DM == "1_10":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$\rho$ - $a_{1}$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")
elif DM == "10_10":
    AcoPlot(bins,r,data1,weight_sm,weight_ps,r"$a_{1}$ - $a_{1}$",f"CP_Tools/plots/Aco_CP_{Run}_{DM}.png")