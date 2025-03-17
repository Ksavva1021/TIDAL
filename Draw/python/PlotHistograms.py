import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
import uproot
import os

plt.style.use(hep.style.CMS)
plt.rcParams.update({"font.size": 16})


class HTT_Histogram:

    def __init__(self, file, category, channel, era, variable, variable_label=None, cp_channel=None, blind=False):
        self.file = uproot.open(file)
        self.file_name = file
        self.directory = os.path.dirname(file)
        self.category = category
        self.channel = channel
        self.cp_channel = cp_channel # specific decay
        self.era = era
        self.blind = blind
        self.variable = variable
        self.variable_label = variable_label

        self.initialize_plotting()
        self.initialize_nodes()

        # get yields
        self.get_backgrounds()
        self.get_data()

    def initialize_plotting(self):
        # get bin information
        self.bin_edges = self.file[f"{self.category}/ZTT"].axis().edges()
        self.bin_centers = self.file[f"{self.category}/ZTT"].axis().centers()
        self.bin_widths = np.diff(self.bin_edges)
        self.step_edges = np.append(self.bin_edges,2*self.bin_edges[-1]-self.bin_edges[-2]) # for outline
        # track height of stacked backgrounds
        self.stacked_block = np.zeros(len(self.bin_centers))
        self.stacked_step = np.zeros(len(self.bin_centers)+2)
        # define colour scheme
        self.colors = {
                "yellow": "#ffa90e",
                "red": "#e42536",
                "pink": "#f6b5fd", #facaff",
                "lightblue": "#3f90da",
                "violet": "#882ebd",
                "purple": "#964a8b",
                "brown": "#a96b59",
                "green": "#b9ac70",
                "grey": "#94a4a2",
                "ash": "#717581",
            }
        # define channel label (inclusive)
        label_map = {
                "tt": "$\\tau_h\\tau_h$",
                "mt": "$\\mu\\tau_h$",
                "et": "$e\\tau_h$",
                "em": "$e\\mu$",
                "ee": "$ee$",
                "mm": "$\\mu\\mu$",
            }
        self.channel_label = label_map[self.channel]
        # define CP channel label
        cp_label_map = {
                "rhorho": "$\\rho^0\\rho^0$",
            }
        if self.cp_channel is not None:
            self.channel_label = cp_label_map[self.cp_channel]


        if self.variable_label is None:
            # set variable label to variable name if not specified
            self.variable_label = self.variable


    def initialize_nodes(self):
        # get total background nodes
        if self.channel == "tt":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W", "ZL", "ZJ"], "color": "red"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        elif self.channel == "mt":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\mu\\mu$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        elif self.channel == "et":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to ee": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        elif self.channel == "em":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\ell\\ell": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        elif self.channel == "ee":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to ee": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        elif self.channel == "mm":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\mu\\mu$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
        # get luminosity
        if self.era == "Run3_2022":
            self.lumi = 7.98
        elif self.era == "Run3_2022EE":
            self.lumi = 26.67
        elif self.era == "Run3_2023":
            self.lumi = 17.79
        elif self.era == "Run3_2023BPix":
            self.lumi = 9.45
        # get color for each background
        for bkg, info in self.backgrounds.items():
            info['color'] = self.colors[info["color"]]  # replace color name with hex code

    def get_counts_errors(self, node):
        # get count and error for a given node
        path = f"{self.category}/{node}"
        histo = self.file[path]
        # We want these to be differential (divide by bin width)
        counts = histo.values()/self.bin_widths
        errors = histo.errors()/self.bin_widths
        return counts, errors

    def get_backgrounds(self):
        # track total background
        total_bkg_counts = np.zeros(len(self.bin_centers))
        total_bkg_sq_errors = np.zeros(len(self.bin_centers))
        # compute each background contribution
        for bkg, info in self.backgrounds.items():
            bkg_counts = np.zeros(len(self.bin_centers))
            bkg_sq_errors = np.zeros(len(self.bin_centers))
            for contribution in info["nodes"]:
                counts, errors = self.get_counts_errors(contribution)
                bkg_counts += counts
                bkg_sq_errors += errors**2
            # store counts and yields for each contribution
            self.backgrounds[bkg]["counts"] = bkg_counts
            self.backgrounds[bkg]["errors"] = np.sqrt(bkg_sq_errors)
            total_bkg_counts += bkg_counts
            total_bkg_sq_errors += bkg_sq_errors
        # store total background
        total_bkg_err = np.sqrt(total_bkg_sq_errors)
        total_bkg_pce_err = np.divide(total_bkg_err, total_bkg_counts, where=total_bkg_counts > 0, out=np.zeros_like(total_bkg_counts))
        self.total_background = {"counts": total_bkg_counts, "errors": total_bkg_err, "pc_error": total_bkg_pce_err}
        return True

    def get_data(self):
        # get counts for the data
        data_counts, data_errors = self.get_counts_errors("data_obs")
        pce_data = np.divide(data_errors, data_counts, where=data_counts > 0, out=np.zeros_like(data_counts))
        self.data = {"counts": data_counts, "errors": data_errors, "pc_error": pce_data}
        # get data/MC ratio
        count_ratio =  np.divide(data_counts, self.total_background["counts"], where=self.total_background["counts"] > 0, out=np.zeros_like(data_counts))
        self.data_MC_ratio = {"counts": count_ratio, "error_data": self.data["pc_error"], "error_MC": self.total_background["pc_error"]}
        return True

    def stack_background(self, ax, name):
        # add background to stacked histo
        counts = self.backgrounds[name]["counts"]
        steps = np.append(np.insert(counts,0,0.0),0.0)
        # plot block and outline
        self.ax.bar(self.bin_centers, counts, width = self.bin_widths, bottom = self.stacked_block,
                color = self.backgrounds[name]["color"], label = rf"{name}")
        self.ax.step(self.step_edges, steps + self.stacked_step, color='black', linewidth = 0.5)
        # Update the bottom for the next process
        self.stacked_block += counts
        self.stacked_step += steps
        return True

    def plot_1D_histo(self, ratio_min=0.5, ratio_max=1.5):
        print("Plotting 1D histogram")
        # plot 1D histogram
        self.fig, (self.ax, self.ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True, figsize=(9, 8))

        ## MAIN PLOT
        # add background contributions
        for bkg in self.backgrounds.keys():
            self.stack_background(self.ax, bkg)
        # add uncertainty on total MC
        self.stacked_block = np.insert(self.stacked_block, len(self.stacked_block), 0)
        self.ax.fill_between(self.bin_edges,
                    self.stacked_block-np.insert(self.total_background["errors"], len(self.total_background["errors"]), 0),
                    self.stacked_block+np.insert(self.total_background["errors"], len(self.total_background["errors"]), 0),
                    step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0, label = "Background Uncertainty")
        if not self.blind:
            # add data
            self.ax.errorbar(self.bin_centers, self.data['counts'], label='Observation', yerr=self.data['errors'], fmt='o', color = 'black', markersize=3, linewidth=0.6)
            self.ax.errorbar(self.bin_centers, self.data['counts'],xerr=self.bin_widths/2, fmt='o', color = 'black', markersize=3, linewidth=0.6) # add width marker

        ## RATIO PLOT
        # add data/MC ratio as points
        if not self.blind:
            self.ax_ratio.errorbar(self.bin_centers, self.data_MC_ratio['counts'], yerr=self.data_MC_ratio['error_data'], xerr=self.bin_widths/2, fmt='o', color = 'black', markersize=3, linewidth=0.6)
        # add MC uncertainty as shaded band
        self.ax_ratio.fill_between(self.bin_edges,
            1 - np.insert(self.data_MC_ratio['error_MC'], len(self.data_MC_ratio['error_MC']), 0),
            1 + np.insert(self.data_MC_ratio['error_MC'], len(self.data_MC_ratio['error_MC']), 0),
            step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0)
        self.ax_ratio.axhline(1, color='black', linestyle=':')
        self.ax_ratio.set_ylabel("Obs/Exp")
        # legends and labels
        self.fig.subplots_adjust(hspace=0.05)
        hep.cms.label(ax=self.ax, label="Preliminary", data=True, lumi=self.lumi, com=13.6, fontsize=16)
        handles, labels = self.ax.get_legend_handles_labels()
        self.ax.legend(handles[::-1], labels[::-1], loc='upper right', frameon=1, framealpha=1, bbox_to_anchor=(0.98, 0.98))
        self.ax.set_ylabel(f"Events/bin")
        self.ax_ratio.set_xlabel(self.variable_label)
        self.ax_ratio.set_ylim(ratio_min, ratio_max)
        self.ax.set_ylim(0, 1.4*np.max(self.stacked_block))
        self.ax.set_xlim(self.bin_edges[0], self.bin_edges[-1])
        self.ax.text(0.035, 0.925, self.channel_label, fontsize=18, fontweight="bold", transform=self.ax.transAxes)

        # Save to pdf and root
        plt.savefig(self.file_name.replace(".root", ".pdf"))
        print(f'Saved histogram to {self.file_name.replace(".root", ".pdf")}')
        plt.savefig(self.file_name.replace(".root", ".png"))
        print(f'Saved histogram to {self.file_name.replace(".root", ".png")}')







if __name__ == "__main__":

    histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/Tests1303/IPHC_FF_VVVLoose/Run3_2022/control/tt/datacard_m_vis_cp_inclusive_tt_Run3_2022.root", "tt_cp_inclusive", "tt", "Run3_2022", "m_vis", variable_label = r"m$_{vis}$ (GeV)")

    histo.plot_1D_histo()


