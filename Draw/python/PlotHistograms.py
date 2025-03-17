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

        self.initialize_plotting()
        self.initialize_nodes()
        self.get_labels()


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



    def initialize_nodes(self):
        # get total background nodes
        if self.channel == "tt":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W", "ZL", "ZJ"], "color": "red"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "\\tau_1"
            self.lep2 = "\\tau_2"
        elif self.channel == "mt":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\mu\\mu$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "\\mu"
            self.lep2 = "\\tau"
        elif self.channel == "et":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to ee": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "e"
            self.lep2 = "\\tau"
        elif self.channel == "em":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\ell\\ell": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "e"
            self.lep2 = "\\mu"
        elif self.channel == "ee":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to ee": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "e_1"
            self.lep2 = "e_2"
        elif self.channel == "mm":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\mu\\mu$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "\\mu_1"
            self.lep2 = "\\mu_2"
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

    def get_labels(self):
        label_map = {
            "m_vis": r"m$_{vis}$ (GeV)",
            "eta_1": rf"$\eta^{{{self.lep1}}}$",
            "eta_2": rf"$\eta^{{{self.lep2}}}$",
            "phi_1": rf"$\phi^{{{self.lep1}}}$",
            "phi_2": rf"$\phi^{{{self.lep2}}}$",
            "pt_1": rf"$p_T^{{{self.lep1}}}$ (GeV)",
            "pt_2": rf"$p_T^{{{self.lep2}}}$ (GeV)",
            "met_pt": r"$p_T^{miss}$ (GeV)",
            "met_phi": r"$\phi^{miss}$",
            "met_dphi_1": rf"$\Delta\phi({self.lep1},$" + r"$p_T^{miss})$",
            "met_dphi_2": rf"$\Delta\phi({self.lep2},$" + r"$p_T^{miss})$",
            "dR": rf"$\Delta R({self.lep1}, {self.lep2})$",
            "dphi": rf"$\Delta\phi({self.lep1}, {self.lep2})$",
            "pt_tt": r"$p_T$" + rf"$^{{{self.lep1}}}$" +rf"$^{{{self.lep2}}}$ (GeV)",
            "pt_vis": r"$p_T^{vis}$ (GeV)",
            "FastMTT_mass": r"$m_{FastMTT}$ (GeV)",
            "FastMTT_mass_constraint": r"$m_{FastMTT}^{constr.}$ (GeV)",
            "mt_1": rf"$m_T({self.lep1},$" + r"$p_T^{miss})$",
            "mt_2": rf"$m_T({self.lep2},$" + r"$p_T^{miss})$",
            "mt_lep": rf"$m_T({self.lep1}, {self.lep2})$",
            "mt_tot": r"$m_T^{tot}$ (GeV)",
            "jpt_1": r"$p_T^{j_1}$ (GeV)",
            "jpt_2": r"$p_T^{j_2}$ (GeV)",
            "jeta_1": r"$\eta^{j_1}$",
            "jeta_2": r"$\eta^{j_2}$",
            "mjj": r"$m_{jj}$ (GeV)",
            "jdeta": r"$\Delta\eta_{jj}$",
            "dijetpt": r"$p_T^{jj}$ (GeV)",
            "n_jets": r"$N_{jets}$",
            "n_bjets": r"$N_{bjets}$",

        }
        self.variable_label = label_map.get(self.variable, self.variable)

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
        self.ax_ratio.axhline(1, color='black', linestyle=':')
        self.fig.subplots_adjust(hspace=0.05)
        # add data/MC ratio as points
        if not self.blind:
            self.ax_ratio.errorbar(self.bin_centers, self.data_MC_ratio['counts'], yerr=self.data_MC_ratio['error_data'], xerr=self.bin_widths/2, fmt='o', color = 'black', markersize=3, linewidth=0.6)
        # add MC uncertainty as shaded band
        self.ax_ratio.fill_between(self.bin_edges,
            1 - np.insert(self.data_MC_ratio['error_MC'], len(self.data_MC_ratio['error_MC']), 0),
            1 + np.insert(self.data_MC_ratio['error_MC'], len(self.data_MC_ratio['error_MC']), 0),
            step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0)

        # legends and labels
        hep.cms.label(ax=self.ax, label="Preliminary", data=True, lumi=self.lumi, com=13.6, fontsize=16)
        handles, labels = self.ax.get_legend_handles_labels()
        self.ax.text(0.035, 0.925, self.channel_label, fontsize=18, fontweight="bold", transform=self.ax.transAxes)
        self.ax.legend(handles[::-1], labels[::-1], loc='upper right', frameon=1, framealpha=1, bbox_to_anchor=(0.98, 0.98))

        # main plot
        if "(GeV)" in self.variable_label:
            self.ax.set_ylabel(f"Events / {self.bin_widths[0]} GeV")
        else:
            self.ax.set_ylabel(f"Events / {self.bin_widths[0]}")
        self.ax.set_ylim(0, 1.5*np.max(self.stacked_block))
        self.ax.set_xlim(self.bin_edges[0], self.bin_edges[-1])

        # ratio plot
        self.ax_ratio.set_ylabel("Obs/Exp")
        self.ax_ratio.set_xlabel(self.variable_label)
        self.ax_ratio.set_ylim(ratio_min, ratio_max)

        # Save to pdf and root
        plt.savefig(self.file_name.replace(".root", ".pdf"))
        print(f'Saved histogram to {self.file_name.replace(".root", ".pdf")}')
        plt.savefig(self.file_name.replace(".root", ".png"))
        print(f'Saved histogram to {self.file_name.replace(".root", ".png")}')







if __name__ == "__main__":

    histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/Tests1303/IPHC_FF_VVVLoose/Run3_2022/control/tt/datacard_m_vis_cp_inclusive_tt_Run3_2022.root", "tt_cp_inclusive", "tt", "Run3_2022", "m_vis", variable_label = r"m$_{vis}$ (GeV)")

    histo.plot_1D_histo()


