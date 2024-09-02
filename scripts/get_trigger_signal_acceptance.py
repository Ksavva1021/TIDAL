import ROOT
from collections import OrderedDict
from prettytable import PrettyTable

input_folder = '/vols/cms/ks1021/offline/HiggsDNA/IC/output/production/Run3_2022/'

VSjet_wp = 5
VSe_wp = 1
VSmu_wp = 3
lepton_iso = "0.15"

t_sel = f"idDeepTau2018v2p5VSjet_X >= {VSjet_wp} && idDeepTau2018v2p5VSmu_X >= {VSmu_wp}  && idDeepTau2018v2p5VSe_X >= {VSe_wp}"
e_sel = f"iso_X < {lepton_iso}"
m_sel = f"iso_X < {lepton_iso}"

charge_sel = "(os==1)"

baseline_sel = {
    "tt":"({sel_1} && {sel_2})".format(sel_1=t_sel.replace("X","1"), sel_2=t_sel.replace("X","2")),
    "mt":"({sel_1} && {sel_2})".format(sel_1=m_sel.replace("X","1"), sel_2=t_sel.replace("X","2")),
    "et":"({sel_1} && {sel_2})".format(sel_1=e_sel.replace("X","1"), sel_2=t_sel.replace("X","2")),
}

files = [
    "GluGluHTo2Tau_UncorrelatedDecay_Filtered",
    "VBFHToTauTau_UncorrelatedDecay_Filtered",
    "DYto2L_M-50_madgraphMLM"
]

test_dict = OrderedDict()

test_dict["tt"] = [
    "trg_doubletau",
    "trg_doubletauandjet",
    "trg_doubletauandjet_2",
    "trg_doubletau || trg_doubletauandjet",
    "trg_doubletau || trg_doubletauandjet_2",
    "trg_doubletau || trg_doubletauandjet || trg_doubletauandjet_2",
]

test_dict["mt"] = [
    "trg_singlemuon",
    "trg_mt_cross",
    "trg_singlemuon || trg_mt_cross",
]

test_dict["et"] = [
    "trg_singleelectron",
    "trg_et_cross",
    "trg_singleelectron || trg_et_cross",
]

output_scores = OrderedDict()

for file in files:
    output_scores[file] = OrderedDict()
    for k, v in test_dict.items():
        fr = ROOT.TFile(f"{input_folder}/{k}/{file}/nominal/merged.root")
        tr = fr.Get('ntuple')
        output_scores[file][k] = OrderedDict()
        for s in v:
            output_scores[file][k][s] = round(float(tr.GetEntries("(("+s+")&&("+baseline_sel[k]+")&&("+charge_sel+"))"))/float(tr.GetEntries("(("+baseline_sel[k]+")&&("+charge_sel+"))")),2)

def highlight_max(data):
    """Highlight the maximum value in each column."""
    max_values = [max(col) for col in zip(*data)]
    highlighted_data = []
    for row in data:
        highlighted_row = []
        for value, max_value in zip(row, max_values):
            if value == max_value:
                highlighted_row.append(f"\033[91m{value}\033[0m")  # Red color
            else:
                highlighted_row.append(str(value))
        highlighted_data.append(highlighted_row)
    return highlighted_data

t = PrettyTable(["Channel", "Selection"] + files)

for k, v in test_dict.items():
    data = []
    for i, s in enumerate(v):
        if i == 0:
            row = [k, s.replace(" ", "")]
        else:
            row = ["", s.replace(" ", "")]
        for f in files:
            row.append(output_scores[f][k][s])
        data.append(row)

    highlighted_data = highlight_max([row[2:] for row in data])
    for i, row in enumerate(data):
        t.add_row(row[:2] + highlighted_data[i])
    if v:
        t.add_row(["", "", "", "", ""], divider=True)

print(t)