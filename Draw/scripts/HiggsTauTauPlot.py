# HiggsTauTauPlot.py

# Methods:
# 1. Basic method using MC samples & SS method for QCD estimation
# 2. Basic method using MC samples & SS method for QCD estimation && WJets shape method (High mT control region)
# 3. ??

import argparse
from collections import OrderedDict
from prettytable import PrettyTable
import copy

import ROOT
from Draw.python import Analysis
from Draw.python import Plotting
from Draw.python.nodes import BuildCutString, GenerateZTT, GenerateZLL, GenerateTop, GenerateVV, GenerateW, GenerateQCD
from Draw.python.HiggsTauTauPlot_utilities import PrintSummary, GetTotals, FixBins

ROOT.TH1.SetDefaultSumw2(True)

parser = argparse.ArgumentParser()
parser.add_argument('--input_folder', type=str, help='Input folder')
parser.add_argument('--output_folder', default='output', help='Output folder')
parser.add_argument('--channel', default='mm', help='Channel to run on')
parser.add_argument('--era', default='2016', help='Era to run on')
parser.add_argument('--parameter_file', type=str, help='Parameter file')
parser.add_argument('--method', default=1, help='Method to run on')
parser.add_argument('--sel', type=str, help='Additional Selection to apply', default='')
parser.add_argument('--var', type=str, help='Variable to plot')
parser.add_argument('--do_ss', action='store_true', help='Do SS')
args = parser.parse_args()

available_channels = ['mm', 'em', 'mt', 'et', 'tt']

if args.channel not in available_channels:
    raise ValueError("Invalid channel. Please choose from: {}".format(available_channels))
available_methods = ["1","2","3"]
if args.method not in available_methods:
    raise ValueError("Invalid method. Please choose from: {}".format(available_methods))

table = PrettyTable()
table.field_names = ['Details', 'Choices']
table.add_row(['Input Folder', args.input_folder])
table.add_row(['Parameter File', args.parameter_file])
table.add_row(['Output Folder', args.output_folder])
table.add_row(['Channel', args.channel])
table.add_row(['Era', args.era])
table.add_row(['Method', args.method])
table.add_row(['Selection', args.sel])
table.add_row(['Variable', args.var])

method = int(args.method)

# ------------------------------------------------------------------------------------------------------------------------
# Define baseline selections and different categories
categories = {}
if args.era in ["Run3_2022"]:
    if args.channel == "mm":
        categories['baseline'] = '(iso_1<0.15 && iso_2<0.15 && (trg_singlemuon && abs(eta_1) < 2.1))'
    if args.channel == "mt":
        categories['baseline'] = '(iso_1 < 0.15 && idDeepTau2018v2p5VSjet_2 >= 5 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_2 >= 4 && (trg_singlemuon && pt_1 >= 25  && abs(eta_1) < 2.1))'
    if args.channel == "et":
        categories['baseline'] = '(iso_1 < 0.15 && idDeepTau2018v2p5VSjet_2 >= 5 && idDeepTau2018v2p5VSe_2 >= 6 && idDeepTau2018v2p5VSmu_2 >= 1 && (trg_singleelectron && pt_1 >= 25))'
    if args.channel == "tt":
        categories['baseline'] = '(idDeepTau2018v2p5VSjet_1 >= 5 && idDeepTau2018v2p5VSjet_2 >= 5 && idDeepTau2018v2p5VSe_1 >= 2 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_1 >= 3 && idDeepTau2018v2p5VSmu_2 >= 3 && (trg_doubletau && pt_1 > 40 && pt_2 > 40))'

categories['inclusive'] = '(1)'
categories['nobtag'] = '(n_bjets==0)'
categories['btag'] = '(n_bjets>=1)'
categories['w_sdb'] = 'mt_1>70.'
categories['w_shape'] = ''
categories['tt_qcd_norm'] = '(idDeepTau2018v2p5VSjet_1 < 5 && idDeepTau2018v2p5VSjet_1 >= 3 && idDeepTau2018v2p5VSjet_2 < 5 && idDeepTau2018v2p5VSjet_2 >= 3 && idDeepTau2018v2p5VSe_1 >=2 && idDeepTau2018v2p5VSe_2 >=2 && idDeepTau2018v2p5VSmu_1 >= 1 && idDeepTau2018v2p5VSmu_2 >= 1 && (trg_doubletau && pt_1 > 40 && pt_2 > 40))'

if args.channel == 'tt':
    categories["inclusive_pipi"]     = "(decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==0 && ip_LengthSig_2>=1.5)"
    categories["inclusive_pirho"]       = "((decayMode_1==1 && decayMode_2==0 && ip_LengthSig_2>=1.5) || (decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==1))"
    categories["inclusive_rhorho"]       = "(decayMode_1==1 && decayMode_2==1)"
    categories["inclusive_a1pi"]     = "((decayMode_1==10 && decayMode_2==0 && ip_LengthSig_2>=1.5) || (decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==10))"
    categories["inclusive_a1rho"]     = "((decayMode_1==10 && decayMode_2==1) || (decayMode_1==1 && decayMode_2==10))"
    categories["inclusive_a1a1"]     = "(decayMode_1==10 && decayMode_2==10)"
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define the samples (Data and MC (Background & Signal))
if args.era in ["Run3_2022"]:
    samples_dict= {}
    # Data Samples
    if args.channel == "et":
        data_samples = ['EGamma_Run2022C', 'EGamma_Run2022D']
    elif args.channel in ["mm","mt"]:
        data_samples = ['SingleMuon_Run2022C','Muon_Run2022C','Muon_Run2022D']
    elif args.channel == "tt":
        data_samples = ['Tau_Run2022C','Tau_Run2022D']

    samples_dict['data_samples'] = data_samples

    # MC Samples
    ztt_samples = ['DYto2L_M-50_madgraphMLM','DYto2L_M-50_madgraphMLM_ext1','DYto2L_M-50_1J_madgraphMLM','DYto2L_M-50_2J_madgraphMLM','DYto2L_M-50_3J_madgraphMLM','DYto2L_M-50_4J_madgraphMLM']
    top_samples = ['TTto2L2Nu','TTto2L2Nu_ext1','TTtoLNu2Q','TTtoLNu2Q_ext1','TTto4Q','TTto4Q_ext1']
    vv_samples = ['WW','WZ','ZZ','ST_t-channel_top_4f_InclusiveDecays','ST_t-channel_antitop_4f_InclusiveDecays','ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1','ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1','ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1','ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1']
    wjets_samples = ['WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1','WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM','WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM']
    if args.channel in ['et','mt','tt']:
        signal_samples = {
            "qqH_sm_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "qqH_ps_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "qqH_mm_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "ggH_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_ps_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_mm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "WH_sm_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "WH_ps_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "WH_mm_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "ZH_sm_htt*": 'ZHToTauTau_UncorrelatedDecay_Filtered',
            "ZH_ps_htt*": 'ZHToTauTau_UncorrelatedDecay_Filtered',
            "ZH_mm_htt*": 'ZHToTauTau_UncorrelatedDecay_Filtered',
        }
    else:
        signal_samples = {}

    samples_dict['ztt_samples'] = ztt_samples
    samples_dict['top_samples'] = top_samples
    samples_dict['vv_samples'] = vv_samples
    samples_dict['wjets_samples'] = wjets_samples
    samples_dict['signal_samples'] = signal_samples
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define gen selections to separate MC samples by Gen Flags

# Hint:
    # Muon Gen Matching:
    # 1 = prompt muon
    # 15 = muon from prompt tau
    # 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched

    # Electron Gen Matching:
    # 1 = prompt electron (including gamma*->mu mu)
    # 15 = electron from prompt tau
    # 22 = prompt photon (likely conversion)
    # 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched

    # Tau Gen Matching:
    # 1 = prompt electron
    # 2 = prompt muon
    # 3 = tau->e decay
    # 4 = tau->mu decay
    # 5 = hadronic tau decay
    # 0 = unknown or unmatched

gen_sels = {}
gen_sels_dict = {}
if args.channel == 'mm':

    gen_sels['ll_sel'] = '(genPartFlav_1==1 & genPartFlav_2==1)'
    gen_sels['tt_sel'] = '(genPartFlav_1==15 & genPartFlav_2==15)'
    gen_sels['j_sel'] = '(!(' + gen_sels['ll_sel'] + ') && !(' + gen_sels['tt_sel'] + '))'

    z_sels = {}
    z_sels['ztt_sel'] = gen_sels['tt_sel']
    z_sels['zl_sel'] = gen_sels['ll_sel']
    z_sels['zj_sel'] = gen_sels['j_sel']

    top_sels = {}
    top_sels['ttt_sel'] = gen_sels['ll_sel']
    top_sels['ttj_sel'] = '!(' + gen_sels['ll_sel'] + ')'

    vv_sels = {}
    vv_sels['vvt_sel'] = gen_sels['ll_sel']
    vv_sels['vvj_sel'] = '!(' + gen_sels['ll_sel'] + ')'

if args.channel in ['mt','et']:

    gen_sels['tt_sel'] = '(genPartFlav_1==15 && genPartFlav_2==5)'
    gen_sels['ll_sel'] = f"(genPartFlav_2!=0 && !{gen_sels['tt_sel']})"
    gen_sels['j_sel'] = '(genPartFlav_2==0)'

    z_sels = {}
    z_sels['ztt_sel'] = gen_sels['tt_sel']
    z_sels['zl_sel'] = gen_sels['ll_sel']
    z_sels['zj_sel'] = gen_sels['j_sel']

    top_sels = {}
    top_sels['ttt_sel'] = '!(' + gen_sels['j_sel'] + ')'
    top_sels['ttj_sel'] = gen_sels['j_sel']

    vv_sels = {}
    vv_sels['vvt_sel'] = '!(' + gen_sels['j_sel'] + ')'
    vv_sels['vvj_sel'] = gen_sels['j_sel']

if args.channel in ['tt']:

    gen_sels['tt_sel'] = '(genPartFlav_1==5 && genPartFlav_2==5)'
    gen_sels['ll_sel'] = f'(!(genPartFlav_1==0 || genPartFlav_2==0) && !{gen_sels["tt_sel"]})'
    gen_sels['j_sel'] = '(genPartFlav_1==0 || genPartFlav_2==0)'

    z_sels = {}
    z_sels['ztt_sel'] = gen_sels['tt_sel']
    z_sels['zl_sel'] = gen_sels['ll_sel']
    z_sels['zj_sel'] = gen_sels['j_sel']

    top_sels = {}
    top_sels['ttt_sel'] = '!(' + gen_sels['j_sel'] + ')'
    top_sels['ttj_sel'] = gen_sels['j_sel']

    vv_sels = {}
    vv_sels['vvt_sel'] = '!(' + gen_sels['j_sel'] + ')'
    vv_sels['vvj_sel'] = gen_sels['j_sel']

gen_sels_dict['z_sels'] = z_sels
gen_sels_dict['top_sels'] = top_sels
gen_sels_dict['vv_sels'] = vv_sels
# ------------------------------------------------------------------------------------------------------------------------

# RunPlotting handles how each process is added to the analysis

def RunPlotting(ana, nodename, samples_dict, gen_sels_dict, systematic='', cat='', cat_data='', categories={}, sel='', add_name='', wt='wt', do_data=True, qcd_factor=1.0, method=1):
    '''
    RunPlotting handles how each process is added to the analysis
    ana: Analysis object
    nodename: name of the node
    samples_dict: dictionary containing the samples
    gen_sels_dict: dictionary containing the gen selections
    systematic: systematic variation
    cat: category to be used for the selection
    cat_data: category to be used for the data selection
    categories: dictionary containing the categories
    sel: additional selection
    add_name: additional name to be added to the process (e.g. ZTT + xxx)
    wt: weight to be applied
    do_data: boolean to decide if data should be added
    qcd_factor: factor to be applied to QCD
    '''

    doZL = True
    doZJ = True
    doTTT = True
    doTTJ = True
    doVVT = True
    doVVJ = True

    if do_data:
        if args.do_ss:
          OSSS = '!os'
        else:
          OSSS = 'os'
        weight= ('weight')
        full_selection = BuildCutString(weight, sel, cat_data, OSSS)
        ana.nodes[nodename].AddNode(ana.SummedFactory('data_obs', data_samples, plot_unmodified, full_selection))

    GenerateZTT(ana, nodename, add_name, samples_dict['ztt_samples'], plot, wt, sel, cat, gen_sels_dict['z_sels'], not args.do_ss)
    GenerateZLL(ana, nodename, add_name, samples_dict['ztt_samples'], plot, wt, sel, cat, gen_sels_dict['z_sels'], not args.do_ss, doZL, doZJ)
    GenerateTop(ana, nodename, add_name, samples_dict['top_samples'], plot, wt, sel, cat, gen_sels_dict['top_sels'], not args.do_ss, doTTT, doTTJ)
    GenerateVV(ana, nodename, add_name, samples_dict['vv_samples'], plot, wt, sel, cat, gen_sels_dict['vv_sels'], not args.do_ss, doVVT, doVVJ)
    GenerateW(ana, nodename, add_name, samples_dict, gen_sels_dict, plot, plot_unmodified, wt, sel, cat, cat_data, categories, method=method, qcd_factor=qcd_factor, get_os=not args.do_ss)
    GenerateQCD(ana, nodename, add_name, samples_dict, gen_sels_dict, systematic, plot, plot_unmodified, wt, sel, cat, cat_data, categories=categories, method=method, qcd_factor=qcd_factor, get_os=not args.do_ss)
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Defining name of the output file
is_2d=False
is_3d=False
var_name = args.var.split('[')[0]
var_name = var_name.split('(')[0]
var_name = var_name.replace('/','_over_')
if var_name.count(',') == 1:
    is_2d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]
if var_name.count(',') == 2:
    is_3d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]+'_vs_'+var_name.split(',')[2]

output_name = f'{args.output_folder}/datacard_{var_name}_{args.channel}_{args.era}.root'
outfile = ROOT.TFile(output_name, 'RECREATE')
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define qcd factor, systematics
if args.channel == 'mm':
    qcd_factor = 1.07
elif args.channel == 'mt':
    qcd_factor = 1.12
elif args.channel == 'et':
    qcd_factor = 1.13
else:
    qcd_factor = 1.0

# set systematics:
# - 1st index sets folder name contaning systematic samples
# - 2nd index sets string to be appended to output histograms
# - 3rd index specifies the weight to be applied
# - 4th lists samples that should be skipped
nodename = args.channel
systematics = OrderedDict()
if args.channel == 'mt':
    systematics['nominal'] = ('nominal','','(w_Zpt_Reweighting*w_DY_soup*w_WJ_soup*w_Pileup*w_Muon_ID*w_Muon_Reco*w_Muon_Isolation*w_Tau_ID*w_Trigger)',[],False)
elif args.channel == 'et':
    systematics['nominal'] = ('nominal','','(w_Zpt_Reweighting*w_DY_soup*w_WJ_soup*w_Pileup*w_Electron_ID*w_Electron_Reco*w_Tau_ID*w_Trigger)',[],False)
elif args.channel == 'mm':
    systematics['nominal'] = ('nominal','','(w_Zpt_Reweighting*w_DY_soup*w_WJ_soup*w_Pileup*w_Muon_ID*w_Muon_Reco*w_Muon_Isolation*w_Trigger)',[],False)
elif args.channel == 'tt':
    systematics['nominal'] = ('nominal','','(w_Zpt_Reweighting*w_DY_soup*w_WJ_soup*w_Pileup*w_Tau_ID*w_Tau_e_FakeRate*w_Tau_mu_FakeRate*w_Trigger)',[],False)
# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
# Loop over systematics & run plotting, etc
samples_to_skip_dict = {}
systematic_suffixes = []
max_systematics_per_pass = 10

while len(systematics) > 0:
    analysis = Analysis.Analysis()
    analysis.nodes.AddNode(Analysis.ListNode(args.channel))
    analysis.remaps = {}

    if args.channel in ["mm","mt"]:
        analysis.remaps['Muon'] = 'data_obs'
    if args.channel == "tt":
        analysis.remaps['Tau'] = 'data_obs'

    previous_systematic_variation = None
    for index, systematic in enumerate(list(systematics.keys())[:max_systematics_per_pass]):
        if previous_systematic_variation is not None and systematics[systematic][0] != previous_systematic_variation:
            continue
        previous_systematic_variation = systematics[systematic][0]
        print("Processing:", systematic)
        print("")

        sel = args.sel
        plot = args.var
        plot_unmodified = plot
        categories_unmodified = copy.deepcopy(categories)
        systematic_folder_name = systematics[systematic][0]
        systematic_suffix = systematics[systematic][1]
        weight = systematics[systematic][2]
        samples_to_skip = systematics[systematic][3]
        is_FFsyst = systematics[systematic][4]

        systematic_suffixes.append(systematic_suffix)

        for sample_name in data_samples:
            analysis.AddSamples(f'{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root', 'ntuple', None, sample_name)

        for sample_name in ztt_samples + top_samples + vv_samples + wjets_samples:
            analysis.AddSamples(f'{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root', 'ntuple', None, sample_name)

        analysis.AddInfo(args.parameter_file, scaleTo='data_obs')

        if systematic == 'nominal':
            do_data = True
        else:
            do_data = False
        RunPlotting(analysis, nodename, samples_dict, gen_sels_dict, systematic, categories['baseline'], categories_unmodified['baseline'], categories_unmodified, sel, systematic_suffix, weight, do_data, qcd_factor, method)

        del systematics[systematic]

    analysis.Run()
    analysis.nodes.Output(outfile)

    FixBins(analysis, nodename, outfile)
    for suffix in systematic_suffixes:
      GetTotals(analysis, nodename, suffix, samples_dict, outfile)
    PrintSummary(analysis, nodename, ['data_obs'], add_names=systematic_suffixes, channel=args.channel, samples_dict=samples_dict)
# ------------------------------------------------------------------------------------------------------------------------

outfile.Close()
plot_file = ROOT.TFile(output_name, 'READ')
titles = Plotting.SetAxisTitles(args.var,args.channel)
x_title = titles[0]
y_title = titles[1]

Plotting.HTTPlot(
  nodename=nodename,
  infile=plot_file,
  channel=args.channel,
  scheme=args.channel,
  ratio_range="0.7,1.3",
  x_title=x_title,
  y_title=y_title,
  plot_name=output_name.replace('.root',''),
  lumi="Run3 2022 - 8.08 fb^{-1} (13.6 TeV)",
)
