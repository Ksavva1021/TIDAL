# HiggsTauTauPlot.py

# Methods:
# 1. Basic method using MC samples & SS method for QCD estimation
# 2. Basic method using MC samples & SS method for QCD estimation && WJets shape method (High mT control region)
# 3. ??

import argparse
from collections import OrderedDict
from prettytable import PrettyTable
import copy
import numpy as np
import ROOT
from Draw.python import Analysis
from Draw.python import Plotting
from Draw.python.nodes import BuildCutString, GenerateZTT, GenerateZLL, GenerateTop, GenerateVV, GenerateW, GenerateQCD, GenerateReweightedCPSignal
from Draw.python.HiggsTauTauPlot_utilities import PrintSummary, GetTotals, FixBins, FindRebinning, RebinHist

ROOT.TH1.SetDefaultSumw2(True)

parser = argparse.ArgumentParser()
parser.add_argument('--input_folder', type=str, help='Input folder')
parser.add_argument('--output_folder', default='output', help='Output folder')
parser.add_argument('--channel', default='mm', help='Channel to run on')
parser.add_argument('--era', default='2016', help='Era to run on')
parser.add_argument('--parameter_file', type=str, help='Parameter file')
parser.add_argument('--method', default=1, help='Method to run on')
parser.add_argument('--category', default='inclusive', help='Category to run on')
parser.add_argument('--sel', type=str, help='Additional Selection to apply', default='')
parser.add_argument('--add_weight', default='', help='Additional weight to apply')
parser.add_argument('--var', type=str, help='Variable to plot')
parser.add_argument('--do_ss', action='store_true', help='Do SS')
parser.add_argument('--blind', action='store_true', help='Blind the plot (remove data)')
parser.add_argument('--masses', default='125', help='Mass points to process, seperated by commas')
parser.add_argument('--datacard_name', help='Override the datacard name')
parser.add_argument('--auto_rebin', action='store_true', help='Automatically rebin histograms')
parser.add_argument('--LO_DY', action='store_true', help='Use LO instead of NLO DY')
args = parser.parse_args()

masses = args.masses.split(',')

available_channels = ['ee', 'mm', 'em', 'mt', 'et', 'tt']

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
#TODO: add option to change triggers
categories = {}
if args.era in ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]:
    if args.channel == "ee":
        categories['baseline'] = '(iso_1<0.15 && iso_2<0.15 && (trg_singleelectron && pt_1 > 31 && abs(eta_1) < 2.1))'
    if args.channel == "mm":
        categories['baseline'] = '(iso_1<0.15 && iso_2<0.15 && (trg_singlemuon && pt_1 > 26 && abs(eta_1) < 2.4))'
    if args.channel == "mt":
        mt_cross_only = '(trg_mt_cross && pt_1 > 21 && pt_1 <= 26 && abs(eta_1) < 2.1 && pt_2 > 32 && abs(eta_2) < 2.1)'
        single_muon_only = '(trg_singlemuon && pt_1 > 26  && abs(eta_1) < 2.4)'
        trg_full = '(%s || %s)' % (mt_cross_only, single_muon_only)
        categories['baseline'] = '(iso_1 < 0.15 && idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)' % trg_full
    if args.channel == "et":
        et_cross_only = '(trg_et_cross && pt_1 > 25 && pt_1 <= 31 && abs(eta_1) < 2.1 && pt_2 > 35 && abs(eta_2) < 2.1)'
        single_electron_only = '(trg_singleelectron && pt_1 > 31 && abs(eta_1) < 2.1 )'
        trg_full = '(%s || %s)' % (et_cross_only, single_electron_only)
        categories['baseline'] = '(iso_1 < 0.15&& idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)' % trg_full
    if args.channel == "tt":
        doubletau_only_trg = '(trg_doubletau && pt_1 > 40 && pt_2 > 40)'
        doubletaujet_only_trg = '(trg_doubletauandjet && pt_1 > 35 && pt_2 > 35 && jpt_1 > 60)' # might need to revise jet cut later on
        trg_full = '(%s || %s)' % (doubletau_only_trg, doubletaujet_only_trg)
        categories['baseline'] = '(abs(eta_1) < 2.1 && abs(eta_2) < 2.1 && idDeepTau2018v2p5VSjet_1 >= 7 && idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_1 >= 2 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_1 >= 4 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)' % trg_full
        categories['tt_qcd_norm'] = categories['baseline'].replace('idDeepTau2018v2p5VSjet_1 >= 7', 'idDeepTau2018v2p5VSjet_1 < 7 && idDeepTau2018v2p5VSjet_1 >= 3')

categories['inclusive'] = '(1)'
categories['nobtag'] = '(n_bjets==0)'
categories['btag'] = '(n_bjets>=1)'
categories['w_sdb'] = 'mt_1>70.'
categories['w_shape'] = ''
categories['aminus_low'] = '(alphaAngle_mu_pi_1 < {} && svfit_Mass < 100 && mt_1<50 && ip_LengthSig_1 > 1)'.format(np.pi/4)
categories['aminus_high'] = '(alphaAngle_mu_pi_1 > {} && svfit_Mass < 100 && mt_1<50 && ip_LengthSig_1 > 1)'.format(np.pi/4)
categories['ip_control'] = '(m_vis > 40 && m_vis < 90 && mt_1 < 40)'

categories['xt_dM0'] = '(decayMode_2 == 0)'
categories['xt_dM1'] = '(decayMode_2 == 1)'
categories['xt_dM10'] = '(decayMode_2 == 10)'
categories['xt_dM11'] = '(decayMode_2 == 11)'

if args.channel == 'tt':


    categories["inclusive_pipi"]     = "(decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==0 && ip_LengthSig_2>=1.5)"
    categories["inclusive_pirho"]       = "((decayMode_1==1 && decayMode_2==0 && ip_LengthSig_2>=1.5) || (decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==1))"
    categories["inclusive_rhorho"]       = "(decayMode_1==1 && decayMode_2==1)"
    categories["inclusive_a1pi"]     = "((decayMode_1==10 && hasRefitSV_1 && decayMode_2==0 && ip_LengthSig_2>=1.5) || (decayMode_1==0 && ip_LengthSig_1>=1.5 && decayMode_2==10 && hasRefitSV_2))"
    categories["inclusive_a1rho"]     = "((decayMode_1==10 && hasRefitSV_1 && decayMode_2==1) || (decayMode_1==1 && decayMode_2==10 && hasRefitSV_2))"
    categories["inclusive_a1a1"]     = "(decayMode_1==10 && decayMode_2==10 && hasRefitSV_1 && hasRefitSV_2)"

    sel_pi = 'decayModePNet_X==0 && ip_LengthSig_X>=1.25'
    sel_rho = 'decayMode_X==1 && decayModePNet_X==1 && pion_E_split_X>0.2'
    sel_a1 = 'decayModePNet_X==10'
    sel_a11pr = 'decayMode_X==1 && decayModePNet_X==2 && pion_E_split_X>0.2'

    sel_pi_1 = sel_pi.replace('X','1')
    sel_pi_2 = sel_pi.replace('X','2')
    sel_rho_1 = sel_rho.replace('X','1')
    sel_rho_2 = sel_rho.replace('X','2')
    sel_a1_1 = sel_a1.replace('X','1')
    sel_a1_2 = sel_a1.replace('X','2')
    sel_a11pr_1 = sel_a11pr.replace('X','1')
    sel_a11pr_2 = sel_a11pr.replace('X','2')

    categories["inclusive_PNet_rhorho"] = '(%(sel_rho_1)s && %(sel_rho_2)s)' % vars() 
    categories["inclusive_PNet_pipi"] = '(%(sel_pi_1)s && %(sel_pi_2)s)' % vars()
    categories["inclusive_PNet_a1a1"] = '(%(sel_a1_1)s && %(sel_a1_2)s)' % vars()
    categories["inclusive_PNet_rhoa11pr"] = '((%(sel_rho_1)s && %(sel_a11pr_2)s) || (%(sel_a11pr_1)s && %(sel_rho_2)s) || (%(sel_a11pr_1)s && %(sel_a11pr_2)s))' % vars()

    categories["inclusive_PNet_pirho"]     = '(%(sel_pi_1)s && %(sel_rho_2)s)' % vars()
    categories["inclusive_PNet_a1rho"]     = '(%(sel_a1_1)s && %(sel_rho_2)s)' % vars()
    categories["inclusive_PNet_a1pi"]      = '(%(sel_a1_1)s && %(sel_pi_2)s)' % vars()
    categories["inclusive_PNet_pia11pr"]   = '(%(sel_pi_1)s && %(sel_a11pr_2)s)' % vars()
    categories["inclusive_PNet_a1a11pr"]   = '(%(sel_a1_1)s && %(sel_a11pr_2)s)' % vars()

    categories["inclusive_PNet_rhopi"]     = '(%(sel_rho_1)s && %(sel_pi_2)s)' % vars()
    categories["inclusive_PNet_rhoa1"]     = '(%(sel_rho_1)s && %(sel_a1_2)s)' % vars()
    categories["inclusive_PNet_pia1"]      = '(%(sel_pi_1)s && %(sel_a1_2)s)' % vars()
    categories["inclusive_PNet_a11prpi"]   = '(%(sel_a11pr_1)s && %(sel_pi_2)s)' % vars()
    categories["inclusive_PNet_a11pra1"]   = '(%(sel_a11pr_1)s && %(sel_a1_2)s)' % vars()

    categories['mva_higgs'] = '(BDT_pred_class==1)'
    categories['mva_fake']  = '(BDT_pred_class==2)'
    categories['mva_tau']   = '(BDT_pred_class==0)'

    tt_channels = ['rhorho','pirho','rhopi','a1rho','rhoa1','a1pi','pia1','a1a1','pipi','pia11pr','a11prpi','rhoa11pr','a1a11pr','a11pra1']
    for c in tt_channels:
        categories["higgs_{}".format(c)] = '({} && {})'.format(categories['mva_higgs'], categories["inclusive_PNet_{}".format(c)])
        categories["tau_{}".format(c)] = '({} && {})'.format(categories['mva_tau'], categories["inclusive_PNet_{}".format(c)])
        categories["fake_{}".format(c)] = '({} && {})'.format(categories['mva_fake'], categories["inclusive_PNet_{}".format(c)])

# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define the samples (Data and MC (Background & Signal))
if args.era in ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]:
    samples_dict= {}
    # Data Samples
    if args.era in ["Run3_2022"]:
        if args.channel in ["ee","et"]:
            data_samples = ['EGamma_Run2022C', 'EGamma_Run2022D']
        elif args.channel in ["mm","mt"]:
            data_samples = ['SingleMuon_Run2022C','Muon_Run2022C','Muon_Run2022D']
        elif args.channel == "tt":
            data_samples = ['Tau_Run2022C','Tau_Run2022D']
    elif args.era in ["Run3_2022EE"]:
        if args.channel in ["ee","et"]:
            data_samples = ['EGamma_Run2022E', 'EGamma_Run2022F', 'EGamma_Run2022G']
        elif args.channel in ["mm","mt"]:
            data_samples = ['Muon_Run2022E','Muon_Run2022F','Muon_Run2022G']
        elif args.channel == "tt":
            data_samples = ['Tau_Run2022E','Tau_Run2022F','Tau_Run2022G']
    elif args.era in ["Run3_2023"]:
        if args.channel in ["ee","et"]:
            data_samples = ['EGamma0_Run2023C_v1', 'EGamma0_Run2023C_v2', 'EGamma0_Run2023C_v3', 'EGamma0_Run2023C_v4', 'EGamma1_Run2023C_v1', 'EGamma1_Run2023C_v2', 'EGamma1_Run2023C_v3', 'EGamma1_Run2023C_v4']
        elif args.channel in ["mm","mt"]:
            data_samples = ['Muon0_Run2023C_v1', 'Muon0_Run2023C_v2', 'Muon0_Run2023C_v3', 'Muon0_Run2023C_v4', 'Muon1_Run2023C_v1', 'Muon1_Run2023C_v2', 'Muon1_Run2023C_v3', 'Muon1_Run2023C_v4']
        elif args.channel == "tt":
            data_samples = ['Tau_Run2023C_v1', 'Tau_Run2023C_v2', 'Tau_Run2023C_v3', 'Tau_Run2023C_v4']
    elif args.era in ["Run3_2023BPix"]:
        if args.channel in ["ee","et"]:
            data_samples = ['EGamma0_Run2023D_v1', 'EGamma0_Run2023D_v2', 'EGamma1_Run2023D_v1', 'EGamma1_Run2023D_v2']
        if args.channel in ["mm","mt"]:
            data_samples = ['Muon0_Run2023D_v1', 'Muon0_Run2023D_v2', 'Muon1_Run2023D_v1', 'Muon1_Run2023D_v2']
        elif args.channel == "tt":
            data_samples = ['Tau_Run2023D_v1', 'Tau_Run2023D_v2']

    samples_dict['data_samples'] = data_samples

    # MC Samples
    if args.LO_DY:
        print('WARNING: Using LO DY samples')
        ztt_samples = ['DYto2L_M_50_madgraphMLM','DYto2L_M_50_madgraphMLM_ext1','DYto2L_M_50_1J_madgraphMLM','DYto2L_M_50_2J_madgraphMLM','DYto2L_M_50_3J_madgraphMLM','DYto2L_M_50_4J_madgraphMLM']
        if args.era in ["Run3_2023", "Run3_2023BPix"]:
            ztt_samples.remove('DYto2L_M_50_madgraphMLM_ext1')
    else:
        ztt_samples = ['DYto2L_M_50_amcatnloFXFX', 'DYto2L_M_50_amcatnloFXFX_ext1', 'DYto2L_M_50_0J_amcatnloFXFX', 'DYto2L_M_50_1J_amcatnloFXFX',
                    'DYto2L_M_50_2J_amcatnloFXFX', 'DYto2L_M_50_PTLL_40to100_1J_amcatnloFXFX', 'DYto2L_M_50_PTLL_100to200_1J_amcatnloFXFX',
                    'DYto2L_M_50_PTLL_200to400_1J_amcatnloFXFX', 'DYto2L_M_50_PTLL_400to600_1J_amcatnloFXFX', 'DYto2L_M_50_PTLL_600_1J_amcatnloFXFX',
                    'DYto2L_M_50_PTLL_40to100_2J_amcatnloFXFX', 'DYto2L_M_50_PTLL_100to200_2J_amcatnloFXFX', 'DYto2L_M_50_PTLL_200to400_2J_amcatnloFXFX',
                    'DYto2L_M_50_PTLL_400to600_2J_amcatnloFXFX', 'DYto2L_M_50_PTLL_600_2J_amcatnloFXFX'] # use NLO samples
        if args.era in ["Run3_2023", "Run3_2023BPix"]:
            ztt_samples.remove('DYto2L_M_50_amcatnloFXFX_ext1')
    top_samples = ['TTto2L2Nu','TTto2L2Nu_ext1','TTtoLNu2Q','TTtoLNu2Q_ext1','TTto4Q','TTto4Q_ext1']
    vv_samples = ['WW','WZ','ZZ','ST_t_channel_top_4f_InclusiveDecays','ST_t_channel_antitop_4f_InclusiveDecays','ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1','ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1','ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1','ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1']
    #vv_samples = ['WW','WZ','ZZ','ST_t-channel_top_4f_InclusiveDecays','ST_t-channel_antitop_4f_InclusiveDecays','ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1','ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1','ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1','ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1']
    wjets_samples = ['WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1','WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM','WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM']

    if args.era in ["Run3_2023", "Run3_2023BPix"]:
        top_samples.remove('TTto2L2Nu_ext1')
        top_samples.remove('TTtoLNu2Q_ext1')
        top_samples.remove('TTto4Q_ext1')
        vv_samples.remove('ST_tW_top_2L2Nu_ext1')
        vv_samples.remove('ST_tW_antitop_2L2Nu_ext1')
        vv_samples.remove('ST_tW_top_LNu2Q_ext1')
        if args.era == "Run3_2023":
            vv_samples.remove('ST_tW_antitop_LNu2Q')
        vv_samples.remove('ST_tW_antitop_LNu2Q_ext1')
        wjets_samples.remove('WtoLNu_madgraphMLM_ext1')

    if args.channel in ['et','mt','tt']:
        signal_samples = {
            "qqH_sm_unfiltered_htt*": 'VBFHToTauTau_UncorrelatedDecay_UnFiltered',
            "qqH_sm_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "qqH_ps_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "qqH_mm_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "Higgs_flat_htt*": ['VBFHToTauTau_UncorrelatedDecay_Filtered', 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
            "qqH_flat_htt*": 'VBFHToTauTau_UncorrelatedDecay_Filtered',
            "ggH_flat_prod_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_sm_prod_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_sm_prod_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_sm_prod_ps_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
            "ggH_sm_prod_mm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
            "ggH_ps_prod_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_ps_prod_ps_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
            "ggH_ps_prod_mm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
            "ggH_mm_prod_sm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "ggH_mm_prod_ps_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
            "ggH_mm_prod_mm_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
            #"ggH_sm_prod_sm_reweight_htt*": ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay', 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
            "WH_sm_unfiltered_htt*" : ['WplusHToTauTau_UncorrelatedDecay_UnFiltered','WminusHToTauTau_UncorrelatedDecay_UnFiltered'],
            "WH_sm_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "WH_ps_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "WH_mm_htt*" : ['WplusHToTauTau_UncorrelatedDecay_Filtered','WminusHToTauTau_UncorrelatedDecay_Filtered'],
            "ZH_sm_unfiltered_htt*": 'ZHToTauTau_UncorrelatedDecay_UnFiltered',
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
if args.channel in ['ee','mm']:

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

def Get1DBinNumFrom2D(h2d,xbin,ybin):
    Nxbins = h2d.GetNbinsX()
    return (ybin-1)*Nxbins + xbin -1

def UnrollHist2D(h2d,inc_y_of=True):
    '''
    Unroll a 2D histogram h2d into a 1d histogram
    inc_y_of = True includes the y over-flow bins
    '''
    if inc_y_of: n = 1
    else: n = 0
    Nbins = (h2d.GetNbinsY()+n)*(h2d.GetNbinsX())
    if isinstance(h2d, ROOT.TH2D): h1d = ROOT.TH1D(h2d.GetName(), '', Nbins, 0, Nbins)
    else: h1d = ROOT.TH1F(h2d.GetName(), '', Nbins, 0, Nbins)
    for i in range(1,h2d.GetNbinsX()+1):
      for j in range(1,h2d.GetNbinsY()+1+n):
        glob_bin = Get1DBinNumFrom2D(h2d,i,j)
        content = h2d.GetBinContent(i,j)
        error = h2d.GetBinError(i,j)
        h1d.SetBinContent(glob_bin+1,content)
        h1d.SetBinError(glob_bin+1,error)
    return h1d

# RunPlotting handles how each process is added to the analysis

def RunPlotting(ana, nodename, samples_dict, gen_sels_dict, systematic='', cat_name='', categories={}, categories_unmodified={}, sel='', add_name='', wt='wt', do_data=True, qcd_factor=1.0, method=1):
    '''
    RunPlotting handles how each process is added to the analysis
    ana: Analysis object
    nodename: name of the node
    samples_dict: dictionary containing the samples
    gen_sels_dict: dictionary containing the gen selections
    systematic: systematic variation
    cat_name: the name of the category to be used for the data selection
    categories: dictionary containing the categories
    categories_unmodified: dictionary containing the categories - unmodified by any systematic variations - to be applied for data only
    sel: additional selection
    add_name: additional name to be added to the process (e.g. ZTT + xxx)
    wt: weight to be applied
    do_data: boolean to decide if data should be added
    qcd_factor: factor to be applied to QCD
    '''

    cat = categories['cat']
    cat_data = categories_unmodified['cat']


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
    GenerateW(ana, nodename, add_name, samples_dict, gen_sels_dict, plot, plot_unmodified, wt, sel, cat_name, categories, categories_unmodified=categories_unmodified, method=method, qcd_factor=qcd_factor, get_os=not args.do_ss)
    GenerateQCD(ana, nodename, add_name, samples_dict, gen_sels_dict, systematic, plot, plot_unmodified, wt, sel, cat_name, categories=categories, categories_unmodified=categories_unmodified, method=method, qcd_factor=qcd_factor, get_os=not args.do_ss)

    # generate correct signal
    # TODO: add scheme or similar flat to determine which ones to use
    GenerateReweightedCPSignal(ana, nodename, add_name, samples_dict['signal_samples'], masses, plot, wt, sel, cat, not args.do_ss)
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

category_name = args.category
if args.datacard_name:
    output_name = f'{args.output_folder}/datacard_{args.datacard_name}_{category_name}_{args.channel}_{args.era}.root'
else:
    output_name = f'{args.output_folder}/datacard_{var_name}_{category_name}_{args.channel}_{args.era}.root'
if args.do_ss:
    output_name = output_name.replace('.root','_ss.root')
outfile = ROOT.TFile(output_name, 'RECREATE')
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define qcd factor, systematics
if args.channel in ['ee','mm']:
    qcd_factor = 1.07
elif args.channel == 'mt':
    qcd_factor = 1.12
elif args.channel == 'et':
    qcd_factor = 1.13
else:
    qcd_factor = 1.0

weight = 'weight'
if args.add_weight:
    weight += '*'+args.add_weight

# set systematics:
# - 1st index sets folder name contaning systematic samples
# - 2nd index sets string to be appended to output histograms
# - 3rd index specifies the weight to be applied
# - 4th lists samples that should be skipped
nodename = args.channel+'_'+args.category
systematics = OrderedDict()
if args.channel == 'mt':
    systematics['nominal'] = ('nominal','',f'({weight})',[],False)
elif args.channel == 'et':
    systematics['nominal'] = ('nominal','',f'({weight})',[],False)
elif args.channel in ['ee','mm']:
    systematics['nominal'] = ('nominal','',f'({weight})',[],False)
elif args.channel == 'tt':
    systematics['nominal'] = ('nominal','',f'({weight})',[],False)
# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
categories['cat'] = '('+categories[args.category]+')*('+categories['baseline']+')'

# Loop over systematics & run plotting, etc
samples_to_skip_dict = {}
systematic_suffixes = []
max_systematics_per_pass = 10

while len(systematics) > 0:
    analysis = Analysis.Analysis()
    analysis.nodes.AddNode(Analysis.ListNode(nodename))
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
        # use plot_unmodified and categories_unmodified in cases where the data and MC get different selections due to a systematic variation
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

        for key, value in signal_samples.items():
            if not isinstance(value, (list,)): value = [value]
            for samp in value:
                for mass in masses:
                    sample_name = samp.replace('*',mass)
                    analysis.AddSamples(f'{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root', 'ntuple', None, sample_name)


        analysis.AddInfo(args.parameter_file, scaleTo='data_obs')


        if systematic == 'nominal':
            do_data = True
        else:
            do_data = False
        RunPlotting(analysis, nodename, samples_dict, gen_sels_dict, systematic, args.category, categories, categories_unmodified, sel, systematic_suffix, weight, do_data, qcd_factor, method)

        del systematics[systematic]

    analysis.Run()
    analysis.nodes.Output(outfile)

    FixBins(analysis, nodename, outfile)
    for suffix in systematic_suffixes:
      GetTotals(analysis, nodename, suffix, samples_dict, outfile)
    PrintSummary(analysis, nodename, ['data_obs'], add_names=systematic_suffixes, channel=args.channel, samples_dict=samples_dict)
# ------------------------------------------------------------------------------------------------------------------------

# unroll 2D histograms into 1D histograms but store both versions

if is_2d:
  x_lines = []
  y_labels = []
  first_hist = True
  # loop over all TH2Ds and for each one unroll to produce TH1D and add to datacard
  directory = outfile.Get(nodename)
  outfile.cd(nodename)
  hists_to_add = []
  for key in directory.GetListOfKeys():
    hist_name = key.GetName()
    hist = directory.Get(hist_name).Clone()

    if not isinstance(hist,ROOT.TDirectory):
      include_of = False

      h1d = UnrollHist2D(hist,include_of)
      hists_to_add.append(h1d)
      if first_hist:
        first_hist=False
        Nxbins = hist.GetNbinsX()
        for i in range(1,hist.GetNbinsY()+1): x_lines.append(Nxbins*i)
        for j in range(1,hist.GetNbinsY()+1): y_labels.append([hist.GetYaxis().GetBinLowEdge(j),hist.GetYaxis().GetBinLowEdge(j+1)])
        if include_of: y_labels.append([hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY()+1),-1])
  for hist in hists_to_add:
      directory.Get(hist.GetName()).Write(hist.GetName()+'_2D') # write a copy of the 2D histogram as this will be overwritten by 1D version
      hist.Write("",ROOT.TObject.kOverwrite)

outfile.Close()
plot_file = ROOT.TFile(output_name, 'READ')

if args.auto_rebin:
    outfile_rebin = ROOT.TFile(output_name.replace(".root","_rebinned.root"), 'RECREATE')
    outfile_rebin.mkdir(nodename)
    outfile_rebin.cd(nodename)
    total_bkghist = plot_file.Get(nodename+'/total_bkg').Clone()
    binning = FindRebinning(total_bkghist,BinThreshold=100,BinUncertFraction=0.5)

    print("New binning:", binning)
    hists_done = []
    for i in  plot_file.Get(nodename).GetListOfKeys():
        if i.GetName() not in hists_done:
            if ".subnodes" not in i.GetName():
                RebinHist(plot_file.Get(nodename+'/'+i.GetName()).Clone(),binning).Write()
                hists_done.append(i.GetName())

    outfile_rebin.Close()
    plot_file = ROOT.TFile(output_name.replace(".root","_rebinned.root"))

titles = Plotting.SetAxisTitles(args.var,args.channel)
x_title = titles[0]
y_title = titles[1]

Plotting.HTTPlot(
  nodename=nodename,
  infile=plot_file,
  channel=args.channel,
  scheme=args.channel,
  ratio_range="0.5,1.5",
  x_title=x_title,
  y_title=y_title,
  plot_name=output_name.replace('.root',''),
  lumi=f"{args.era}",
  #lumi="Run3 2022 - 8.08 fb^{-1} (13.6 TeV)",
  blind=args.blind,
  log_y=False,
)
