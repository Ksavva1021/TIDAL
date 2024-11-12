# Plotting.py

import ROOT as R
import math
import numpy as np
from array import array
from math import log10, floor
import re
import json
import types


COL_STORE = []


def SetAxisTitles(plot, channel):
    if "[" in plot:
        isVarBins = True
        var = plot.split("[")[0]
    else:
        isVarBins = False
        var = plot.split("(")[0]
        if len(plot.split("(")) > 2:
            var = "".join(plot.split("(")[:-1])

    chan_label = "#tau#tau"
    lep1_label = "#tau"
    lep2_label = "#tau"
    if channel == "et":
        chan_label = "e#tau"
        lep1_label = "e"
        lep2_label = "#tau"
    elif channel == "mt":
        chan_label = "#mu#tau"
        lep1_label = "#mu"
        lep2_label = "#tau"
    if channel == "em":
        chan_label = "e#mu"
        lep1_label = "e"
        lep2_label = "#mu"
    elif channel == "tt":
        chan_label = "#tau#tau"
        lep1_label = "#tau_{1}"
        lep2_label = "#tau_{2}"
    elif channel == "zee" or channel == "ee":
        chan_label = "ee"
        lep1_label = "e_{1}"
        lep2_label = "e_{2}"
    elif channel == "zmm" or channel == "mm":
        chan_label = "#mu#mu"
        lep1_label = "#mu_{1}"
        lep2_label = "#mu_{2}"

    bin_width = ""
    # if not isVarBins:
    #    binning = plot.split('(')[1].split(')')[0].split(',')
    #    binning = map(float,binning)
    #    bin_width = str(round((binning[2]-binning[1])/binning[0],1))

    titles = {}
    titles["iso_1"] = [
        "I^{" + lep1_label + "}_{rel}",
        "Events / " + bin_width + " GeV",
        "I^{" + lep1_label + "}_{rel}",
    ]
    titles["iso_2"] = [
        "I^{" + lep2_label + "}_{rel}",
        "Events / " + bin_width + " GeV",
        "I^{" + lep2_label + "}_{rel}",
    ]
    titles["iso_2_V2p5"] = [
        "I^{" + lep2_label + "}_{rel}",
        "Events / " + bin_width + " GeV",
        "I^{" + lep2_label + "}_{rel}",
    ]
    titles["pt_1"] = [
        "p_{T}^{" + lep1_label + "} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dp_{T}^{" + lep1_label + "} (1/GeV)",
    ]
    titles["pt_2"] = [
        "p_{T}^{" + lep2_label + "} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dp_{T}^{" + lep2_label + "} (1/GeV)",
    ]
    titles["met"] = [
        "E_{T}^{miss} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dE_{T}^{miss} (1/GeV)",
    ]
    titles["eta_1"] = [
        "#eta_{" + lep1_label + "}",
        "Events / " + bin_width,
        "dN/d#eta_{" + lep1_label + "}",
    ]
    titles["eta_2"] = [
        "#eta_{" + lep2_label + "}",
        "Events / " + bin_width,
        "dN/d#eta_{" + lep2_label + "}",
    ]
    titles["mt_tot"] = [
        "m_{T}^{tot} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}^{tot} (1/GeV)",
    ]
    titles["mt_1"] = [
        "m_{T}(p_{T}^{" + lep1_label + "},p_{T}^{miss}) (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}(p_{T}^{" + lep1_label + "},p_{T}^{miss}) (1/GeV)",
    ]
    titles["mt_2"] = [
        "m_{T}(p_{T}^{" + lep2_label + "},p_{T}^{miss}) (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}(p_{T}^{" + lep2_label + "},p_{T}^{miss}) (1/GeV)",
    ]
    titles["mt_lep"] = [
        "m_{T}(p_{T}^{" + lep1_label + "},p_{T}^{" + lep2_label + "}) (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}(p_{T}^{" + lep1_label + "},p_{T}^{" + lep2_label + "}) (1/GeV)",
    ]
    titles["m_vis"] = [
        "m_{" + chan_label + "} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{" + chan_label + "} (1/GeV)",
    ]
    titles["m_sv"] = [
        "m_{#tau#tau} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{#tau#tau} (1/GeV)",
    ]
    titles["svfit_mass"] = [
        "m_{#tau#tau} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{#tau#tau} (1/GeV)",
    ]
    titles["mjj"] = [
        "m_{jj} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{jj} (1/GeV)",
    ]
    titles["dphi"] = [
        "#Delta#phi(" + lep1_label + "," + lep2_label + ")",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
        "",
    ]
    titles["met_dphi_1"] = [
        "#Delta#phi(" + lep1_label + ",E_{T}^{miss})",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
    ]
    titles["met_dphi_2"] = [
        "#Delta#phi(" + lep2_label + ",E_{T}^{miss})",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
    ]
    titles["pt_1+pt_2+jpt_1+met"] = [
        "S_{T}^{MET} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dS_{T}^{MET} (1/GeV)",
    ]  # Excess
    if channel in ["zee", "zmm"]:
        titles["pt_tt"] = [
            "p_{T}^{" + chan_label + "} (GeV)",
            "Events / " + bin_width + " GeV",
            "dN/dp_{T}^{" + chan_label + "} (1/GeV)",
        ]
    else:
        titles["pt_tt"] = [
            "p_{T}^{#tau#tau} (GeV)",
            "Events / " + bin_width + " GeV",
            "dN/dp_{#tau#tau}^{tot} (1/GeV)",
        ]
    titles["n_jets"] = ["N_{jets}", "Events", "dN/dN_{jets}"]
    titles["n_deepbjets"] = ["N_{b-jets}", "Events", "dN/dN_{b-jets}"]
    titles["n_prebjets"] = ["N_{pre b-jets}", "Events", "dN/dN_{pre b-jets}"]
    titles["n_bjets"] = ["N_{b-jets}", "Events", "dN/dN_{b-jets}"]
    titles["n_btag"] = ["N_{b-tag}^{tight}", "Events", "dN/dN_{b-tag}^{tight}"]
    titles["n_loose_btag"] = ["N_{b-tag}^{loose}", "Events", "dN/dN_{b-tag}^{loose}"]
    titles["pzeta"] = [
        "D_{#zeta} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dD_{#zeta} (1/GeV)",
    ]
    titles["jdeta"] = [
        "|#Delta#eta_{jj}|",
        "Events / " + bin_width,
        "dN/d|#Delta#eta_{jj}|",
    ]
    titles["sjdphi"] = [
        "#Delta#phi_{jj}",
        "Events / " + bin_width,
        "dN/d#Delta#phi_{jj}",
    ]
    titles["D0"] = ["D_{0}", "Events / " + bin_width, "dN/dD_{0}"]
    titles["DCP"] = ["D_{CP}", "Events / " + bin_width, "dN/dD_{CP}"]
    titles["D0star"] = ["D_{0}^{*}", "Events / " + bin_width, "dN/dD_{0}^{*}"]
    titles["jeta_1"] = ["#eta_{j_{1}}", "Events / " + bin_width, "dN/d#eta_{j_{1}}"]
    titles["jeta_2"] = ["#eta_{j_{2}}", "Events / " + bin_width, "dN/d#eta_{j_{2}}"]
    titles["jpt_1"] = [
        "P_{T}^{j_{1}} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dP_{T}^{j_{1}} (1/GeV)",
    ]
    titles["jpt_2"] = [
        "P_{T}^{j_{2}} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dP_{T}^{j_{2}} (1/GeV)",
    ]
    titles["IC_lowMjj_Oct05_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["IC_highMjj_Oct05_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["aco_angle_mod"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_1"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_2"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_3"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_4"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_5"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["aco_angle_6"] = [
        "#phi#mbox{*}_{CP}",
        "Events / " + bin_width,
        "dN/d#phi#mbox{*}_{CP}",
    ]
    titles["IC_Feb13_fix1_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["IC_Mar26_fix2_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["IC_Apr02_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["IC_Vienna_fix_max_score"] = ["NN Score", "Events", "dN/d(NN Score)"]
    titles["IC_keras_sm4_max_score"] = ["NN Score", "Events", "dN/d(NN Score)"]
    titles["IC_keras_sm5_max_score"] = ["NN Score", "Events", "dN/d(NN Score)"]
    titles["IC_keras_sm6_max_score"] = ["NN Score", "Events", "dN/d(NN Score)"]
    titles["IC_Nov25_tauspinner_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)"]
    titles["tau_decay_mode_2"] = ["#tau decay mode", "Events", "Events"]
    # hacky for now
    titles[
        "(dphi_jtt<0.)*(dphi_jtt+2*3.14159265359)+(dphi_jtt>0)*(dphi_jtt)-3.14159265359)"
    ] = ["#Delta#phi_{jet, Z} - #pi", "Events", "Events"]
    if channel == "tt":
        titles["tau_decay_mode_1"] = ["Lead #tau decay mode", "Events", "Events"]
        titles["tau_decay_mode_2"] = ["Sub-lead #tau decay mode", "Events", "Events"]

    if var not in titles:
        if not isVarBins:
            return [var, "Events"]
        else:
            return [var, "dN/d" + var]
    else:
        if not isVarBins:
            return [titles[var][0], titles[var][1]]
        else:
            return [titles[var][0], titles[var][2]]


def SetAxisTitles2D(plot, channel):
    if "[" in plot:
        var = plot.split("[")[0]
    else:
        var = plot.split("(")[0]

    isVarBins = "[" in plot
    yvar = var.split(",")[0]
    xvar = var.split(",")[1]

    chan_label = "#tau#tau"
    lep1_label = "#tau"
    lep2_label = "#tau"
    if channel == "et":
        chan_label = "e#tau"
        lep1_label = "e"
        lep2_label = "#tau"
    elif channel == "mt":
        chan_label = "#mu#tau"
        lep1_label = "#mu"
        lep2_label = "#tau"
    if channel == "em":
        chan_label = "e#mu"
        lep1_label = "e"
        lep2_label = "#mu"
    elif channel == "tt":
        chan_label = "#tau#tau"
        lep1_label = "#tau_{1}"
        lep2_label = "#tau_{2}"
    elif channel == "zee":
        chan_label = "ee"
        lep1_label = "e_{1}"
        lep2_label = "e_{2}"
    elif channel == "ee":
        chan_label = "ee"
        lep1_label = "e_{1}"
        lep2_label = "e_{2}"
    elif channel == "zmm":
        chan_label = "#mu#mu"
        lep1_label = "#mu_{1}"
        lep2_label = "#mu_{2}"
    elif channel == "mm":
        chan_label = "#mu#mu"
        lep1_label = "#mu_{1}"
        lep2_label = "#mu_{2}"

    bin_width = ""
    if not isVarBins:
        binning = plot.split("(")[1].split(")")[0].split(",")
        binning = map(float, binning)
        bin_width = str(round((binning[2] - binning[1]) / binning[0], 1))

    titles = {}
    titles["pt_1"] = [
        "p_{T}^{" + lep1_label + "} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dP_{T}^{" + lep1_label + "} (1/GeV)",
        "GeV",
    ]
    titles["pt_2"] = [
        "p_{T}^{" + lep2_label + "} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dP_{T}^{" + lep2_label + "} (1/GeV)",
        "GeV",
    ]
    titles["met"] = [
        "E_{T}^{miss} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dE_{T}^{miss} (1/GeV)",
        "GeV",
    ]
    titles["eta_1"] = [
        "#eta_{" + lep1_label + "}",
        "Events / " + bin_width,
        "dN/d#eta_{" + lep1_label + "}",
        "",
    ]
    titles["eta_2"] = [
        "#eta_{" + lep2_label + "}",
        "Events / " + bin_width,
        "dN/d#eta_{" + lep2_label + "}",
        "",
    ]
    titles["mt_tot"] = [
        "M_{T}^{tot} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dM_{T}^{tot} (1/GeV)",
        "GeV",
    ]
    titles["mt_1"] = [
        "m_{T} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T} (1/GeV)",
        "GeV",
    ]
    titles["m_vis"] = [
        "m_{" + chan_label + "}^{vis} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{" + chan_label + "}^{vis} (1/GeV)",
        "GeV",
    ]
    titles["m_sv"] = [
        "m_{#tau#tau} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{#tau#tau} (1/GeV)",
        "GeV",
    ]
    titles["svfit_mass"] = [
        "m_{#tau#tau} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{#tau#tau} (1/GeV)",
        "GeV",
    ]
    titles["mjj"] = [
        "m_{jj} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{jj} (1/GeV)",
        "GeV",
    ]
    titles["tau_decay_mode_2"] = ["#tau decay mode", "Events", "Events", ""]
    titles["mva_dm_2"] = ["#tau MVA decay mode", "Events", "Events", ""]
    titles["mt_1"] = [
        "m_{T}^{2} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}^{2} (1/GeV)",
        "GeV",
    ]
    titles["mt_lep"] = [
        "m_{T}^{lep} (GeV)",
        "Events / " + bin_width + " GeV",
        "dN/dm_{T}^{lep} (1/GeV)",
        "GeV",
    ]
    titles["dphi"] = [
        "#Delta#phi(lep1,lep2)",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
        "",
    ]
    titles["met_dphi_1"] = [
        "#Delta#phi(lep1,E_{T}^{miss})",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
        "",
    ]
    titles["met_dphi_2"] = [
        "#Delta#phi(lep2,E_{T}^{miss})",
        "Events / " + bin_width,
        "dN/d#Delta#phi",
        "",
    ]

    if channel == "tt":
        titles["tau_decay_mode_1"] = ["Lead #tau decay mode", "Events", "Events", ""]
        titles["tau_decay_mode_2"] = [
            "Sub-lead #tau decay mode",
            "Events",
            "Events",
            "",
        ]
        titles["mva_dm_1"] = ["Lead MVA #tau MVA decay mode", "Events", "Events", ""]
        titles["mva_dm_2"] = ["Sub-lead #tau MVA decay mode", "Events", "Events", ""]

    titles["sjdphi"] = ["#Delta#phi_{jj}", "Events", "dN/d#Delta#phi_{jj}", ""]
    titles["D0"] = ["D_{0}", "Events", "dN/dD_{0}", ""]
    titles["DCP"] = ["D_{CP}", "Events", "dN/dD_{CP}", ""]
    titles["D0star"] = ["D_{0}^{*}", "Events", "dN/dD_{0}^{*}", ""]

    if channel in ["zee", "zmm"]:
        titles["pt_tt"] = [
            "p_{T}^{" + chan_label + "} (GeV)",
            "Events / " + bin_width + " GeV",
            "dN/dp_{T}^{" + chan_label + "} (1/GeV)",
            "GeV",
        ]
    else:
        titles["pt_tt"] = [
            "p_{T}^{#tau#tau} (GeV)",
            "Events / " + bin_width + " GeV",
            "dN/dp_{T}^{#tau#tau} (1/GeV)",
            "GeV",
        ]
    titles["n_jets"] = ["N_{jets}", "Events", "dN/dN_{jets}", ""]
    titles["n_bjets"] = ["N_{b-jets}", "Events", "dN/dN_{b-jets}", ""]
    titles["n_deepbjets"] = ["N_{b-jets}", "Events", "dN/dN_{b-jets}", ""]
    titles["IC_lowMjj_Sep25_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)", ""]
    titles["IC_highMjj_Oct05_max_score"] = [
        "MVA Score",
        "Events",
        "dN/d(MVA Score)",
        "",
    ]
    titles["aco_angle_mod"] = [
        "#phi#mbox{*}_{CP}",
        "Events",
        "dN/d#phi#mbox{*}_{CP}",
        "",
    ]
    titles["aco_angle_1"] = ["#phi#mbox{*}_{CP}", "Events", "dN/d#phi#mbox{*}_{CP}", ""]
    titles["IC_Feb13_fix1_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)", ""]
    titles["IC_Mar26_fix2_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)", ""]
    titles["IC_Apr02_max_score"] = ["MVA Score", "Events", "dN/d(MVA Score)", ""]
    titles["IC_Vienna_fix_max_score"] = ["NN Score", "Events", "dN/d(NN Score)", ""]
    titles["IC_keras_sm4_max_score"] = ["NN Score", "Events", "dN/d(NN Score)", ""]
    titles["IC_keras_sm5_max_score"] = ["NN Score", "Events", "dN/d(NN Score)", ""]
    titles["IC_keras_sm6_max_score"] = ["NN Score", "Events", "dN/d(NN Score)", ""]
    titles["IC_Nov13_tauspinner_v1_max_score"] = [
        "MVA Score",
        "Events",
        "dN/d(MVA Score)",
        "",
    ]

    if xvar not in titles:
        if not isVarBins:
            x_titles = [xvar, "Events"]
        else:
            x_titles = [xvar, "dN/d" + xvar]
    else:
        if not isVarBins:
            x_titles = [titles[xvar][0], titles[xvar][1]]
        else:
            x_titles = [titles[xvar][0], titles[xvar][2]]
    if yvar not in titles:
        unit = ""
        if not isVarBins:
            y_titles = [yvar, "Events", unit]
        else:
            y_titles = [yvar, "dN/d" + yvar, unit]
    else:
        unit = titles[yvar][3]
        if not isVarBins:
            y_titles = [titles[yvar][0], titles[yvar][1], unit]
        else:
            y_titles = [titles[yvar][0], titles[yvar][2], unit]

    return [x_titles, y_titles]


def SetTDRStyle():
    """Sets the PubComm recommended style

    Just a copy of <http://ghm.web.cern.ch/ghm/plots/MacroExample/tdrstyle.C>
    @sa ModTDRStyle() to use this style with some additional customisation.
    """
    # For the canvas:
    R.gStyle.SetCanvasBorderMode(0)
    R.gStyle.SetCanvasColor(R.kWhite)
    R.gStyle.SetCanvasDefH(600)  # Height of canvas
    R.gStyle.SetCanvasDefW(600)  # Width of canvas
    R.gStyle.SetCanvasDefX(0)  # POsition on screen
    R.gStyle.SetCanvasDefY(0)

    # For the Pad:
    R.gStyle.SetPadBorderMode(0)
    # R.gStyle.SetPadBorderSize(Width_t size = 1)
    R.gStyle.SetPadColor(R.kWhite)
    R.gStyle.SetPadGridX(False)
    R.gStyle.SetPadGridY(False)
    R.gStyle.SetGridColor(0)
    R.gStyle.SetGridStyle(3)
    R.gStyle.SetGridWidth(1)

    # For the frame:
    R.gStyle.SetFrameBorderMode(0)
    R.gStyle.SetFrameBorderSize(1)
    R.gStyle.SetFrameFillColor(0)
    R.gStyle.SetFrameFillStyle(0)
    R.gStyle.SetFrameLineColor(1)
    R.gStyle.SetFrameLineStyle(1)
    R.gStyle.SetFrameLineWidth(1)

    # For the histo:
    # R.gStyle.SetHistFillColor(1)
    # R.gStyle.SetHistFillStyle(0)
    R.gStyle.SetHistLineColor(1)
    R.gStyle.SetHistLineStyle(0)
    R.gStyle.SetHistLineWidth(1)
    # R.gStyle.SetLegoInnerR(Float_t rad = 0.5)
    # R.gStyle.SetNumberContours(Int_t number = 20)

    R.gStyle.SetEndErrorSize(2)
    # R.gStyle.SetErrorMarker(20)
    # R.gStyle.SetErrorX(0.)

    R.gStyle.SetMarkerStyle(20)

    # For the fit/function:
    R.gStyle.SetOptFit(1)
    R.gStyle.SetFitFormat("5.4g")
    R.gStyle.SetFuncColor(2)
    R.gStyle.SetFuncStyle(1)
    R.gStyle.SetFuncWidth(1)

    # For the date:
    R.gStyle.SetOptDate(0)
    # R.gStyle.SetDateX(Float_t x = 0.01)
    # R.gStyle.SetDateY(Float_t y = 0.01)

    # For the statistics box:
    R.gStyle.SetOptFile(0)
    R.gStyle.SetOptStat(0)
    # To display the mean and RMS:   SetOptStat('mr')
    R.gStyle.SetStatColor(R.kWhite)
    R.gStyle.SetStatFont(42)
    R.gStyle.SetStatFontSize(0.025)
    R.gStyle.SetStatTextColor(1)
    R.gStyle.SetStatFormat("6.4g")
    R.gStyle.SetStatBorderSize(1)
    R.gStyle.SetStatH(0.1)
    R.gStyle.SetStatW(0.15)
    # R.gStyle.SetStatStyle(Style_t style = 1001)
    # R.gStyle.SetStatX(Float_t x = 0)
    # R.gStyle.SetStatY(Float_t y = 0)

    # Margins:
    R.gStyle.SetPadTopMargin(0.05)
    R.gStyle.SetPadBottomMargin(0.13)
    R.gStyle.SetPadLeftMargin(0.16)
    R.gStyle.SetPadRightMargin(0.02)

    # For the Global title:
    R.gStyle.SetOptTitle(0)
    R.gStyle.SetTitleFont(42)
    R.gStyle.SetTitleColor(1)
    R.gStyle.SetTitleTextColor(1)
    R.gStyle.SetTitleFillColor(10)
    R.gStyle.SetTitleFontSize(0.05)
    # R.gStyle.SetTitleH(0); # Set the height of the title box
    # R.gStyle.SetTitleW(0); # Set the width of the title box
    # R.gStyle.SetTitleX(0); # Set the position of the title box
    # R.gStyle.SetTitleY(0.985); # Set the position of the title box
    # R.gStyle.SetTitleStyle(Style_t style = 1001)
    # R.gStyle.SetTitleBorderSize(2)

    # For the axis titles:
    R.gStyle.SetTitleColor(1, "XYZ")
    R.gStyle.SetTitleFont(42, "XYZ")
    R.gStyle.SetTitleSize(0.06, "XYZ")
    # Another way to set the size?
    # R.gStyle.SetTitleXSize(Float_t size = 0.02)
    # R.gStyle.SetTitleYSize(Float_t size = 0.02)
    R.gStyle.SetTitleXOffset(0.9)
    R.gStyle.SetTitleYOffset(1.25)
    # R.gStyle.SetTitleOffset(1.1, 'Y'); # Another way to set the Offset

    # For the axis labels:

    R.gStyle.SetLabelColor(1, "XYZ")
    R.gStyle.SetLabelFont(42, "XYZ")
    R.gStyle.SetLabelOffset(0.007, "XYZ")
    R.gStyle.SetLabelSize(0.05, "XYZ")

    # For the axis:

    R.gStyle.SetAxisColor(1, "XYZ")
    R.gStyle.SetStripDecimals(True)
    R.gStyle.SetTickLength(0.03, "XYZ")
    R.gStyle.SetNdivisions(510, "XYZ")
    R.gStyle.SetPadTickX(1)
    R.gStyle.SetPadTickY(1)

    # Change for log plots:
    R.gStyle.SetOptLogx(0)
    R.gStyle.SetOptLogy(0)
    R.gStyle.SetOptLogz(0)

    # Postscript options:
    R.gStyle.SetPaperSize(20.0, 20.0)
    # R.gStyle.SetLineScalePS(Float_t scale = 3)
    # R.gStyle.SetLineStyleString(Int_t i, const char* text)
    # R.gStyle.SetHeaderPS(const char* header)
    # R.gStyle.SetTitlePS(const char* pstitle)

    # R.gStyle.SetBarOffset(Float_t baroff = 0.5)
    # R.gStyle.SetBarWidth(Float_t barwidth = 0.5)
    # R.gStyle.SetPaintTextFormat(const char* format = 'g')
    # R.gStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
    # R.gStyle.SetTimeOffset(Double_t toffset)
    # R.gStyle.SetHistMinimumZero(kTRUE)

    R.gStyle.SetHatchesLineWidth(5)
    R.gStyle.SetHatchesSpacing(0.05)


def ModTDRStyle(width=600, height=600, t=0.06, b=0.12, l=0.16, r=0.04):
    """Modified version of the tdrStyle

    Args:
        width (int): Canvas width in pixels
        height (int): Canvas height in pixels
        t (float): Pad top margin [0-1]
        b (float): Pad bottom margin [0-1]
        l (float): Pad left margin [0-1]
        r (float): Pad right margin [0-1]
    """
    SetTDRStyle()

    # Set the default canvas width and height in pixels
    R.gStyle.SetCanvasDefW(width)
    R.gStyle.SetCanvasDefH(height)

    # Set the default margins. These are given as fractions of the pad height
    # for `Top` and `Bottom` and the pad width for `Left` and `Right`. But we
    # want to specify all of these as fractions of the shortest length.
    def_w = float(R.gStyle.GetCanvasDefW())
    def_h = float(R.gStyle.GetCanvasDefH())

    scale_h = (def_w / def_h) if (def_h > def_w) else 1.0
    scale_w = (def_h / def_w) if (def_w > def_h) else 1.0

    def_min = def_h if (def_h < def_w) else def_w

    R.gStyle.SetPadTopMargin(t * scale_h)
    # default 0.05
    R.gStyle.SetPadBottomMargin(b * scale_h)
    # default 0.13
    R.gStyle.SetPadLeftMargin(l * scale_w)
    # default 0.16
    R.gStyle.SetPadRightMargin(r * scale_w)
    # default 0.02
    # But note the new CMS style sets these:
    # 0.08, 0.12, 0.12, 0.04

    # Set number of axis tick divisions
    R.gStyle.SetNdivisions(510, "XYZ")  # default 510

    # Some marker properties not set in the default tdr style
    R.gStyle.SetMarkerColor(R.kBlack)
    R.gStyle.SetMarkerSize(1.0)

    R.gStyle.SetLabelOffset(0.007, "YZ")
    # This is an adhoc adjustment to scale the x-axis label
    # offset when we stretch plot vertically
    # Will also need to increase if first x-axis label has more than one digit
    R.gStyle.SetLabelOffset(0.005 * (3.0 - 2.0 / scale_h), "X")

    # In this next part we do a slightly involved calculation to set the axis
    # title offsets, depending on the values of the TPad dimensions and
    # margins. This is to try and ensure that regardless of how these pad
    # values are set, the axis titles will be located towards the edges of the
    # canvas and not get pushed off the edge - which can often happen if a
    # fixed value is used.
    title_size = 0.05
    title_px = title_size * def_min
    label_size = 0.04
    R.gStyle.SetTitleSize(title_size, "XYZ")
    R.gStyle.SetLabelSize(label_size, "XYZ")

    R.gStyle.SetTitleXOffset(
        0.5 * scale_h * (1.2 * (def_h * b * scale_h - 0.6 * title_px)) / title_px
    )
    R.gStyle.SetTitleYOffset(
        0.5 * scale_w * (1.2 * (def_w * l * scale_w - 0.6 * title_px)) / title_px
    )

    # Only draw ticks where we have an axis
    R.gStyle.SetPadTickX(0)
    R.gStyle.SetPadTickY(0)
    R.gStyle.SetTickLength(0.02, "XYZ")

    R.gStyle.SetLegendBorderSize(0)
    R.gStyle.SetLegendFont(42)
    R.gStyle.SetLegendFillColor(0)
    R.gStyle.SetFillColor(0)

    R.gROOT.ForceStyle()


def SetBirdPalette():
    nRGBs = 9
    stops = array(
        "d", [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]
    )
    red = array(
        "d", [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    )
    green = array(
        "d", [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    )
    blue = array(
        "d", [0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    )
    R.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, 1)


def SetDeepSeaPalette():
    nRGBs = 9
    stops = array(
        "d", [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]
    )
    red = array(
        "d",
        reversed(
            [
                0.0 / 255.0,
                9.0 / 255.0,
                13.0 / 255.0,
                17.0 / 255.0,
                24.0 / 255.0,
                32.0 / 255.0,
                27.0 / 255.0,
                25.0 / 255.0,
                29.0 / 255.0,
            ]
        ),
    )
    green = array(
        "d",
        reversed(
            [
                0.0 / 255.0,
                0.0 / 255.0,
                0.0 / 255.0,
                2.0 / 255.0,
                37.0 / 255.0,
                74.0 / 255.0,
                113.0 / 255.0,
                160.0 / 255.0,
                221.0 / 255.0,
            ]
        ),
    )
    blue = array(
        "d",
        reversed(
            [
                28.0 / 255.0,
                42.0 / 255.0,
                59.0 / 255.0,
                78.0 / 255.0,
                98.0 / 255.0,
                129.0 / 255.0,
                154.0 / 255.0,
                184.0 / 255.0,
                221.0 / 255.0,
            ]
        ),
    )
    R.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, 1)


def SetCorrMatrixPalette():
    R.TColor.CreateGradientColorTable(
        3,
        array("d", [0.00, 0.50, 1.00]),
        array("d", [1.00, 1.00, 0.00]),
        array("d", [0.70, 1.00, 0.34]),
        array("d", [0.00, 1.00, 0.82]),
        255,
        1.0,
    )


def CreateTransparentColor(color, alpha):
    adapt = R.gROOT.GetColor(color)
    new_idx = R.gROOT.GetListOfColors().GetLast() + 1
    trans = R.TColor(
        new_idx, adapt.GetRed(), adapt.GetGreen(), adapt.GetBlue(), "", alpha
    )
    COL_STORE.append(trans)
    trans.SetName("userColor%i" % new_idx)
    return new_idx


def Set(obj, **kwargs):
    for key, value in kwargs.iteritems():
        if value is None:
            getattr(obj, "Set" + key)()
        elif isinstance(value, (list, tuple)):
            getattr(obj, "Set" + key)(*value)
        else:
            getattr(obj, "Set" + key)(value)


def OnePad():
    pad = R.TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
    pad.SetTicks(1)
    pad.Draw()
    pad.cd()
    result = [pad]
    return result


def TwoPadSplit(split_point, gap_low, gap_high):
    upper = R.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
    upper.SetBottomMargin(split_point + gap_high)
    upper.SetFillStyle(4000)
    upper.SetTicks(1)
    upper.Draw()
    lower = R.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
    lower.SetTopMargin(1 - split_point + gap_low)
    lower.SetFillStyle(4000)
    lower.Draw()
    upper.cd()
    result = [upper, lower]
    return result


def MultiRatioSplit(split_points, gaps_low, gaps_high):
    """Create a set of TPads split vertically on the TCanvas

    This is a generalisation of the two pad main/ratio split but for the case
    of multiple ratio pads.

    Args:

        split_points (list[float]): Height of each ratio pad as a fraction of the
        canvas height. Pads will be created from the bottom of the frame
        upwards. The final, main pad will occupy however much space remains,
        therefore the size of this list should be [number of pads] - 1.
        gaps_low (list[float]): Gaps between ratio pad frames created on the
        lower pad side at each boundary. Give  a list of zeroes for no gap
        between pad frames. Should be the same length as `split_points`.1
        gaps_high (list[float]): Gaps between ratio pad frames created on the
        upper pad side at each boundary. Give a list of zeroes for no gap
        between pad frames.

    Returns:
        list[TPad]: List of TPads, indexed from top to bottom on the canvas.
    """
    pads = []
    for i in range(len(split_points) + 1):
        pad = R.TPad("pad%i" % i, "", 0.0, 0.0, 1.0, 1.0)
        if i > 0:
            pad.SetBottomMargin(sum(split_points[0:i]) + gaps_high[i - 1])
        if i < len(split_points):
            pad.SetTopMargin(1.0 - sum(split_points[0 : i + 1]) + gaps_low[i])
        pad.SetFillStyle(4000)
        pad.Draw()
        pads.append(pad)
    pads.reverse()
    return pads


def TwoPadSplitColumns(split_point, gap_left, gap_right):
    left = R.TPad("left", "left", 0.0, 0.0, 1.0, 1.0)
    left.SetRightMargin(1 - split_point + gap_right)
    left.SetFillStyle(4000)
    left.Draw()
    right = R.TPad("right", "right", 0.0, 0.0, 1.0, 1.0)
    right.SetLeftMargin(split_point + gap_left)
    right.SetFillStyle(4000)
    right.Draw()
    left.cd()
    result = [left, right]
    return result


def SetupTwoPadSplitAsRatio(pads, upper, lower, y_title, y_centered, y_min, y_max):
    if lower.GetXaxis().GetTitle() == "":
        lower.GetXaxis().SetTitle(upper.GetXaxis().GetTitle())
    upper.GetXaxis().SetTitle("")
    upper.GetXaxis().SetLabelSize(0)
    upper_h = 1.0 - pads[0].GetTopMargin() - pads[0].GetBottomMargin()
    lower_h = 1.0 - pads[1].GetTopMargin() - pads[1].GetBottomMargin()
    lower.GetYaxis().SetTickLength(R.gStyle.GetTickLength() * upper_h / lower_h)
    pads[1].SetTickx(1)
    pads[1].SetTicky(1)
    lower.GetYaxis().SetTitle(y_title)
    lower.GetYaxis().CenterTitle(y_centered)
    if y_max > y_min:
        lower.SetMinimum(y_min)
        lower.SetMaximum(y_max)


def CreateAxisHist(src, at_limits=True):
    backup = R.gPad
    tmp = R.TCanvas()
    tmp.cd()
    src.Draw("AP")
    result = src.GetHistogram().Clone("tmp")
    if at_limits:
        min = 0.0
        max = 0.0
        x = R.Double(0.0)
        y = R.Double(0.0)
        src.GetPoint(0, x, y)
        min = float(x)
        max = float(x)
        for i in range(1, src.GetN()):
            src.GetPoint(i, x, y)
            if x < min:
                min = float(x)
            if x > max:
                max = float(x)
        result.GetXaxis().SetLimits(min, max)
    R.gPad = backup
    return result


def CreateAxisHists(n, src, at_limits):
    res = []
    h = CreateAxisHist(src, at_limits)
    for i in range(n):
        res.append(h.Clone("tmp%i" % i))
    return res


def GetAxisHist(pad):
    pad_obs = pad.GetListOfPrimitives()
    if pad_obs is None:
        return None
    obj = None
    for obj in pad_obs:
        if obj.InheritsFrom(R.TH1.Class()):
            return obj
        if obj.InheritsFrom(R.TMultiGraph.Class()):
            return obj.GetHistogram()
        if obj.InheritsFrom(R.TGraph.Class()):
            return obj.GetHistogram()
        if obj.InheritsFrom(R.THStack.Class()):
            return obj.GetHistogram()
    return None


def TFileIsGood(filename):
    """Performs a series of tests on a TFile to ensure that it can be opened
    without errors

    Args:
        filename: `str` The name of the TFile to check

    Returns:
        `bool` True if the file can opened, is not a zombie, and if ROOT did
        not need to try and recover the contents
    """
    fin = R.TFile(filename)
    if not fin:
        return False
    if fin and not fin.IsOpen():
        return False
    elif fin and fin.IsOpen() and fin.IsZombie():
        fin.Close()
        return False
    elif fin and fin.IsOpen() and fin.TestBit(R.TFile.kRecovered):
        fin.Close()
        # don't consider a recovered file to be OK
        return False
    else:
        fin.Close()
        return True


def MakeTChain(files, tree):
    chain = R.TChain(tree)
    for f in files:
        chain.Add(f)
    return chain


def Get(file, obj):
    R.TH1.AddDirectory(False)
    f_in = R.TFile(file)
    res = R.gDirectory.Get(obj)
    f_in.Close()
    return res


def ParamFromFilename(filename, param):
    if len(re.findall(param + "\.\d+\.\d+", filename)):
        num1 = re.findall(param + "\.\d+\.\d+", filename)[0].replace(param + ".", "")
        return float(num1)
    elif len(re.findall(param + "\.\d+", filename)):
        num1 = re.findall(param + "\.\d+", filename)[0].replace(param + ".", "")
        return int(num1)
    else:
        print("Error: parameter " + param + " not found in filename")


def TGraphFromTree(tree, xvar, yvar, selection):
    tree.Draw(xvar + ":" + yvar, selection, "goff")
    gr = R.TGraph(tree.GetSelectedRows(), tree.GetV1(), tree.GetV2())
    return gr


def TGraph2DFromTree(tree, xvar, yvar, zvar, selection):
    tree.Draw(xvar + ":" + yvar + ":" + zvar, selection, "goff")
    gr = R.TGraph2D(tree.GetSelectedRows(), tree.GetV1(), tree.GetV2(), tree.GetV3())
    return gr


def RocCurveFrom1DHists(h_x, h_y, cut_is_greater_than):
    backup = R.TH1.AddDirectoryStatus()
    R.TH1.AddDirectory(False)
    x_den = h_x.Clone()
    x_num = h_x.Clone()
    x_err = R.Double(0.0)
    x_int = h_x.IntegralAndError(0, h_x.GetNbinsX() + 1, x_err)
    for i in range(1, h_x.GetNbinsX() + 1):
        x_part_err = R.Double(0.0)
        x_part_int = (
            h_x.IntegralAndError(i, h_x.GetNbinsX() + 1, x_part_err)
            if cut_is_greater_than
            else h_x.IntegralAndError(0, i, x_part_err)
        )
        x_den.SetBinContent(i, x_int)
        x_den.SetBinError(i, x_err)
        x_num.SetBinContent(i, x_part_int)
        x_num.SetBinError(i, x_part_err)
    y_den = h_y.Clone()
    y_num = h_y.Clone()
    y_err = R.Double(0.0)
    y_int = h_y.IntegralAndError(0, h_y.GetNbinsX() + 1, y_err)
    for i in range(1, h_y.GetNbinsX() + 1):
        y_part_err = R.Double(0.0)
        y_part_int = (
            h_y.IntegralAndError(i, h_y.GetNbinsX() + 1, y_part_err)
            if cut_is_greater_than
            else h_y.IntegralAndError(0, i, y_part_err)
        )
        y_den.SetBinContent(i, y_int)
        y_den.SetBinError(i, y_err)
        y_num.SetBinContent(i, y_part_int)
        y_num.SetBinError(i, y_part_err)
    # x_den.Print('all')
    # x_num.Print('all')
    # y_den.Print('all')
    # y_num.Print('all')
    x_gr = R.TGraphAsymmErrors(x_num, x_den)
    y_gr = R.TGraphAsymmErrors(y_num, y_den)

    res = y_gr.Clone()
    for i in range(0, res.GetN()):
        res.GetX()[i] = x_gr.GetY()[i]
        res.GetEXlow()[i] = x_gr.GetEYlow()[i]
        res.GetEXhigh()[i] = x_gr.GetEYhigh()[i]
    res.Sort()
    R.TH1.AddDirectory(backup)
    return res


def TH2FromTGraph2D(
    graph, method="BinEdgeAligned", force_x_width=None, force_y_width=None
):
    """Build an empty TH2 from the set of points in a TGraph2D

    There is no unique way to define a TH2 binning given an arbitrary
    TGraph2D, therefore this function supports multiple named methods:

     - `BinEdgeAligned` simply takes the sets of x- and y- values in the
       TGraph2D and uses these as the bin edge arrays in the TH2. The
       implication of this is that when filling the bin contents interpolation
       will be required when evaluating the TGraph2D at the bin centres.
     - `BinCenterAligned` will try to have the TGraph2D points at the bin
       centers, but this will only work completely correctly when the input
       point spacing is regular. The algorithm first identifies the bin width
       as the smallest interval between points on each axis. The start
       position of the TH2 axis is then defined as the lowest value in the
       TGraph2D minus half this width, and the axis continues with regular
       bins until the graph maximum is passed.

    Args:
        graph (TGraph2D): Should have at least two unique x and y values,
            otherwise we can't define any bins
        method (str): The binning algorithm to use
        force_x_width (bool): Override the derived x-axis bin width in the
            CenterAligned method
        force_y_width (bool): Override the derived y-axis bin width in the
            CenterAligned method

    Raises:
        RuntimeError: If the method name is not recognised

    Returns:
        TH2F: The exact binning of the TH2F depends on the chosen method
    """
    x_vals = set()
    y_vals = set()

    for i in range(graph.GetN()):
        x_vals.add(graph.GetX()[i])
        y_vals.add(graph.GetY()[i])

    x_vals = sorted(x_vals)
    y_vals = sorted(y_vals)
    if method == "BinEdgeAligned":
        h_proto = R.TH2F(
            "prototype",
            "",
            len(x_vals) - 1,
            array("d", x_vals),
            len(y_vals) - 1,
            array("d", y_vals),
        )
    elif method == "BinCenterAligned":
        x_widths = []
        y_widths = []
        for i in range(1, len(x_vals)):
            x_widths.append(x_vals[i] - x_vals[i - 1])
        for i in range(1, len(y_vals)):
            y_widths.append(y_vals[i] - y_vals[i - 1])
        x_min = min(x_widths) if force_x_width is None else force_x_width
        y_min = min(y_widths) if force_y_width is None else force_y_width
        x_bins = int(((x_vals[-1] - (x_vals[0] - 0.5 * x_min)) / x_min) + 0.5)
        y_bins = int(((y_vals[-1] - (y_vals[0] - 0.5 * y_min)) / y_min) + 0.5)
        print(
            "[TH2FromTGraph2D] x-axis binning: (%i, %g, %g)"
            % (
                x_bins,
                x_vals[0] - 0.5 * x_min,
                x_vals[0] - 0.5 * x_min + x_bins * x_min,
            )
        )
        print(
            "[TH2FromTGraph2D] y-axis binning: (%i, %g, %g)"
            % (
                y_bins,
                y_vals[0] - 0.5 * y_min,
                y_vals[0] - 0.5 * y_min + y_bins * y_min,
            )
        )
        # Use a number slightly smaller than 0.49999 because the TGraph2D interpolation
        # is fussy about evaluating on the boundary
        h_proto = R.TH2F(
            "prototype",
            "",
            x_bins,
            x_vals[0] - 0.49999 * x_min,
            x_vals[0] - 0.50001 * x_min + x_bins * x_min,
            y_bins,
            y_vals[0] - 0.49999 * y_min,
            y_vals[0] - 0.50001 * y_min + y_bins * y_min,
        )
    else:
        raise RuntimeError("[TH2FromTGraph2D] Method %s not supported" % method)
    h_proto.SetDirectory(0)
    return h_proto


def MakeErrorBand(LowerGraph, UpperGraph):
    errorBand = R.TGraphAsymmErrors()
    lower_list = []
    upper_list = []
    for i in range(LowerGraph.GetN()):
        lower_list.append((float(LowerGraph.GetX()[i]), float(LowerGraph.GetY()[i])))
        upper_list.append((float(UpperGraph.GetX()[i]), float(UpperGraph.GetY()[i])))
    lower_list = sorted(set(lower_list))
    upper_list = sorted(set(upper_list))
    for i in range(LowerGraph.GetN()):
        errorBand.SetPoint(i, lower_list[i][0], lower_list[i][1])
        errorBand.SetPointEYlow(i, lower_list[i][1] - lower_list[i][1])
        errorBand.SetPointEYhigh(i, upper_list[i][1] - lower_list[i][1])
    return errorBand


def LimitTGraphFromJSON(js, label):
    xvals = []
    yvals = []
    for key in js:
        xvals.append(float(key))
        yvals.append(js[key][label])
    graph = R.TGraph(len(xvals), array("d", xvals), array("d", yvals))
    graph.Sort()
    return graph


def LimitTGraphFromJSONFile(jsfile, label):
    with open(jsfile) as jsonfile:
        js = json.load(jsonfile)
    return LimitTGraphFromJSON(js, label)


def ToyTGraphFromJSON(js, label):
    xvals = []
    yvals = []
    if isinstance(label, types.StringTypes):
        for entry in js[label]:
            xvals.append(float(entry))
            yvals.append(1.0)
    else:
        if len(label) == 1:
            return ToyTGraphFromJSON(js, label[0])
        else:
            return ToyTGraphFromJSON(js[label[0]], label[1:])
    graph = R.TGraph(len(xvals), array("d", xvals), array("d", yvals))
    graph.Sort()
    return graph
    # hist = R.TH1F("toy", "toy", 100, min(xvals), max(xvals))
    # for xval in xvals:
    # hist.AddBinContent(hist.GetXaxis().FindBin(xval))
    # return hist


def ToyTGraphFromJSONFile(jsfile, label):
    with open(jsfile) as jsonfile:
        js = json.load(jsonfile)
    return ToyTGraphFromJSON(js, label)


def LimitBandTGraphFromJSON(js, central, lo, hi):
    xvals = []
    yvals = []
    yvals_lo = []
    yvals_hi = []
    for key in js:
        xvals.append(float(key))
        yvals.append(js[key][central])
        yvals_lo.append(js[key][central] - js[key][lo])
        yvals_hi.append(js[key][hi] - js[key][central])
    graph = R.TGraphAsymmErrors(
        len(xvals),
        array("d", xvals),
        array("d", yvals),
        array("d", [0]),
        array("d", [0]),
        array("d", yvals_lo),
        array("d", yvals_hi),
    )
    graph.Sort()
    return graph


def StandardLimitsFromJSONFile(json_file, draw=["obs", "exp0", "exp1", "exp2"]):
    graphs = {}
    data = {}
    with open(json_file) as jsonfile:
        data = json.load(jsonfile)
    if "obs" in draw:
        graphs["obs"] = LimitTGraphFromJSON(data, "obs")
    if "exp0" in draw or "exp" in draw:
        graphs["exp0"] = LimitTGraphFromJSON(data, "exp0")
    if "exp1" in draw or "exp" in draw:
        graphs["exp1"] = LimitBandTGraphFromJSON(data, "exp0", "exp-1", "exp+1")
    if "exp2" in draw or "exp" in draw:
        graphs["exp2"] = LimitBandTGraphFromJSON(data, "exp0", "exp-2", "exp+2")
    return graphs


def bestFit(tree, x, y, cut):
    nfind = tree.Draw(y + ":" + x, cut + "deltaNLL == 0")
    gr0 = R.TGraph(1)
    if nfind == 0:
        gr0.SetPoint(0, -999, -999)
    else:
        grc = R.gROOT.FindObject("Graph").Clone()
        if grc.GetN() > 1:
            grc.Set(1)
        gr0.SetPoint(0, grc.GetXmax(), grc.GetYmax())
    gr0.SetMarkerStyle(34)
    gr0.SetMarkerSize(2.0)
    return gr0


def treeToHist2D(t, x, y, name, cut, xmin, xmax, ymin, ymax, xbins, ybins):
    t.Draw(
        "2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)"
        % (y, x, name, xbins, xmin, xmax, ybins, ymin, ymax),
        cut + "deltaNLL != 0",
        "PROF",
    )
    prof = R.gROOT.FindObject(name + "_prof")
    h2d = R.TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax)
    for ix in range(1, xbins + 1):
        for iy in range(1, ybins + 1):
            z = prof.GetBinContent(ix, iy)
            if (z != z) or (z > 4294967295):  # protect against NANs
                z = 0
            h2d.SetBinContent(ix, iy, z)
    h2d.GetXaxis().SetTitle(x)
    h2d.GetYaxis().SetTitle(y)
    h2d.SetDirectory(0)
    h2d = NewInterpolate(h2d)
    return h2d


def makeHist1D(name, xbins, graph):
    len_x = graph.GetX()[graph.GetN() - 1] - graph.GetX()[0]
    binw_x = (len_x * 0.5 / (float(xbins) - 1.0)) - 1e-5
    hist = R.TH1F(
        name, "", xbins, graph.GetX()[0], graph.GetX()[graph.GetN() - 1] + binw_x
    )
    return hist


def makeHist2D(name, xbins, ybins, graph2d):
    len_x = graph2d.GetXmax() - graph2d.GetXmin()
    binw_x = (len_x * 0.5 / (float(xbins) - 1.0)) - 1e-5
    len_y = graph2d.GetYmax() - graph2d.GetYmin()
    binw_y = (len_y * 0.5 / (float(ybins) - 1.0)) - 1e-5
    hist = R.TH2F(
        name,
        "",
        xbins,
        graph2d.GetXmin() - binw_x,
        graph2d.GetXmax() + binw_x,
        ybins,
        graph2d.GetYmin() - binw_y,
        graph2d.GetYmax() + binw_y,
    )
    return hist


def makeVarBinHist2D(name, xbins, ybins):
    # create new arrays in which bin low edge is adjusted to make measured
    # points at the bin centres
    xbins_new = [None] * (len(xbins) + 1)
    for i in range(len(xbins) - 1):
        if i == 0 or i == 1:
            xbins_new[i] = xbins[i] - ((xbins[i + 1] - xbins[i]) / 2) + 1e-5
        else:
            xbins_new[i] = xbins[i] - ((xbins[i + 1] - xbins[i]) / 2)
    xbins_new[len(xbins) - 1] = xbins[len(xbins) - 2] + (
        (xbins[len(xbins) - 2] - xbins[len(xbins) - 3]) / 2
    )
    xbins_new[len(xbins)] = (
        xbins[len(xbins) - 1]
        + ((xbins[len(xbins) - 1] - xbins[len(xbins) - 2]) / 2)
        - 1e-5
    )

    ybins_new = [None] * (len(ybins) + 1)
    for i in range(len(ybins) - 1):
        if i == 0 or i == 1:
            ybins_new[i] = ybins[i] - ((ybins[i + 1] - ybins[i]) / 2) + 1e-5
        else:
            ybins_new[i] = ybins[i] - ((ybins[i + 1] - ybins[i]) / 2)
    ybins_new[len(ybins) - 1] = ybins[len(ybins) - 2] + (
        (ybins[len(ybins) - 2] - ybins[len(ybins) - 3]) / 2
    )
    ybins_new[len(ybins)] = (
        ybins[len(ybins) - 1]
        + ((ybins[len(ybins) - 1] - ybins[len(ybins) - 2]) / 2)
        - 1e-5
    )
    hist = R.TH2F(
        name,
        "",
        len(xbins_new) - 1,
        array("d", xbins_new),
        len(ybins_new) - 1,
        array("d", ybins_new),
    )
    return hist


def GraphDifference(graph1, graph2, relative):
    xvals = []
    yvals = []
    if graph1.GetN() != graph2.GetN():
        return graph1
    for i in range(graph1.GetN()):
        xvals.append(graph1.GetX()[i])
        if relative:
            yvals.append(
                2
                * abs(graph1.GetY()[i] - graph2.GetY()[i])
                / (graph1.GetY()[i] + graph2.GetY()[i])
            )
        else:
            yvals.append(
                2
                * (graph1.GetY()[i] - graph2.GetY()[i])
                / (graph1.GetY()[i] + graph2.GetY()[i])
            )
    diff_graph = R.TGraph(len(xvals), array("d", xvals), array("d", yvals))
    diff_graph.Sort()
    return diff_graph


def GraphDivide(num, den):
    res = num.Clone()
    for i in range(num.GetN()):
        if den.Eval(res.GetX()[i]) == 0:
            res.GetY()[i] = 0
        else:
            res.GetY()[i] = res.GetY()[i] / den.Eval(res.GetX()[i])
    if type(res) is R.TGraphAsymmErrors:
        for i in range(num.GetN()):
            if den.Eval(res.GetX()[i]) == 0:
                res.GetEYhigh()[i] = 0
                res.GetEYlow()[i] = 0
            else:
                res.GetEYhigh()[i] = res.GetEYhigh()[i] / den.Eval(res.GetX()[i])

    return res


def GraphDivideErrors(num, den):
    res = num.Clone()
    for i in range(num.GetN()):
        if type(res) is R.TGraphAsymmErrors:
            if den.Eval(res.GetX()[i]) == 0:
                res.GetEYhigh()[i] = 0
                res.GetEYlow()[i] = 0
            else:
                if res.GetY()[i] < 1e-100 or den.GetY()[i] < 1e-100:
                    res.GetEYhigh()[i] = 0
                    res.GetEYlow()[i] = 0
                else:
                    res.GetEYhigh()[i] = math.sqrt(
                        (res.GetEYhigh()[i] / res.GetY()[i]) ** 2
                        + (den.GetEYhigh()[i] / den.GetY()[i]) ** 2
                    )
                    res.GetEYlow()[i] = math.sqrt(
                        (res.GetEYlow()[i] / res.GetY()[i]) ** 2
                        + (den.GetEYlow()[i] / den.GetY()[i]) ** 2
                    )
        if den.Eval(res.GetX()[i]) == 0:
            res.GetY()[i] = 0
        else:
            res.GetY()[i] = res.GetY()[i] / den.Eval(res.GetX()[i])
    return res


def MakeRatioHist(num, den, num_err, den_err, ZeroFix=False):
    """Make a new ratio TH1 from numerator and denominator TH1s with optional
    error propagation

    Args:
        num (TH1): Numerator histogram
        den (TH1): Denominator histogram
        num_err (bool): Propagate the error in the numerator TH1
        den_err (bool): Propagate the error in the denominator TH1

    Returns:
        TH1: A new TH1 containing the ratio
    """
    result = num.Clone()
    if not num_err:
        for i in range(1, result.GetNbinsX() + 1):
            result.SetBinError(i, 0.0)
    den_fix = den.Clone()
    if not den_err:
        for i in range(1, den_fix.GetNbinsX() + 1):
            den_fix.SetBinError(i, 0.0)
    result.Divide(den_fix)
    if ZeroFix:
        for i in range(1, result.GetNbinsX() + 1):
            if result.GetBinContent(i) == 0 and den_fix.GetBinContent(i) == 0:
                result.SetBinContent(i, 1)
    return result


def MakePurityHist(sighist, stack, category):
    """Make a new purity TH1 from stack(bkgs TH1s) and sighist TH1

    Args:
        hist (TH1): histogram

    Returns:
        TH1: A new TH1 containing the bin purity
    """
    total = sighist.Clone()
    bkghists = stack.GetHists()
    hist_dict = {}
    for bkghist in bkghists:
        # define a dict to use later
        hist_dict[bkghist.GetName()] = bkghist
        total.Add(bkghist.Clone())
    fraction = 0.0
    fraction_hist = sighist.Clone()
    for i in range(0, fraction_hist.GetNbinsX() + 1):
        fraction_hist.SetBinContent(i, 0.0)
        fraction_hist.SetBinError(i, 0.0)
    for i in range(1, total.GetNbinsX() + 1):
        total_content = total.GetBinContent(i)

        if category == "NN_higgs":
            fraction = sighist.GetBinContent(i) / total_content
        elif category == "NN_zttEmbed":
            fraction = hist_dict["EmbedZTT"].GetBinContent(i) / total_content
        elif category == "NN_jetFakes":
            fraction = hist_dict["jetFakes"].GetBinContent(i) / total_content
        fraction_hist.SetBinContent(i, fraction)
        fraction_hist.SetBinError(i, np.sqrt(total_content))

    return fraction_hist


def RemoveGraphXDuplicates(graph):
    for i in range(graph.GetN() - 1):
        if graph.GetX()[i + 1] == graph.GetX()[i]:
            # print 'Removing duplicate point (%f, %f)' % (graph.GetX()[i+1], graph.GetY()[i+1])
            graph.RemovePoint(i + 1)
            RemoveGraphXDuplicates(graph)
            break


def ApplyGraphYOffset(graph, y_off):
    for i in range(graph.GetN() - 1):
        graph.GetY()[i] = graph.GetY()[i] + y_off


def RemoveGraphYAll(graph, val):
    for i in range(graph.GetN()):
        if graph.GetY()[i] == val:
            print(
                "[RemoveGraphYAll] Removing point (%f, %f)"
                % (graph.GetX()[i], graph.GetY()[i])
            )
            graph.RemovePoint(i)
            RemoveGraphYAll(graph, val)
            break


def RemoveSmallDelta(graph, val):
    for i in range(graph.GetN()):
        diff = abs(graph.GetY()[i])
        if diff < val:
            print(
                "[RemoveSmallDelta] Removing point (%f, %f)"
                % (graph.GetX()[i], graph.GetY()[i])
            )
            graph.RemovePoint(i)
            RemoveSmallDelta(graph, val)
            break


def RemoveGraphYAbove(graph, val):
    for i in range(graph.GetN()):
        if graph.GetY()[i] > val:
            # print 'Removing point (%f, %f)' % (graph.GetX()[i],
            # graph.GetY()[i])
            graph.RemovePoint(i)
            RemoveGraphYAbove(graph, val)
            break


def ImproveMinimum(graph, func, doIt=False):
    fit_x = 0.0
    fit_y = 999.0
    fit_i = 0
    for i in range(graph.GetN()):
        if graph.GetY()[i] < fit_y:
            fit_i = i
            fit_x = graph.GetX()[i]
            fit_y = graph.GetY()[i]
    if fit_i == 0 or fit_i == (graph.GetN() - 1):
        if doIt:
            min_x = graph.GetX()[fit_i]
            min_y = graph.GetY()[fit_i]
            for i in range(graph.GetN()):
                before = graph.GetY()[i]
                graph.GetY()[i] -= min_y
                after = graph.GetY()[i]
                print("Point %i, before=%f, after=%f" % (i, before, after))
        return (fit_x, fit_y)
    search_min = fit_i - 2 if fit_i >= 2 else fit_i - 1
    search_max = fit_i + 2 if fit_i + 2 < graph.GetN() else fit_i + 1
    min_x = func.GetMinimumX(graph.GetX()[search_min], graph.GetX()[search_max])
    min_y = func.Eval(min_x)
    print("[ImproveMinimum] Fit minimum was (%f, %f)" % (fit_x, fit_y))
    print("[ImproveMinimum] Better minimum was (%f, %f)" % (min_x, min_y))
    if doIt:
        for i in range(graph.GetN()):
            before = graph.GetY()[i]
            graph.GetY()[i] -= min_y
            after = graph.GetY()[i]
            print("Point %i, before=%f, after=%f" % (i, before, after))
        graph.Set(graph.GetN() + 1)
        graph.SetPoint(graph.GetN() - 1, min_x, 0)
        graph.Sort()
    return (min_x, min_y)


def FindCrossingsWithSpline(graph, func, yval):
    crossings = []
    intervals = []
    current = None
    for i in range(graph.GetN() - 1):
        if (graph.GetY()[i] - yval) * (graph.GetY()[i + 1] - yval) < 0.0:
            cross = func.GetX(yval, graph.GetX()[i], graph.GetX()[i + 1])
            if (graph.GetY()[i] - yval) > 0.0 and current is None:
                current = {
                    "lo": cross,
                    "hi": graph.GetX()[graph.GetN() - 1],
                    "valid_lo": True,
                    "valid_hi": False,
                }
            if (graph.GetY()[i] - yval) < 0.0 and current is None:
                current = {
                    "lo": graph.GetX()[0],
                    "hi": cross,
                    "valid_lo": False,
                    "valid_hi": True,
                }
                intervals.append(current)
                current = None
            if (graph.GetY()[i] - yval) < 0.0 and current is not None:
                current["hi"] = cross
                current["valid_hi"] = True
                intervals.append(current)
                current = None
            # print 'Crossing between: (%f, %f) -> (%f, %f) at %f' %
            # (graph.GetX()[i], graph.GetY()[i], graph.GetX()[i+1],
            # graph.GetY()[i+1], cross)
            crossings.append(cross)
    if current is not None:
        intervals.append(current)
    if len(intervals) == 0:
        current = {
            "lo": graph.GetX()[0],
            "hi": graph.GetX()[graph.GetN() - 1],
            "valid_lo": False,
            "valid_hi": False,
        }
        intervals.append(current)
    print(intervals)
    return intervals
    # return crossings


def ReZeroTGraph(gr, doIt=False):
    fit_x = 0.0
    fit_y = 0.0
    for i in range(gr.GetN()):
        if gr.GetY()[i] == 0.0:
            fit_x = gr.GetX()[i]
            fit_y = gr.GetY()[i]
            break
    min_x = 0.0
    min_y = 0.0
    for i in range(gr.GetN()):
        if gr.GetY()[i] < min_y:
            min_y = gr.GetY()[i]
            min_x = gr.GetX()[i]
    if min_y < fit_y:
        print("[ReZeroTGraph] Fit minimum was (%f, %f)" % (fit_x, fit_y))
        print("[ReZeroTGraph] Better minimum was (%f, %f)" % (min_x, min_y))
        if doIt:
            for i in range(gr.GetN()):
                before = gr.GetY()[i]
                gr.GetY()[i] -= min_y
                after = gr.GetY()[i]
                print("Point %i, before=%f, after=%f" % (i, before, after))
    return min_y


def RemoveNearMin(graph, val, spacing=None):
    # assume graph is sorted:
    n = graph.GetN()
    if n < 5:
        return
    if spacing is None:
        spacing = (graph.GetX()[n - 1] - graph.GetX()[0]) / float(n - 2)
        # print '[RemoveNearMin] Graph has spacing of %.3f' % spacing
    bf_i = None
    for i in range(graph.GetN()):
        if graph.GetY()[i] == 0.0:
            bf = graph.GetX()[i]
            bf_i = i
            # print '[RemoveNearMin] Found best-fit at %.3f' % bf
            break
    if bf_i is None:
        print("[RemoveNearMin] No minimum found!")
        return
    for i in range(graph.GetN()):
        if i == bf_i:
            continue
        if abs(graph.GetX()[i] - bf) < (val * spacing):
            print(
                "[RemoveNearMin] Removing point (%f, %f) close to minimum at %f"
                % (graph.GetX()[i], graph.GetY()[i], bf)
            )
            graph.RemovePoint(i)
            RemoveNearMin(graph, val, spacing)
            break


def SortGraph(Graph):
    sortedGraph = R.TGraph()
    graph_list = []
    for i in range(Graph.GetN()):
        graph_list.append((float(Graph.GetX()[i]), float(Graph.GetY()[i])))
    graph_list = sorted(set(graph_list))
    for i in range(Graph.GetN()):
        sortedGraph.SetPoint(i, graph_list[i][0], graph_list[i][1])
    return sortedGraph


def FixTopRange(pad, fix_y, fraction):
    hobj = GetAxisHist(pad)
    ymin = hobj.GetMinimum()
    hobj.SetMaximum((fix_y - fraction * ymin) / (1.0 - fraction))
    if R.gPad.GetLogy():
        if ymin == 0.0:
            print("Cannot adjust log-scale y-axis range if the minimum is zero!")
            return
        maxval = (math.log10(fix_y) - fraction * math.log10(ymin)) / (1 - fraction)
        maxval = math.pow(10, maxval)
        hobj.SetMaximum(maxval)


def FixBothRanges(pad, fix_y_lo, frac_lo, fix_y_hi, frac_hi):
    """Adjusts y-axis range such that a lower and a higher value are located a
    fixed fraction of the frame height away from a new minimum and maximum
    respectively.

    This function is useful in conjunction with GetPadYMax which returns the
    maximum or minimum y value of all histograms and graphs drawn on the pad.

    In the example below, the minimum and maximum values found via this function
    are used as the `fix_y_lo` and `fix_y_hi` arguments, and the spacing fractions
    as 0.15 and 0.30 respectively.

    @code
    FixBothRanges(pad, GetPadYMin(pad), 0.15, GetPadYMax(pad), 0.30)
    @endcode

    ![](figures/FixBothRanges.png)

    Args:
        pad (TPad): A TPad on which histograms and graphs have already been drawn
        fix_y_lo (float): The y value which will end up a fraction `frac_lo` above
                          the new axis minimum.
        frac_lo (float): A fraction of the y-axis height
        fix_y_hi (float): The y value which will end up a fraction `frac_hi` below
                         from the new axis maximum.
        frac_hi (float): A fraction of the y-axis height
    """
    hobj = GetAxisHist(pad)
    ymin = fix_y_lo
    ymax = fix_y_hi
    if R.gPad.GetLogy():
        if ymin == 0.0:
            print("Cannot adjust log-scale y-axis range if the minimum is zero!")
            return
        ymin = math.log10(ymin)
        ymax = math.log10(ymax)
    fl = frac_lo
    fh = frac_hi

    ymaxn = (
        (1.0 / (1.0 - (fh * fl / ((1.0 - fl) * (1.0 - fh)))))
        * (1.0 / (1.0 - fh))
        * (ymax - fh * ymin)
    )
    yminn = (ymin - fl * ymaxn) / (1.0 - fl)
    if R.gPad.GetLogy():
        yminn = math.pow(10, yminn)
        ymaxn = math.pow(10, ymaxn)
    hobj.SetMinimum(yminn)
    hobj.SetMaximum(ymaxn)


def GetPadYMaxInRange(pad, x_min, x_max, do_min=False):
    pad_obs = pad.GetListOfPrimitives()
    if pad_obs is None:
        return 0.0
    h_max = -99999.0
    h_min = +99999.0
    for obj in pad_obs:
        if obj.InheritsFrom(R.TH1.Class()):
            hobj = obj
            for j in range(1, hobj.GetNbinsX() + 1):
                if (
                    hobj.GetBinLowEdge(j) + hobj.GetBinWidth(j) < x_min
                    or hobj.GetBinLowEdge(j) > x_max
                ):
                    continue
                if hobj.GetBinContent(j) + hobj.GetBinError(j) > h_max:
                    h_max = hobj.GetBinContent(j) + hobj.GetBinError(j)
                if (hobj.GetBinContent(j) - hobj.GetBinError(j) < h_min) and not do_min:
                    # If we're looking for the minimum don't count TH1s
                    # because we probably only care about graphs
                    h_min = hobj.GetBinContent(j) - hobj.GetBinError(j)
        elif obj.InheritsFrom(R.TGraphAsymmErrors.Class()):
            gobj = obj
            n = gobj.GetN()
            for k in range(0, n):
                x = gobj.GetX()[k]
                y = gobj.GetY()[k]
                if x < x_min or x > x_max:
                    continue
                if (y + gobj.GetEYhigh()[k]) > h_max:
                    h_max = y + gobj.GetEYhigh()[k]
                if (y - gobj.GetEYlow()[k]) < h_min:
                    h_min = y - gobj.GetEYlow()[k]
        elif obj.InheritsFrom(R.TGraphErrors.Class()):
            gobj = obj
            n = gobj.GetN()
            for k in range(0, n):
                x = gobj.GetX()[k]
                y = gobj.GetY()[k]
                if x < x_min or x > x_max:
                    continue
                if (y + gobj.GetEY()[k]) > h_max:
                    h_max = y + gobj.GetEY()[k]
                if (y - gobj.GetEY()[k]) < h_min:
                    h_min = y - gobj.GetEY()[k]
        elif obj.InheritsFrom(R.TGraph.Class()):
            gobj = obj
            n = gobj.GetN()
            for k in range(0, n):
                x = gobj.GetX()[k]
                y = gobj.GetY()[k]
                if x < x_min or x > x_max:
                    continue
                if y > h_max:
                    h_max = y
                if y < h_min:
                    h_min = y
    return h_max if do_min is False else h_min


def GetPadYMax(pad, do_min=False):
    pad_obs = pad.GetListOfPrimitives()
    if pad_obs is None:
        return 0.0
    xmin = GetAxisHist(pad).GetXaxis().GetXmin()
    xmax = GetAxisHist(pad).GetXaxis().GetXmax()
    return GetPadYMaxInRange(pad, xmin, xmax, do_min)


def GetPadYMin(pad):
    return GetPadYMax(pad, True)


def FixOverlay():
    R.gPad.GetFrame().Draw()
    R.gPad.RedrawAxis()


def StandardAxes(xaxis, yaxis, var, units):
    width = xaxis.GetBinWidth(1)
    w_label = "%.1f" % width
    if units == "":
        xaxis.SetTitle(var)
        yaxis.SetTitle("Events / " + w_label)
    else:
        xaxis.SetTitle(var + " (" + units + ")")
        yaxis.SetTitle("Events / " + w_label + " " + units)


def DrawCMSLogo(
    pad,
    cmsText,
    extraText,
    iPosX,
    relPosX,
    relPosY,
    relExtraDY,
    extraText2="",
    cmsTextSize=0.8,
):
    """Blah

    Args:
        pad (TYPE): Description
        cmsText (TYPE): Description
        extraText (TYPE): Description
        iPosX (TYPE): Description
        relPosX (TYPE): Description
        relPosY (TYPE): Description
        relExtraDY (TYPE): Description
        extraText2 (str): Description
        cmsTextSize (float): Description

    Returns:
        TYPE: Description
    """
    pad.cd()
    cmsTextFont = 62  # default is helvetic-bold

    writeExtraText = len(extraText) > 0
    writeExtraText2 = len(extraText2) > 0
    extraTextFont = 52

    # text sizes and text offsets with respect to the top frame
    # in unit of the top margin size
    lumiTextOffset = 0.2
    # cmsTextSize = 0.8
    # float cmsTextOffset    = 0.1;  // only used in outOfFrame version

    # ratio of 'CMS' and extra text size
    extraOverCmsTextSize = 0.76

    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True

    alignY_ = 3
    alignX_ = 2
    if iPosX / 10 == 0:
        alignX_ = 1
    if iPosX == 0:
        alignX_ = 1
    if iPosX == 0:
        alignY_ = 1
    if iPosX / 10 == 1:
        alignX_ = 1
    if iPosX / 10 == 2:
        alignX_ = 2
    if iPosX / 10 == 3:
        alignX_ = 3
    # if (iPosX == 0): relPosX = 0.14
    align_ = 10 * alignX_ + alignY_

    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()

    latex = R.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(R.kBlack)

    extraTextSize = extraOverCmsTextSize * cmsTextSize
    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / (
        float(pad.GetWw()) * pad.GetAbsWNDC()
    )
    if pad_ratio < 1.0:
        pad_ratio = 1.0

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)

    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif iPosX % 10 == 2:
        posX_ = l + 0.5 * (1 - l - r)
    elif iPosX % 10 == 3:
        posX_ = 1 - r - relPosX * (1 - l - r)

    posY_ = 1 - t - relPosY * (1 - t - b)
    if not outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_ - 0.04, posY_, cmsText)
        if writeExtraText:
            latex.SetTextFont(extraTextFont)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extraTextSize * t * pad_ratio)
            latex.DrawLatex(posX_, posY_ - relExtraDY * cmsTextSize * t, extraText)
            if writeExtraText2:
                latex.DrawLatex(
                    posX_, posY_ - 1.8 * relExtraDY * cmsTextSize * t, extraText2
                )
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r)
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)


def PositionedLegend(width, height, pos, offset):
    o = offset
    w = width
    h = height
    l = R.gPad.GetLeftMargin()
    t = R.gPad.GetTopMargin()
    b = R.gPad.GetBottomMargin()
    r = R.gPad.GetRightMargin()
    if pos == 1:
        return R.TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, "", "NBNDC")
    if pos == 2:
        c = l + 0.5 * (1 - l - r)
        return R.TLegend(
            c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, "", "NBNDC"
        )
    if pos == 3:
        return R.TLegend(
            1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, "", "NBNDC"
        )
    if pos == 4:
        return R.TLegend(l + o, b + o, l + o + w, b + o + h, "", "NBNDC")
    if pos == 5:
        c = l + 0.5 * (1 - l - r)
        return R.TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, "", "NBNDC")
    if pos == 6:
        return R.TLegend(1 - r - o - w, b + o, 1 - r - o, b + o + h, "", "NBNDC")
    if pos == 7:
        return R.TLegend(1 - o - w, 1 - t - o - h, 1 - o, 1 - t - o, "", "NBNDC")


def DrawHorizontalLine(pad, line, yval):
    axis = GetAxisHist(pad)
    xmin = axis.GetXaxis().GetXmin()
    xmax = axis.GetXaxis().GetXmax()
    line.DrawLine(xmin, yval, xmax, yval)


def DrawVerticalLine(pad, line, xval):
    axis = GetAxisHist(pad)
    ymin = axis.GetMinimum()
    ymax = axis.GetMaximum()
    line.DrawLine(xval, ymin, xval, ymax)


def DrawTitle(pad, text, align, scale=1):
    pad_backup = R.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / (
        float(pad.GetWw()) * pad.GetAbsWNDC()
    )
    if pad_ratio < 1.0:
        pad_ratio = 1.0

    textSize = 0.6
    textOffset = 0.2

    latex = R.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(R.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t
    if align == 1:
        latex.SetTextAlign(11)
    if align == 1:
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
    if align == 2:
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
    if align == 3:
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()


def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a - b) <= max(abs_tol, rel_tol * max(abs(a), abs(b)))


def StyleLimitBand(graph_dict, overwrite_style_dict=None):
    style_dict = {
        "obs": {"LineWidth": 2},
        "exp0": {"LineWidth": 2, "LineColor": R.kRed},
        "exp1": {"FillColor": R.kGreen},
        "exp2": {"FillColor": R.kYellow},
    }
    if overwrite_style_dict is not None:
        for key in overwrite_style_dict:
            if key in style_dict:
                style_dict[key].update(overwrite_style_dict[key])
            else:
                style_dict[key] = overwrite_style_dict[key]
    for key in graph_dict:
        Set(graph_dict[key], **style_dict[key])


def DrawLimitBand(
    pad,
    graph_dict,
    draw=["exp2", "exp1", "exp0", "obs"],
    draw_legend=None,
    legend=None,
    legend_overwrite=None,
):
    legend_dict = {
        "obs": {"Label": "Observed", "LegendStyle": "LP", "DrawStyle": "PLSAME"},
        "exp0": {"Label": "Expected", "LegendStyle": "L", "DrawStyle": "LSAME"},
        "exp1": {
            "Label": "#pm1#sigma Expected",
            "LegendStyle": "F",
            "DrawStyle": "3SAME",
        },
        "exp2": {
            "Label": "#pm2#sigma Expected",
            "LegendStyle": "F",
            "DrawStyle": "3SAME",
        },
    }
    if legend_overwrite is not None:
        for key in legend_overwrite:
            if key in legend_dict:
                legend_dict[key].update(legend_overwrite[key])
            else:
                legend_dict[key] = legend_overwrite[key]
    pad.cd()
    for key in draw:
        if key in graph_dict:
            graph_dict[key].Draw(legend_dict[key]["DrawStyle"])
    if legend is not None:
        if draw_legend is None:
            draw_legend = reversed(draw)
        for key in draw_legend:
            if key in graph_dict:
                legend.AddEntry(
                    graph_dict[key],
                    legend_dict[key]["Label"],
                    legend_dict[key]["LegendStyle"],
                )


def contourFromTH2(h2in, threshold, minPoints=10, frameValue=1000.0):
    # // http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
    contoursList = [threshold]
    contours = array("d", contoursList)
    # if (h2in.GetNbinsX() * h2in.GetNbinsY()) > 10000: minPoints = 50
    # if (h2in.GetNbinsX() * h2in.GetNbinsY()) <= 100: minPoints = 10

    h2 = frameTH2D(h2in, threshold, frameValue)

    h2.SetContour(1, contours)

    # Draw contours as filled regions, and Save points
    # backup = R.gPad # doesn't work in pyroot, backup behaves like a ref to gPad
    canv = R.TCanvas("tmp", "tmp")
    canv.cd()
    h2.Draw("CONT Z LIST")
    R.gPad.Update()  # Needed to force the plotting and retrieve the contours in

    conts = R.gROOT.GetListOfSpecials().FindObject("contours")
    contLevel = None

    if conts is None or conts.GetSize() == 0:
        print("*** No Contours Were Extracted!")
        return None
    ret = R.TList()
    for i in range(conts.GetSize()):
        contLevel = conts.At(i)
        print(">> Contour %d has %d Graphs" % (i, contLevel.GetSize()))
        for j in range(contLevel.GetSize()):
            gr1 = contLevel.At(j)
            print("\t Graph %d has %d points" % (j, gr1.GetN()))
            if gr1.GetN() > minPoints:
                ret.Add(gr1.Clone())
            # // break;
    # backup.cd()
    canv.Close()
    return ret


def frameTH2D(hist, threshold, frameValue=1000):
    # Now supports variable-binned histograms First adds a narrow frame (1% of
    # of bin widths) around the outside with same values as the real edge. Then
    # adds another frame another frame around this one filled with some chosen
    # value that will make the contours close

    # Get lists of the bin edges
    x_bins = [hist.GetXaxis().GetBinLowEdge(x) for x in range(1, hist.GetNbinsX() + 2)]
    y_bins = [hist.GetYaxis().GetBinLowEdge(y) for y in range(1, hist.GetNbinsY() + 2)]

    # New bin edge arrays will need an extra four values
    x_new = [0.0] * (len(x_bins) + 4)
    y_new = [0.0] * (len(y_bins) + 4)

    # Calculate bin widths at the edges
    xw1 = x_bins[1] - x_bins[0]
    xw2 = x_bins[-1] - x_bins[-2]
    yw1 = y_bins[1] - y_bins[0]
    yw2 = y_bins[-1] - y_bins[-2]

    # Set the edges of the outer framing bins and the adjusted
    # edge of the real edge bins
    x_new[0] = x_bins[0] - 2 * xw1 * 0.02
    x_new[1] = x_bins[0] - 1 * xw1 * 0.02
    x_new[-1] = x_bins[-1] + 2 * xw2 * 0.02
    x_new[-2] = x_bins[-1] + 1 * xw2 * 0.02
    y_new[0] = y_bins[0] - 2 * yw1 * 0.02
    y_new[1] = y_bins[0] - 1 * yw1 * 0.02
    y_new[-1] = y_bins[-1] + 2 * yw2 * 0.02
    y_new[-2] = y_bins[-1] + 1 * yw2 * 0.02

    # Copy the remaining bin edges from the hist
    for i in range(0, len(x_bins)):
        x_new[i + 2] = x_bins[i]
    for i in range(0, len(y_bins)):
        y_new[i + 2] = y_bins[i]

    # print x_new
    # print y_new

    framed = R.TH2D(
        "%s framed" % hist.GetName(),
        "%s framed" % hist.GetTitle(),
        len(x_new) - 1,
        array("d", x_new),
        len(y_new) - 1,
        array("d", y_new),
    )
    framed.SetDirectory(0)

    for x in range(1, framed.GetNbinsX() + 1):
        for y in range(1, framed.GetNbinsY() + 1):
            if x == 1 or x == framed.GetNbinsX() or y == 1 or y == framed.GetNbinsY():
                # This is a a frame bin
                framed.SetBinContent(x, y, frameValue)
            else:
                # adjust x and y if we're in the first frame so as to copy the output
                # values from the real TH2
                ux = x
                uy = y
                if x == 2:
                    ux += 1
                elif x == (len(x_new) - 2):
                    ux -= 1
                if y == 2:
                    uy += 1
                elif y == (len(y_new) - 2):
                    uy -= 1
                framed.SetBinContent(x, y, hist.GetBinContent(ux - 2, uy - 2))
    return framed


def fastFillTH2(hist2d, graph, initalValue=99999, interpolateMissing=False):
    for x in range(1, hist2d.GetNbinsX() + 1):
        for y in range(1, hist2d.GetNbinsY() + 1):
            hist2d.SetBinContent(x, y, initalValue)
    # for i in range(graph.GetN()):
    # hist2d.Fill(graph.GetX()[i],graph.GetY()[i],graph.GetZ()[i])
    for i in range(graph.GetN()):
        xbin = hist2d.GetXaxis().FindBin(graph.GetX()[i])
        ybin = hist2d.GetYaxis().FindBin(graph.GetY()[i])
        if isclose(
            hist2d.GetXaxis().GetBinCenter(xbin), graph.GetX()[i], rel_tol=1e-2
        ) and isclose(
            hist2d.GetYaxis().GetBinCenter(ybin), graph.GetY()[i], rel_tol=1e-2
        ):
            hist2d.SetBinContent(xbin, ybin, graph.GetZ()[i])
    interpolated = 0
    if interpolateMissing:
        for x in range(1, hist2d.GetNbinsX() + 1):
            for y in range(1, hist2d.GetNbinsY() + 1):
                if hist2d.GetBinContent(x, y) == initalValue:
                    interpolated += 1
                    hist2d.SetBinContent(
                        x,
                        y,
                        graph.Interpolate(
                            hist2d.GetXaxis().GetBinCenter(x),
                            hist2d.GetYaxis().GetBinCenter(y),
                        ),
                    )


def fillTH2(hist2d, graph):
    for x in range(1, hist2d.GetNbinsX() + 1):
        for y in range(1, hist2d.GetNbinsY() + 1):
            xc = hist2d.GetXaxis().GetBinCenter(x)
            yc = hist2d.GetYaxis().GetBinCenter(y)
            val = graph.Interpolate(xc, yc)
            hist2d.SetBinContent(x, y, val)


# Functions 'NewInterpolate' and 'rebin' are taken, translated and modified into python from:
# https://indico.cern.ch/event/256523/contribution/2/attachments/450198/624259/07JUN2013_cawest.pdf
# http://hep.ucsb.edu/people/cawest/interpolation/interpolate.h
def NewInterpolate(hist):
    histCopy = hist.Clone()

    # make temporary histograms to store the results of both steps
    hist_step1 = histCopy.Clone()
    hist_step1.Reset()
    hist_step2 = histCopy.Clone()
    hist_step2.Reset()

    nBinsX = histCopy.GetNbinsX()
    nBinsY = histCopy.GetNbinsY()

    xMin = 1
    yMin = 1
    xMax = histCopy.GetNbinsX() + 1
    yMax = histCopy.GetNbinsY() + 1

    for i in range(1, nBinsX + 1):
        for j in range(1, nBinsY + 1):
            # do not extrapolate outside the scan
            if (i < xMin) or (i > xMax) or (j < yMin) or (j > yMax):
                continue
            binContent = histCopy.GetBinContent(i, j)
            binContentNW = histCopy.GetBinContent(i + 1, j + 1)
            binContentSE = histCopy.GetBinContent(i - 1, j - 1)
            binContentNE = histCopy.GetBinContent(i + 1, j - 1)
            binContentSW = histCopy.GetBinContent(i - 1, j + 1)
            binContentUp = histCopy.GetBinContent(i, j + 1)
            binContentDown = histCopy.GetBinContent(i, j - 1)
            binContentLeft = histCopy.GetBinContent(i - 1, j)
            binContentRight = histCopy.GetBinContent(i + 1, j)
            nFilled = 0
            if binContentNW > 0:
                nFilled += 1
            if binContentSE > 0:
                nFilled += 1
            if binContentNE > 0:
                nFilled += 1
            if binContentSW > 0:
                nFilled += 1
            if binContentUp > 0:
                nFilled += 1
            if binContentDown > 0:
                nFilled += 1
            if binContentRight > 0:
                nFilled += 1
            if binContentLeft > 0:
                nFilled += 1
            # if we are at an empty bin and there are neighbors
            # in specified direction with non-zero entries
            if (binContent == 0) and (nFilled > 1):
                # average over non-zero entries
                binContent = (
                    binContentNW
                    + binContentSE
                    + binContentNE
                    + binContentSW
                    + binContentUp
                    + binContentDown
                    + binContentRight
                    + binContentLeft
                ) / nFilled
                hist_step1.SetBinContent(i, j, binContent)

        # add result of interpolation
    histCopy.Add(hist_step1)

    for i in range(1, nBinsX):
        for j in range(1, nBinsY):
            if (i < xMin) or (i > xMax) or (j < yMin) or (j > yMax):
                continue
            binContent = histCopy.GetBinContent(i, j)
            # get entries for "Swiss Cross" average
            binContentUp = histCopy.GetBinContent(i, j + 1)
            binContentDown = histCopy.GetBinContent(i, j - 1)
            binContentLeft = histCopy.GetBinContent(i - 1, j)
            binContentRight = histCopy.GetBinContent(i + 1, j)
            nFilled = 0
            if binContentUp > 0:
                nFilled += 1
            if binContentDown > 0:
                nFilled += 1
            if binContentRight > 0:
                nFilled += 1
            if binContentLeft > 0:
                nFilled += 1
            if (binContent == 0) and (nFilled > 0):
                # only average over non-zero entries
                binContent = (
                    binContentUp + binContentDown + binContentRight + binContentLeft
                ) / nFilled
                hist_step2.SetBinContent(i, j, binContent)
    # add "Swiss Cross" average
    histCopy.Add(hist_step2)

    return histCopy


def rebin(hist):
    histName = hist.GetName()
    histName += "_rebin"

    # bin widths are needed so as to not shift histogram by half a bin with each rebinning
    # assume constant binning
    #  binWidthX = hist.GetXaxis().GetBinWidth(1)
    #  binWidthY = hist.GetYaxis().GetBinWidth(1)

    #  histRebinned = R.TH2F(histName, histName, 2*hist.GetNbinsX(), hist.GetXaxis().GetXmin()+binWidthX/4, hist.GetXaxis().GetXmax()+binWidthX/4, 2*hist.GetNbinsY(), hist.GetYaxis().GetXmin()+binWidthY/4, hist.GetYaxis().GetXmax()+binWidthY/4)
    histRebinned = R.TH2F(
        histName,
        histName,
        2 * hist.GetNbinsX() - 1,
        hist.GetXaxis().GetXmin(),
        hist.GetXaxis().GetXmax(),
        2 * hist.GetNbinsY() - 1,
        hist.GetYaxis().GetXmin(),
        hist.GetYaxis().GetXmax(),
    )

    # copy results from previous histogram
    for iX in range(1, hist.GetNbinsX() + 1):
        for iY in range(1, hist.GetNbinsY() + 1):
            binContent = hist.GetBinContent(iX, iY)
            histRebinned.SetBinContent(2 * iX - 1, 2 * iY - 1, binContent)
    histRebinned.SetMaximum(hist.GetMaximum())
    histRebinned.SetMinimum(hist.GetMinimum())

    # use interpolation to re-fill histogram
    histRebinnedInterpolated = NewInterpolate(histRebinned)

    return histRebinnedInterpolated


def higgsConstraint(model, higgstype):
    higgsBand = R.TGraph2D()
    masslow = 150
    masshigh = 500
    massstep = 10
    n = 0
    for mass in range(masslow, masshigh, massstep):
        myfile = open(
            "../../HiggsAnalysis/HiggsToTauTau/data/Higgs125/"
            + model
            + "/higgs_"
            + str(mass)
            + ".dat",
            "r",
        )
        for line in myfile:
            tanb = (line.split())[0]
            mh = float((line.split())[1])
            mH = float((line.split())[3])
            if higgstype == "h":
                higgsBand.SetPoint(n, mass, float(tanb), mh)
            elif higgstype == "H":
                higgsBand.SetPoint(n, mass, float(tanb), mH)
            n = n + 1
        myfile.close()
    return higgsBand


def getOverlayMarkerAndLegend(
    legend, entries, options, borderSize=2.0 / 3, markerStyle="P"
):
    borderLegend = legend.Clone()
    borderLegend.Clear()
    graphs = []
    for i in range(legend.GetNRows()):
        if i in entries:
            graph = entries[i].Clone()
            options[i]["MarkerSize"] = graph.GetMarkerSize() * borderSize
            Set(graph, **options[i])
            borderLegend.AddEntry(graph, " ", markerStyle)
            graphs.append(graph)
        else:
            borderLegend.AddEntry("", " ", "")
    borderLegend.SetFillStyle(0)
    borderLegend.SetFillColor(0)
    return (borderLegend, graphs)


def signalComp(leg, plots, colour, stacked):
    return dict(
        [
            ("leg_text", leg),
            ("plot_list", plots),
            ("colour", colour),
            ("in_stack", stacked),
        ]
    )


def backgroundComp(leg, plots, colour):
    return dict([("leg_text", leg), ("plot_list", plots), ("colour", colour)])


def createAxisHists(n, src, xmin=0, xmax=499):
    result = []
    for i in range(0, n):
        res = src.Clone()
        res.Reset()
        res.SetTitle("")
        res.SetName("axis%(i)d" % vars())
        res.SetAxisRange(xmin, xmax)
        res.SetStats(0)
        result.append(res)
    return result


def PassAutoBlindMetric(s, b, epsilon=0.09, metric=0.5):
    if b <= 0:
        return True
    y = s / math.sqrt(b + (epsilon * b) ** 2)
    return y >= metric


def Norm2DBins(h):
    for i in range(1, h.GetNbinsX() + 1):
        bin_lab = h.GetXaxis().GetBinLabel(i)
        scale = float(bin_lab.split("-")[1]) - float(bin_lab.split("-")[0])
        content = h.GetBinContent(i)
        error = h.GetBinError(i)
        h.SetBinContent(i, content / scale)
        h.SetBinError(i, error / scale)


def ConvertHistToDiscreteBins(hist, bin_labels):
    # get current binning
    curr_binning = []
    for r in range(1, hist.GetNbinsX() + 2):
        curr_binning.append(hist.GetBinLowEdge(r))

    do_geq = False
    do_eq = False
    for i in bin_labels:
        if ">=" in i:
            do_geq = True
        if "==" in i:
            do_eq = True

    # if do_geq add extra bin so >= bin included in plot
    if do_geq:
        new_binning = []
        for r in curr_binning:
            if r == curr_binning[-1]:
                new_binning.append(r)
                new_binning.append(2 * curr_binning[-1] - curr_binning[-2])
            else:
                new_binning.append(r)
        new_bins = array("f", map(float, new_binning))
        h_clone = hist.Clone()
        hist = R.TH1F(h_clone.GetName(), "", len(new_bins) - 1, new_bins)
        for r in range(1, hist.GetNbinsX() + 2):
            hist.SetBinContent(r, h_clone.GetBinContent(r))
            hist.SetBinError(r, h_clone.GetBinError(r))

        bin_names = []
        for i in bin_labels:
            bin_names.append(str(i).replace(">=", "#geq"))

    if do_eq:
        h_clone = hist.Clone()
        new_bins = range(1, len(bin_labels) + 1)
        new_bins_array = array("f", range(1, len(bin_labels) + 2))
        hist = R.TH1F(h_clone.GetName(), "", len(new_bins_array) - 1, new_bins_array)
        for i in new_bins:
            hist.SetBinContent(
                i,
                h_clone.GetBinContent(
                    h_clone.FindBin(int(bin_labels[i - 1].replace("==", "")))
                ),
            )
            hist.SetBinError(
                i,
                h_clone.GetBinError(
                    h_clone.FindBin(int(bin_labels[i - 1].replace("==", "")))
                ),
            )

        bin_names = []
        for i in bin_labels:
            bin_names.append(str(i).replace("==", ""))

    return hist, bin_names


def NonZeroMinimum(h):
    min_value = 9999
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) > 0 and h.GetBinContent(i) < min_value:
            min_value = h.GetBinContent(i)
    return min_value


def HTTPlot(
    nodename,
    infile=None,
    signal_scale=1,
    signal_mass="",
    FF=False,
    norm_bins=True,
    channel="mt",
    blind=False,
    x_blind_min=-1e5,
    x_blind_max=1e5,
    ratio=True,
    threePads=False,
    ratio_log_y=False,
    log_y=False,
    log_x=False,
    ratio_range="0.7,1.3",
    custom_x_range=False,
    x_axis_max=4000,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=4000,
    y_axis_min=0,
    x_title="",
    y_title="",
    extra_pad=0,
    signal_scheme="run2_mssm",
    do_custom_uncerts=False,
    add_stat_to_syst=False,
    add_flat_uncert=False,
    uncert_title="background uncertainty",
    lumi="35.9",
    plot_name="htt_plot",
    custom_uncerts_up_name="total_bkg_custom_uncerts_up",
    custom_uncerts_down_name="total_bkg_custom_uncerts_down",
    scheme="mt",
    embedding=False,
    vbf_background=False,
    split_sm_scheme=False,
    ggh_scheme="powheg",
    cat="",
    split_taus=False,
    auto_blind=True,
    discrete_x_axis=False,
    discrete_x_labels=None,
    qcd_ff_closure=False,
    w_ff_closure=False,
    bkg_comp=False,
    plot_signals=[],
    qcdname="QCD",  # Only used with tt for specific analysis, not implemented fully
):

    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    # Define signal schemes here
    sig_schemes = {}
    sig_schemes["sm_ggH"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggH_ph_htt"],
        False,
    )
    sig_schemes["sm_ggH_JHU"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggHsm_htt"],
        False,
    )
    sig_schemes["sm_qqH"] = (
        str(int(signal_scale)) + "#times SM qqH#rightarrow#tau#tau",
        ["qqH_htt"],
        False,
    )
    sig_schemes["sm_VH"] = (
        str(int(signal_scale)) + "#times SM VH#rightarrow#tau#tau",
        ["WminusH_htt", "WplusH_htt", "ZH_htt"],
        False,
    )
    sig_schemes["sm_default"] = (
        str(int(signal_scale))
        + "#times SM H("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH_ph_htt", "qqH_htt"],
        True,
    )
    sig_schemes["smsummer16"] = (
        str(int(signal_scale))
        + "#times SM H("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH_htt", "qqH_htt", "WminusH_htt", "WplusH_htt", "ZH_htt"],
        False,
    )
    sig_schemes["smsummer17"] = (
        str(int(signal_scale))
        + "#times SM H("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH_htt", "qqH_htt"],
        False,
    )
    sig_schemes["run2_mssm"] = (
        "("
        + str(int(signal_scale))
        + " pb) gg#phi("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH"],
        False,
    )
    sig_schemes["run2_mssm_bbH"] = (
        "("
        + str(int(signal_scale))
        + " pb) bb#phi("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["bbH"],
        False,
    )
    sig_schemes["sm_cp"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggH_sm_htt"],
        False,
    )
    sig_schemes["sm_ps"] = (
        str(int(signal_scale)) + "#times PS ggH#rightarrow#tau#tau",
        ["ggH_ps_htt"],
        False,
    )
    sig_schemes["sm_cp_decays_ggh"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggH_sm_htt"],
        False,
    )
    sig_schemes["sm_cp_decays_qqh"] = (
        str(int(signal_scale)) + "#times SM qqH#rightarrow#tau#tau",
        ["qqH_sm_htt"],
        False,
    )
    sig_schemes["sm_cp_decays"] = (
        str(int(signal_scale)) + "#times SM H#rightarrow#tau#tau",
        ["ggH_sm_htt", "qqH_sm_htt"],
        False,
    )
    sig_schemes["sm_cp_decays_ps"] = (
        str(int(signal_scale)) + "#times PS H#rightarrow#tau#tau",
        ["ggH_ps_htt", "qqH_ps_htt"],
        False,
    )
    sig_schemes["sm_cp_decays_mm"] = (
        str(int(signal_scale)) + "#times MM H#rightarrow#tau#tau",
        ["ggH_mm_htt", "qqH_mm_htt"],
        False,
    )
    sig_schemes["sm_18"] = (
        str(int(signal_scale)) + "#times SM H#rightarrow#tau#tau",
        ["ggH_ph_htt", "qqH_htt"],
        False,
    )

    plot_signals_dict = {
        "vlq_betaRd33_0_mU2_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=2 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_0_mU2_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=2 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_0_mU2_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=2 TeV, g_{U}=3 (XS)",
        "vlq_betaRd33_0_mU3_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=3 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_0_mU3_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=3 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_0_mU3_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=3 TeV, g_{U}=3 (XS)",
        "vlq_betaRd33_0_mU4_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=4 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_0_mU4_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=4 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_0_mU4_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=0, m_{U}=4 TeV, g_{U}=3 (XS)",
        "vlq_betaRd33_minus1_mU2_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=2 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_minus1_mU2_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=2 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_minus1_mU2_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=2 TeV, g_{U}=3 (XS)",
        "vlq_betaRd33_minus1_mU3_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=3 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_minus1_mU3_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=3 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_minus1_mU3_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=3 TeV, g_{U}=3 (XS)",
        "vlq_betaRd33_minus1_mU4_gU1": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=4 TeV, g_{U}=1 (XS)",
        "vlq_betaRd33_minus1_mU4_gU2": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=4 TeV, g_{U}=2 (XS)",
        "vlq_betaRd33_minus1_mU4_gU3": str(int(signal_scale))
        + "#times VLQ #beta_{R}^{b#tau}=-1, m_{U}=4 TeV, g_{U}=3 (XS)",
    }

    ModTDRStyle(r=0.04, l=0.14)
    R.TGaxis.SetExponentOffset(-0.06, 0.01, "y")

    if ("sm" in signal_scheme and "mssm" not in signal_scheme) or True:
        background_schemes = {
            "mt": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL", "ZJ"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau",
                    ["ZTT", "EWKZ"],
                    R.TColor.GetColor(248, 206, 104),
                ),
            ],
            "et": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrowee", ["ZL", "ZJ"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau",
                    ["ZTT", "EWKZ"],
                    R.TColor.GetColor(248, 206, 104),
                ),
            ],
            "tt": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", [qcdname], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak",
                    ["VVT", "VVJ", "W", "ZL", "ZJ"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau",
                    ["ZTT", "EWKZ"],
                    R.TColor.GetColor(248, 206, 104),
                ),
            ],
            "em": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak", ["VVJ", "VVT", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrowll", ["ZLL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau",
                    ["ZTT", "EWKZ"],
                    R.TColor.GetColor(248, 206, 104),
                ),
            ],
            "zm": [
                backgroundComp(
                    "Misidentified #mu", ["QCD"], R.TColor.GetColor(250, 202, 255)
                ),
                backgroundComp("t#bar{t}", ["TT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VV", "W", "ZJ"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
            ],
            "zmm": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu",
                    ["ZL", "ZJ", "EWKZ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "mm": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu",
                    ["ZL", "ZJ", "EWKZ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "zee": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow ee",
                    ["ZL", "ZJ", "ZTT", "EWKZ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "ee": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow ee",
                    ["ZL", "ZJ", "ZTT", "EWKZ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "w": [backgroundComp("W", ["W"], R.TColor.GetColor(222, 90, 106))],
            "w_shape": [
                backgroundComp(
                    "W loosened shape", ["W_shape"], R.TColor.GetColor(222, 90, 106)
                )
            ],
            "qcd": [backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255))],
            "qcd_shape": [
                backgroundComp(
                    "QCD loosened shape",
                    ["QCD_shape"],
                    R.TColor.GetColor(250, 202, 255),
                )
            ],
        }
    else:
        background_schemes = {
            "mt": [
                backgroundComp("TTT", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("TTJ", ["TTJ"], R.TColor.GetColor(59, 49, 196)),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp("VVT", ["VVT"], R.TColor.GetColor(255, 178, 102)),
                backgroundComp("VVJ", ["VVJ"], R.TColor.GetColor(203, 124, 19)),
                backgroundComp("W", ["W"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp("ZL", ["ZL"], R.TColor.GetColor(100, 192, 232)),
                backgroundComp("ZJ", ["ZJ"], R.TColor.GetColor(14, 97, 132)),
                backgroundComp("ZTT", ["ZTT"], R.TColor.GetColor(248, 206, 104)),
            ],
            "et": [
                backgroundComp("TTT", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("TTJ", ["TTJ"], R.TColor.GetColor(59, 49, 196)),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp("VVT", ["VVT"], R.TColor.GetColor(255, 178, 102)),
                backgroundComp("VVJ", ["VVJ"], R.TColor.GetColor(203, 124, 19)),
                backgroundComp("W", ["W"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp("ZL", ["ZL"], R.TColor.GetColor(100, 192, 232)),
                backgroundComp("ZJ", ["ZJ"], R.TColor.GetColor(14, 97, 132)),
                backgroundComp("ZTT", ["ZTT"], R.TColor.GetColor(248, 206, 104)),
            ],
            "tt": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak",
                    ["VVT", "VVJ", "W", "ZL", "ZJ"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "em": [
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak", ["VVJ", "VVT", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrowll", ["ZLL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "zm": [
                backgroundComp(
                    "Misidentified #mu", ["QCD"], R.TColor.GetColor(250, 202, 255)
                ),
                backgroundComp("t#bar{t}", ["TT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VV", "W", "ZJ"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
            ],
            "zmm": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#mu#mu",
                    ["ZL", "ZJ", "ZTT"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "zee": [
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
                ),
                backgroundComp(
                    "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow ee",
                    ["ZL", "ZJ", "ZTT"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
            "dy": [
                backgroundComp(
                    "Z#rightarrowll", ["ZLL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "w": [backgroundComp("W", ["W"], R.TColor.GetColor(222, 90, 106))],
            "w_shape": [
                backgroundComp(
                    "W loosened shape", ["W_shape"], R.TColor.GetColor(222, 90, 106)
                )
            ],
            "qcd": [backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255))],
            "qcd_shape": [
                backgroundComp(
                    "QCD loosened shape",
                    ["QCD_shape"],
                    R.TColor.GetColor(250, 202, 255),
                )
            ],
            "ff_comp": [
                backgroundComp(
                    "t#bar{t} j#rightarrow#tau",
                    ["TTJ"],
                    R.TColor.GetColor(155, 152, 204),
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak j#rightarrow#tau",
                    ["VVJ", "W", "ZJ"],
                    R.TColor.GetColor(222, 90, 106),
                ),
            ],
        }
    if channel == "zee" or channel == "zmm" or channel == "mm" or channel == "ee":
        background_schemes["dy"] = [
            backgroundComp("DY", ["ZLL"], R.TColor.GetColor(100, 192, 232))
        ]
    if FF and not w_ff_closure and not qcd_ff_closure:
        background_schemes = {
            "mt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "et": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrowee", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "tt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VVT", "ZL"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "jet#rightarrow#tau_{h}",
                    ["jetFakes", "Wfakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                #                    backgroundComp("jet#rightarrow#tau_{h}",["jetFakes"],R.TColor.GetColor(192,232,100)),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "ff_comp": [
                backgroundComp(
                    "t#bar{t} jet#rightarrow#tau_{h}",
                    ["TTJ"],
                    R.TColor.GetColor(155, 152, 204),
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak jet#rightarrow#tau_{h}",
                    ["VVJ", "W"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow ll jet#rightarrow#tau_{h}",
                    ["ZJ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
        }
    elif FF and qcd_ff_closure:
        background_schemes = {
            "mt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "QCD jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Other jet#rightarrow#tau",
                    ["W_res", "VVJ_res", "ZJ_res", "TTJ_res"],
                    R.TColor.GetColor(250, 202, 255),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "et": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrowee", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "QCD jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Other jet#rightarrow#tau",
                    ["W_res", "VVJ_res", "ZJ_res", "TTJ_res"],
                    R.TColor.GetColor(250, 202, 255),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "tt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VVT", "ZL"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "jet#rightarrow#tau_{h}",
                    ["jetFakes", "Wfakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                #                    backgroundComp("jet#rightarrow#tau_{h}",["jetFakes"],R.TColor.GetColor(192,232,100)),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "ff_comp": [
                backgroundComp(
                    "t#bar{t} jet#rightarrow#tau_{h}",
                    ["TTJ"],
                    R.TColor.GetColor(155, 152, 204),
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak jet#rightarrow#tau_{h}",
                    ["VVJ", "W"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow ll jet#rightarrow#tau_{h}",
                    ["ZJ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
        }
    elif FF and w_ff_closure:
        background_schemes = {
            "mt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "W jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Other jet#rightarrow#tau",
                    ["QCD_res", "VVJ_res", "ZJ_res", "TTJ_res"],
                    R.TColor.GetColor(250, 202, 255),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "et": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrowee", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "W jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                backgroundComp(
                    "Other jet#rightarrow#tau",
                    ["QCD_res", "VVJ_res", "ZJ_res", "TTJ_res"],
                    R.TColor.GetColor(250, 202, 255),
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "tt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VVT", "ZL"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "jet#rightarrow#tau_{h}",
                    ["jetFakes", "Wfakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
                #                    backgroundComp("jet#rightarrow#tau_{h}",["jetFakes"],R.TColor.GetColor(192,232,100)),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
            ],
            "ff_comp": [
                backgroundComp(
                    "t#bar{t} jet#rightarrow#tau_{h}",
                    ["TTJ"],
                    R.TColor.GetColor(155, 152, 204),
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak jet#rightarrow#tau_{h}",
                    ["VVJ", "W"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow ll jet#rightarrow#tau_{h}",
                    ["ZJ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
        }

    if embedding:
        background_schemes["zmm"] = [
            backgroundComp(
                "#mu#rightarrow#mu embedding",
                ["EmbedZL"],
                R.TColor.GetColor(100, 192, 232),
            )
        ]
        for chan in ["em", "et", "mt", "tt", "zmm", "zee"]:
            if chan not in background_schemes:
                continue
            schemes = background_schemes[chan]
            print("HI",chan,schemes)
            for bkg in schemes:
                if (
                    chan != "zmm"
                    and chan != "zee"
                    and bkg["leg_text"] == "Z#rightarrow#tau#tau"
                ):
                    bkg["plot_list"] = ["EmbedZTT"]
                    bkg["leg_text"] = "#mu#rightarrow#tau embedding"
                if chan == "zee" and bkg["leg_text"] == "Z#rightarrow ee":
                    bkg["plot_list"] = ["EmbedZL", "ZJ"]

    total_datahist = infile.Get(nodename + "/data_obs").Clone()
    if scheme == "w_shape":
        total_datahist = infile.Get(nodename + "/W").Clone()
    if scheme == "qcd_shape":
        total_datahist = infile.Get(nodename + "/QCD").Clone()
    if scheme == "ff_comp":
        total_datahist = infile.Get(nodename + "/jetFakes").Clone()

    blind_datahist = total_datahist.Clone()
    total_datahist.SetMarkerStyle(20)
    blind_datahist.SetMarkerStyle(20)
    blind_datahist.SetLineColor(1)
    blind_datahist.SetMarkerSize(0.5)

    # Blinding by hand using requested range, set to 200-4000 by default:
    if blind:
        for i in range(0, total_datahist.GetNbinsX()):
            low_edge = total_datahist.GetBinLowEdge(i + 1)
            high_edge = low_edge + total_datahist.GetBinWidth(i + 1)
            if (low_edge > float(x_blind_min) and low_edge < float(x_blind_max)) or (
                high_edge > float(x_blind_min) and high_edge < float(x_blind_max)
            ):
                blind_datahist.SetBinContent(i + 1, 0)
                blind_datahist.SetBinError(i + 1, 0)

    # Create stacked plot for the backgrounds
    bkg_histos = []
    for i, t in enumerate(background_schemes[scheme]):
        plots = t["plot_list"]
        h = R.TH1F()
        for j, k in enumerate(plots):
            if not (
                isinstance(infile.Get(nodename + "/" + k), R.TH1D)
                or isinstance(infile.Get(nodename + "/" + k), R.TH1F)
            ):
                continue

            if h.GetEntries() == 0:
                h = infile.Get(nodename + "/" + k).Clone()
                h.SetName(k)
            else:
                h.Add(infile.Get(nodename + "/" + k).Clone())

        if discrete_x_axis:
            h, new_discrete_x_labels = ConvertHistToDiscreteBins(h, discrete_x_labels)

        h.SetFillColor(t["colour"])
        h.SetLineColor(R.kBlack)
        h.SetMarkerSize(0)
        if norm_bins:
            h.Scale(1.0, "width")
        if h.GetName() == "":
            continue
        bkg_histos.append(h)

    stack = R.THStack("hs", "")
    bkghist = R.TH1F()
    for hists in bkg_histos:
        stack.Add(hists.Clone())
        if bkghist.GetEntries() == 0:
            bkghist = hists.Clone()
        else:
            bkghist.Add(hists.Clone())

    c1 = R.TCanvas()
    c1.cd()

    if ratio and not threePads:
        pads = TwoPadSplit(0.29, 0.025, 0.025)
    elif ratio and threePads:
        pads = MultiRatioSplit([0.25, 0.15], [0.0, 0.01], [0.01, 0.01])
    else:
        pads = OnePad()
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if log_x:
        pads[0].SetLogx(1)
    if custom_x_range:
        if x_axis_max > bkghist.GetXaxis().GetXmax():
            x_axis_max = bkghist.GetXaxis().GetXmax()
    if ratio:
        if log_x:
            pads[1].SetLogx(1)
        if not threePads:
            axish = createAxisHists(
                2,
                bkghist,
                bkghist.GetXaxis().GetXmin(),
                bkghist.GetXaxis().GetXmax() - 0.0001,
            )
            axish[1].GetXaxis().SetTitle(x_title)
            axish[1].GetXaxis().SetLabelSize(0.03)
            axish[1].GetXaxis().SetTitleSize(0.04)
            axish[1].GetYaxis().SetNdivisions(4)
            if scheme == "w_shape" or scheme == "qcd_shape" or scheme == "ff_comp":
                axish[1].GetYaxis().SetTitle("Ratio")
            else:
                axish[1].GetYaxis().SetTitle("Obs/Exp")
            axish[1].GetYaxis().SetTitleOffset(1.6)
            axish[1].GetYaxis().SetTitleSize(0.04)
            axish[1].GetYaxis().SetLabelSize(0.03)

            # centre X axis bin labels for the following variables
            if "decay mode" in x_title or x_title in [
                "N_{jets}",
                "N_{b-jets}",
                "N_{b-tag}^{tight}",
                "N_{b-tag}^{loose}",
            ]:
                axish[0].GetXaxis().SetNdivisions(axish[1].GetNbinsX() + 1)
                axish[1].GetXaxis().SetNdivisions(axish[1].GetNbinsX() + 1)
                axish[1].GetXaxis().CenterLabels(True)

            axish[0].GetXaxis().SetTitleSize(0)
            axish[0].GetXaxis().SetLabelSize(0)
            if log_x:
                axish[1].GetXaxis().SetMoreLogLabels()
                axish[1].GetXaxis().SetNoExponent()
            if custom_x_range:
                axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
                axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            if custom_y_range:
                axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
            if discrete_x_axis:
                for i in range(0, len(new_discrete_x_labels)):
                    axish[1].GetXaxis().SetBinLabel(i + 1, new_discrete_x_labels[i])
                    axish[1].GetXaxis().SetLabelSize(0.05)
                    axish[1].GetYaxis().SetTitleOffset(1.2)
        elif threePads:
            axish = createAxisHists(
                3,
                bkghist,
                bkghist.GetXaxis().GetXmin(),
                bkghist.GetXaxis().GetXmax() - 0.01,
            )

            if bkg_comp:
                axish[1].GetYaxis().SetTitle("Bkg Frac")
            else:
                axish[1].GetYaxis().SetTitle("Obs/Exp")
            axish[1].GetYaxis().SetTitleOffset(1.3)
            axish[1].GetYaxis().SetTitleSize(0.03)
            axish[1].GetYaxis().SetLabelSize(0.03)
            axish[1].GetYaxis().SetNdivisions(4)

            axish[2].GetXaxis().SetTitle(x_title)
            axish[2].GetXaxis().SetLabelSize(0.03)
            axish[2].GetXaxis().SetTitleSize(0.04)
            axish[2].GetYaxis().SetNdivisions(4)
            if bkg_comp:
                axish[2].GetYaxis().SetTitle("Obs/Exp")
            else:
                axish[2].GetYaxis().SetTitle("Purity")
            axish[2].GetYaxis().SetTitleOffset(1.2)
            axish[2].GetYaxis().SetTitleSize(0.03)
            axish[2].GetYaxis().SetLabelSize(0.03)

            axish[0].GetXaxis().SetTitleSize(0)
            axish[0].GetXaxis().SetLabelSize(0)
            axish[1].GetXaxis().SetTitleSize(0)
            axish[1].GetXaxis().SetLabelSize(0)
            if custom_x_range:
                axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
                axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
                axish[2].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            if custom_y_range:
                axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)

    else:
        axish = createAxisHists(
            1,
            bkghist,
            bkghist.GetXaxis().GetXmin(),
            bkghist.GetXaxis().GetXmax() - 0.01,
        )
        axish[0].GetXaxis().SetTitle(x_title)
        axish[0].GetXaxis().SetTitleSize(0.04)
        axish[0].GetXaxis().SetLabelSize(0.03)
        if log_x:
            axish[0].GetXaxis().SetMoreLogLabels()
            axish[0].GetXaxis().SetNoExponent()
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)
    if not ratio:
        axish[0].GetXaxis().SetLabelSize(0.03)
    if custom_y_range and log_y:
        axish[0].SetMinimum(0.1)
        axish[0].SetMaximum(
            10
            ** (
                (1 + extra_pad)
                * (
                    math.log10(
                        1.1 * bkghist.GetMaximum() - math.log10(axish[0].GetMinimum())
                    )
                )
            )
        )

    if not custom_y_range:
        if log_y:
            axish[0].SetMinimum(0.1)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * bkghist.GetMaximum()
                            - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            axish[0].SetMinimum(0)
            axish[0].SetMaximum(1.1 * (1 + extra_pad) * bkghist.GetMaximum())
    axish[0].Draw()

    # Draw uncertainty band
    bkghist.SetFillColor(CreateTransparentColor(12, 0.4))
    bkghist.SetLineColor(CreateTransparentColor(12, 0.4))
    bkghist.SetMarkerSize(0)
    bkghist.SetMarkerColor(CreateTransparentColor(12, 0.4))

    sighist = R.TH1F()
    # for single signal scheme
    if signal_mass != "":
        sig_scheme = sig_schemes[signal_scheme]
        for i in sig_scheme[1]:
            h = infile.Get(nodename + "/" + i + signal_mass).Clone()
            if sighist.GetEntries() == 0:
                sighist = h
            else:
                sighist.Add(h)
        sighist.SetLineColor(R.kRed)
        sighist.SetLineWidth(3)
        sighist.Scale(signal_scale)
        if norm_bins:
            sighist.Scale(1.0, "width")
        if not split_sm_scheme:
            if sig_scheme[2]:
                stack.Add(sighist.Clone())
                if not custom_y_range:
                    axish[0].SetMaximum(1.1 * (1 + extra_pad) * stack.GetMaximum())
            stack.Draw("histsame")
            if not sig_scheme[2]:
                sighist.Draw("histsame")
    elif not split_sm_scheme:
        stack.Draw("histsame")

    colours = [2, 3, 4, 5, 6, 7, 8]
    h = []
    h_ratio = []
    pads[0].cd()
    if plot_signals != [""]:
        for i in range(0, len(plot_signals)):
            h.append(infile.Get(nodename + "/" + plot_signals[i]).Clone())
            if norm_bins:
                h[i].Scale(1.0, "width")
            h[i].SetLineColor(colours[i])
            h[i].SetLineWidth(2)
            h[i].Scale(signal_scale)
            h[i].Draw("histsame")

    # separate out signal into powheg ggH and qqH,
    # JHU ggH, and VH

    if split_sm_scheme:
        sighists = dict()

        signal_split_schemes = ""
        if signal_scheme == "run2_mssm":
            signal_split_schemes = ["run2_mssm", "run2_mssm_bbH"]
        elif ggh_scheme == "powheg":
            signal_split_schemes = ["sm_ggH", "sm_qqH"]
            # signal_split_schemes = ['sm_cp_decays','sm_cp_decays_ps']
        elif ggh_scheme == "tauspinner":
            signal_split_schemes = [
                "sm_cp_decays",
                "sm_cp_decays_ps",
                "sm_cp_decays_mm",
            ]
            # signal_split_schemes = ['sm_cp_decays_ggh','sm_cp_decays_qqh']
        elif ggh_scheme == "JHU":
            signal_split_schemes = ["sm_ggH_JHU", "sm_qqH", "sm_VH"]
        if ggh_scheme == "madgraph":
            signal_split_schemes = ["sm_cp", "ps_sm"]
            # signal_split_schemes = ['sm_cp','sm_ps','sm_mm']

        for index, split_scheme in enumerate(signal_split_schemes):
            sighists[split_scheme] = R.TH1F()
            if signal_mass != "":
                sig_scheme = sig_schemes[split_scheme]
                for i in sig_scheme[1]:
                    h = infile.Get(nodename + "/" + i + signal_mass).Clone()
                    if sighists[split_scheme].GetEntries() == 0:
                        sighists[split_scheme] = h
                    else:
                        sighists[split_scheme].Add(h)

                if split_scheme in [
                    "sm_cp",
                    "sm_ggH",
                    "sm_cp_decays",
                    "sm_cp_decays_ggh",
                    "run2_mssm",
                ]:
                    sighists[split_scheme].SetLineColor(R.kRed)
                elif split_scheme in ["sm_qqH", "sm_cp_decays_qqh", "run2_mssm_bbH"]:
                    sighists[split_scheme].SetLineColor(R.TColor.GetColor(51, 51, 230))
                elif split_scheme in ["sm_ps", "sm_cp_decays_ps"]:
                    sighists[split_scheme].SetLineColor(R.kGreen + 3)
                elif split_scheme in ["sm_mm", "sm_cp_decays_mm"]:
                    sighists[split_scheme].SetLineColor(R.kOrange - 5)

                sighists[split_scheme].SetLineWidth(3)
                sighists[split_scheme].Scale(signal_scale)
                if norm_bins:
                    sighists[split_scheme].Scale(1.0, "width")
                if sig_scheme[2]:
                    stack.Add(sighists[split_scheme].Clone())
                    if not custom_y_range:
                        axish[0].SetMaximum(1.1 * (1 + extra_pad) * stack.GetMaximum())
                if index == 0:
                    stack.Draw("histsame")
                if not sig_scheme[2]:
                    sighists[split_scheme].Draw("histsame")

            else:
                stack.Draw("histsame")

    # Auto blinding option
    if auto_blind:
        for i in range(1, total_datahist.GetNbinsX() + 1):
            b = bkghist.GetBinContent(i)
            s = sighist.GetBinContent(i) / signal_scale
            if PassAutoBlindMetric(s, b, metric=0.5):
                blind_datahist.SetBinContent(i, 0)
                blind_datahist.SetBinError(i, 0)

    if norm_bins:
        blind_datahist.Scale(1.0, "width")
        total_datahist.Scale(1.0, "width")

    ## Add another signal mass point
    # sighist2 = R.TH1F()
    # sighist2= infile.Get(nodename+'/ggH350').Clone()
    # sighist2.SetLineColor(R.kRed)
    # sighist2.SetLineWidth(3)
    # sighist2.Scale(signal_scale)
    ## Add another signal mass point
    # sighist3 = R.TH1F()
    # sighist3= infile.Get(nodename+'/ggH700').Clone()
    # sighist3.SetLineColor(R.kGreen+2)
    # sighist3.SetLineWidth(3)
    # sighist3.Scale(signal_scale)
    # if norm_bins:
    #    sighist2.Scale(1.0,"width")
    #    sighist3.Scale(1.0,"width")
    # sighist2.Draw("histsame")
    # sighist3.Draw("histsame")
    error_hist = bkghist.Clone()
    if do_custom_uncerts:
        bkg_uncert_up = infile.Get(nodename + "/" + custom_uncerts_up_name).Clone()
        bkg_uncert_down = infile.Get(nodename + "/" + custom_uncerts_down_name).Clone()

        if discrete_x_axis:
            bkg_uncert_up, _ = ConvertHistToDiscreteBins(
                bkg_uncert_up, discrete_x_labels
            )
            bkg_uncert_down, _ = ConvertHistToDiscreteBins(
                bkg_uncert_up, discrete_x_labels
            )

        if norm_bins:
            bkg_uncert_up.Scale(1.0, "width")
            bkg_uncert_down.Scale(1.0, "width")

        for i in range(1, bkg_uncert_up.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            bin_up = bkg_uncert_up.GetBinContent(i)
            bin_down = bkg_uncert_down.GetBinContent(i)
            error = abs(bin_up - bin_down) / 2
            band_center = max(bin_up, bin_down) - error
            if add_stat_to_syst:
                error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinContent(i, band_center)
            error_hist.SetBinError(i, error)

    if add_flat_uncert > 0:
        for i in range(1, error_hist.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            error = add_flat_uncert * error_hist.GetBinContent(i)
            error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinError(i, error)

    if discrete_x_axis:
        blind_datahist, _ = ConvertHistToDiscreteBins(blind_datahist, discrete_x_labels)

    error_hist.Draw("e2same")
    if not blind:
        blind_datahist.Draw("E same")
    axish[0].Draw("axissame")

    # Setup legend
    legend = PositionedLegend(0.2,0.2,3,0.01)
    # legend = PositionedLegend(0.4, 0.15, 3, 0.005)  # when showing plots of signal
    legend.SetTextFont(42)
    legend.SetTextSize(0.019)
    # legend.SetTextSize(0.018)  # when showing plots of signal
    legend.SetFillColor(0)
    if scheme == "w_shape" or scheme == "qcd_shape":
        legend.AddEntry(blind_datahist, "un-loosened shape", "PE")
    elif scheme == "ff_comp":
        legend.AddEntry(blind_datahist, "FF jet#rightarrow#tau_{h}", "PE")
    elif not blind:
        legend.AddEntry(blind_datahist, "Observation", "PE")
    # Drawn on legend in reverse order looks better
    bkg_histos.reverse()
    background_schemes[scheme].reverse()
    for legi, hists in enumerate(bkg_histos):
        legend.AddEntry(hists, background_schemes[scheme][legi]["leg_text"], "f")
    if do_custom_uncerts and uncert_title != "":
        legend.AddEntry(error_hist, uncert_title, "f")
    else:
        legend.AddEntry(error_hist, "#splitline{Background}{Uncertainty}", "f")
    if plot_signals != [""]:
        for i in range(0, len(plot_signals)):
            if plot_signals[i] in plot_signals_dict.keys():
                legend.AddEntry(h[i], plot_signals_dict[plot_signals[i]], "l")
            else:
                legend.AddEntry(h[i], plot_signals[i], "l")

    if signal_mass != "":
        if not split_sm_scheme:
            legend.AddEntry(sighist, sig_schemes[signal_scheme][0], "l")
        else:
            for split_scheme in signal_split_schemes:
                legend.AddEntry(
                    sighists[split_scheme], sig_schemes[split_scheme][0], "l"
                )
                # if split_scheme != "sm_cp":
                #     ks_score =  sighists[split_scheme].KolmogorovTest(sighists["sm_cp"],"DNOU")
                #     legend.AddEntry(R.TObject(), "K-S ({}) = {:.3f}".format(split_scheme,ks_score), "")

    ## Add a second signal mass
    # legend.AddEntry(sighist2,str(int(signal_scale))+"#times gg#phi(350 GeV)#rightarrow#tau#tau","l")
    # legend.AddEntry(sighist3,str(int(signal_scale))+"#times gg#phi(700 GeV)#rightarrow#tau#tau","l")
    if scheme == "qcd_shape" or scheme == "w_shape":
        ks_score = error_hist.KolmogorovTest(total_datahist)
        legend.AddEntry(R.TObject(), "K-S probability = %.3f" % ks_score, "")
    legend.Draw("same")
    if channel == "em":
        channel_label = "e_{}#mu_{}"
    if channel == "et":
        channel_label = "e_{}#tau_{h}"
    if channel == "mt":
        channel_label = "#mu_{}#tau_{h}"
    if channel == "tt":
        channel_label = "#tau_{h}#tau_{h}"
    if channel == "zmm" or channel == "mm":
        channel_label = "Z#rightarrow#mu#mu"
    if channel == "zee" or channel == "ee":
        channel_label = "Z#rightarrow ee"
    if "MVA" in x_title or "NN" in x_title:
        channel_label += " {}".format(cat)
    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(0.145, 0.955, channel_label)

    # CMS and lumi labels
    if not custom_y_range:
        FixTopRange(pads[0], GetPadYMax(pads[0]), extra_pad if extra_pad > 0 else 0.30)
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.15, 0.05, 1.0, "", 1.0)
    # DrawCMSLogo(pads[0], 'CMS', '', 11, 0.045, 0.05, 1.0, '', 1.0)
    DrawTitle(pads[0], lumi, 3)

    if ratio:
        if bkg_comp:
            pad_shift = 1
        else:
            pad_shift = 0
        ratio_bkghist = MakeRatioHist(error_hist.Clone(), bkghist.Clone(), True, False)
        blind_ratio = MakeRatioHist(
            blind_datahist.Clone(), bkghist.Clone(), True, False
        )
        pads[1 + pad_shift].cd()
        pads[1 + pad_shift].SetGrid(0, 1)
        axish[1 + pad_shift].Draw("axis")
        if not ratio_log_y and ratio_range != "auto":
            if ratio_range == "0,2":
                ratio_range = "0.01,1.99"
            axish[1 + pad_shift].SetMinimum(float(ratio_range.split(",")[0]))
            axish[1 + pad_shift].SetMaximum(float(ratio_range.split(",")[1]))
        if plot_signals != [""]:
            for i in range(0, len(plot_signals)):
                h_ratio.append(h[i].Clone())
                h_ratio[i].Scale(1 / signal_scale)
                h_ratio[i].Add(bkghist)
                h_ratio[i].Divide(bkghist)
                h_ratio[i].SetLineColor(colours[i])
                h_ratio[i].SetLineWidth(2)
                h_ratio[i].Draw("histsame")
        ratio_bkghist.SetMarkerSize(0)
        ratio_bkghist.Draw("e2same")
        if ratio_range == "auto":
            max_list, min_list = [], []
            max_list.append(ratio_bkghist.GetMaximum())
            min_list.append(NonZeroMinimum(ratio_bkghist))
            # max_list.append(blind_ratio.GetMaximum())
            # min_list.append(blind_ratio.GetMinimum())
            for i in h_ratio:
                max_list.append(i.GetMaximum())
                min_list.append(NonZeroMinimum(i))
            if min(min_list) <= 0:
                min_val = 0.01
            else:
                min_val = min(min_list)
            if max(max_list) > 100:
                max_val = 99.9
            else:
                max_val = max(max_list)
            print(min_list, min(min_list), min_val)
            print(max_list, min(max_list), max_val)
            axish[1 + pad_shift].SetMinimum(
                round(min_val, -int(floor(log10(abs(min_val)))))
            )
            if round(max_val, -int(floor(log10(abs(max_val))))) < max_val:
                axish[1 + pad_shift].SetMaximum(
                    round(max_val, -int(floor(log10(abs(max_val)))))
                    + 10 ** int(floor(log10(abs(max_val))))
                )
            else:
                axish[1 + pad_shift].SetMaximum(
                    round(max_val, -int(floor(log10(abs(max_val)))))
                )

        if not blind:
            blind_ratio.DrawCopy("e0same")
        if ratio_log_y:
            pads[1].SetLogy(1)
            axish[1 + pad_shift].SetMinimum(0.9)
            axish[1 + pad_shift].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * ratio_bkghist.GetMaximum()
                            - math.log10(axish[1 + pad_shift].GetMinimum())
                        )
                    )
                )
            )
        pads[1].RedrawAxis("G")

        ## lines below will show seperate lines indicating systematic up and down bands
        if do_custom_uncerts:
            # bkg_uncert_up.SetLineColor(R.TColor.GetColor("#1f78b4"))
            # bkg_uncert_down.SetLineColor(R.TColor.GetColor("#ff7f00"))
            # if discrete_x_axis:
            #  bkg_uncert_up_clone = bkg_uncert_up.Clone()
            #  bkg_uncert_down_clone = bkg_uncert_down.Clone()
            #  bkg_uncert_up = R.TH1F(bkg_uncert_up_clone.GetName(),"",len(new_bins)-1, new_bins)
            #  bkg_uncert_down = R.TH1F(bkg_uncert_down_clone.GetName(),"",len(new_bins)-1, new_bins)
            #  for r in range(1,bkg_uncert_up_clone.GetNbinsX()+2):
            #     bkg_uncert_up.SetBinContent(r,bkg_uncert_up_clone.GetBinContent(r))
            #     bkg_uncert_up.SetBinError(r,bkg_uncert_up_clone.GetBinError(r))
            #     bkg_uncert_down.SetBinContent(r,bkg_uncert_down_clone.GetBinContent(r))
            #     bkg_uncert_down.SetBinError(r,bkg_uncert_down_clone.GetBinError(r))

            bkg_uncert_up.SetLineColor(CreateTransparentColor(12, 0.4))
            bkg_uncert_down.SetLineColor(CreateTransparentColor(12, 0.4))
            bkg_uncert_up.SetLineWidth(0)
            bkg_uncert_down.SetLineWidth(0)
            bkg_uncert_up = MakeRatioHist(bkg_uncert_up, bkghist.Clone(), True, False)
            bkg_uncert_down = MakeRatioHist(
                bkg_uncert_down, bkghist.Clone(), True, False
            )
            bkg_uncert_up.Draw("histsame")
            bkg_uncert_down.Draw("histsame")

        pads[1 + pad_shift].RedrawAxis("G")

        # Draw uncertainty band
        # these lines will highlight bins with bb uncertainty > 90%
        large_uncert_hist = ratio_bkghist.Clone()
        no_large_uncerts = True
        for i in range(1, large_uncert_hist.GetNbinsX() + 1):
            if large_uncert_hist.GetBinError(i) <= 0.9:
                large_uncert_hist.SetBinError(i, 0)
            else:
                no_large_uncerts = False
        large_uncert_hist.SetFillColor(CreateTransparentColor(2, 0.4))
        large_uncert_hist.SetLineColor(CreateTransparentColor(2, 0.4))
        large_uncert_hist.SetMarkerSize(0)
        large_uncert_hist.SetMarkerColor(CreateTransparentColor(2, 0.4))
        # if not no_large_uncerts: large_uncert_hist.Draw("e2same")

    if threePads:
        if bkg_comp:
            # Create stacked plot for the backgrounds
            fracstack = R.THStack("hs", "")
            for hists in reversed(bkg_histos):
                hists.Divide(bkghist)
                fracstack.Add(hists)
            pads[1].cd()
            if log_x:
                pads[1].SetLogx(1)
                pads[2].SetLogx(1)
                axish[2].GetXaxis().SetMoreLogLabels()
                axish[2].GetXaxis().SetNoExponent()

            pads[1].SetGrid(0, 1)
            axish[1].Draw("axis")
            axish[1].SetMinimum(0.0)
            axish[1].SetMaximum(1.0)
            fracstack.Draw("histsame")
        else:
            purity_hist = MakePurityHist(sighist.Clone(), stack, cat)
            pads[2].cd()
            pads[2].SetGrid(0, 1)
            axish[2].Draw("axis")
            axish[2].SetMinimum(0.0)
            axish[2].SetMaximum(1.0)
            purity_hist.Draw("histsame")

    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.SetCanvasSize(1100, 1000)
    c1.SetWindowSize(1100, 1000 + (1100 - c1.GetWw()))
    c1.SaveAs(plot_name + ".png")
    c1.Close()


def CompareSysts(hists=[], plot_name="plot", label=""):
    legend_titles = ["nominal", "up", "down"]
    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    colourlist = [R.kBlue, R.kRed, R.kGreen + 3]
    hs = R.THStack("hs", "")
    hist_count = 0

    for hist in hists:
        h = hist
        h.SetFillColor(0)
        h.SetLineWidth(3)
        h.SetLineColor(colourlist[hist_count])
        h.SetMarkerSize(0)
        hs.Add(h)
        hist_count += 1

    c1 = R.TCanvas()
    c1.cd()

    if hists[0].Integral() > 0:
        uncert_up = hists[1].Integral() / hists[0].Integral()
        uncert_down = hists[2].Integral() / hists[0].Integral()
    else:
        uncert_up = 0
        uncert_down = 0

    pads = TwoPadSplit(0.49, 0.01, 0.01)

    axish = createAxisHists(
        2, hists[0], hists[0].GetXaxis().GetXmin(), hists[0].GetXaxis().GetXmax() - 0.01
    )
    # axish[1].GetXaxis().SetTitle("")
    axish[1].GetXaxis().SetLabelSize(0.03)
    axish[1].GetYaxis().SetNdivisions(4)
    axish[1].GetYaxis().SetTitle("Ratio")
    axish[1].GetYaxis().SetTitleOffset(1.6)
    axish[1].GetYaxis().SetTitleSize(0.04)
    axish[1].GetYaxis().SetLabelSize(0.03)

    axish[0].GetXaxis().SetTitleSize(0)
    axish[0].GetXaxis().SetLabelSize(0)
    axish[0].GetYaxis().SetTitle("")
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)

    axish[0].SetMinimum(0)
    axish[0].SetMaximum(1.1 * hs.GetMaximum("nostack"))
    axish[0].Draw()

    hs.Draw("nostack hist same")

    axish[0].Draw("axissame")

    # Setup legend
    legend = PositionedLegend(0.3, 0.2, 3, 0.03)
    legend.SetTextFont(42)
    legend.SetTextSize(0.022)
    legend.SetFillColor(0)

    for legi, hist in enumerate(hists):
        legend.AddEntry(hist, legend_titles[legi], "l")
    legend.Draw("same")

    # CMS label and title
    FixTopRange(pads[0], axish[0].GetMaximum(), 0.2)

    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.145, 0.955, label)

    latex3 = R.TLatex()
    latex3.SetNDC()
    latex3.SetTextAngle(0)
    latex3.SetTextColor(R.kBlack)
    latex3.SetTextSize(0.028)
    lnN_uncert = "lnN uncert = ^{+ %.3f}_{- %.3f}" % (uncert_up, uncert_down)
    latex3.DrawLatex(0.17, 0.855, lnN_uncert)

    ks_hist_nom = hists[0].Clone()
    for i in range(1, ks_hist_nom.GetNbinsX() + 1):
        ks_hist_nom.SetBinError(i, 0)
    ks_up = hists[1].KolmogorovTest(ks_hist_nom)
    ks_down = hists[2].KolmogorovTest(ks_hist_nom)

    latex4 = R.TLatex()
    latex4.SetNDC()
    latex4.SetTextAngle(0)
    latex4.SetTextColor(R.kBlack)
    latex4.SetTextSize(0.028)

    latex4.DrawLatex(0.65, 0.7, "KS up = %.3f" % ks_up)
    latex4.DrawLatex(0.65, 0.65, "KS down = %.3f" % ks_down)

    ratio_hs = R.THStack("ratio_hs", "")
    hist_count = 0
    pads[1].cd()
    pads[1].SetGrid(0, 1)
    axish[1].Draw("axis")

    for hist in hists:
        h = hist.Clone()
        h.SetFillColor(0)
        h.SetLineWidth(3)
        h.SetLineColor(colourlist[hist_count])
        h.SetMarkerSize(0)
        h.SetMarkerColor(h.GetLineColor())
        for i in range(1, h.GetNbinsX() + 1):
            nom_cont = hists[0].GetBinContent(i)
            old_cont = hist.GetBinContent(i)
            old_error = hist.GetBinError(i)
            if nom_cont == 0:
                new_cont = 0
                new_error = 0
            else:
                new_cont = old_cont / nom_cont
                new_error = old_error / nom_cont
            h.SetBinContent(i, new_cont)
            h.SetBinError(i, new_error)
        ratio_hs.Add(h)
        hist_count += 1
    ratio_hs.Draw("nostack l same")

    ratio_min = 0.0
    for i in range(1, hists[0].GetNbinsX() + 1):
        bin_contents = []
        if hists[0].GetBinContent(i) > 0:
            if hists[1].GetBinContent(i) > 0:
                bin_contents.append(
                    hists[1].GetBinContent(i) / hists[0].GetBinContent(i)
                )
            if hists[2].GetBinContent(i) > 0:
                bin_contents.append(
                    hists[2].GetBinContent(i) / hists[0].GetBinContent(i)
                )
        bin_min = 0
        if len(bin_contents) > 0:
            bin_min = min(bin_contents)

        if (bin_min < ratio_min and bin_min > 0) or ratio_min == 0:
            ratio_min = bin_min

    axish[1].SetMinimum(max(1.1 * (ratio_min - 1.0) + 1.0, 0.5))
    axish[1].SetMaximum(min(1.1 * (ratio_hs.GetMaximum("nostack") - 1.0) + 1.0, 1.5))

    pads[1].RedrawAxis("G")
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.Close()


def CompareHists(
    hists=[],
    legend_titles=[],
    title="",
    ratio=True,
    log_y=False,
    log_x=False,
    ratio_range="0.7,1.3",
    custom_x_range=False,
    x_axis_max=4000,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=4000,
    y_axis_min=0,
    x_title="",
    y_title="",
    extra_pad=0,
    norm_hists=False,
    plot_name="plot",
    label="",
    norm_bins=True,
    uncert_hist=None,
    uncert_title="",
    ReweightPlot=False,
    output_file=None,
):

    objects = []
    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    colourlist = [
        R.kBlue,
        R.kRed,
        R.kGreen + 3,
        R.kBlack,
        R.kYellow + 2,
        R.kOrange,
        R.kCyan + 3,
        R.kMagenta + 2,
        R.kViolet - 5,
        R.kGray,
    ]
    if ReweightPlot:
        colourlist = [
            R.kBlack,
            R.kBlue,
            R.kRed,
            R.kGreen + 3,
            R.kYellow + 2,
            R.kOrange,
            R.kCyan + 3,
            R.kMagenta + 2,
            R.kViolet - 5,
            R.kGray,
        ]

    hs = R.THStack("hs", "")
    hist_count = 0
    legend_hists = []
    if isinstance(uncert_hist, (list,)):
        for i in uncert_hist:
            if i is None:
                continue
            if norm_bins and i is not None:
                i.Scale(1.0, "width")
    else:
        if norm_bins and uncert_hist is not None:
            uncert_hist.Scale(1.0, "width")

    for hist in hists:
        # print hist.GetName()
        # for bin_ in range(1,hist.GetNbinsX()):
        #     print hist.GetBinContent(bin_)
        #     print np.sqrt(hist.GetBinContent(bin_))
        if norm_hists:
            hist.Scale(1.0 / hist.Integral(0, hist.GetNbinsX() + 1))
        if norm_bins:
            hist.Scale(1.0, "width")
        h = hist.Clone()
        objects.append(h)
        h.SetFillColor(0)
        h.SetLineWidth(3)
        h.SetLineColor(colourlist[hist_count])
        h.SetMarkerColor(colourlist[hist_count])
        h.SetMarkerSize(0)
        hs.Add(h)
        hist_count += 1
        o = h.Clone()
        objects.append(o)
        legend_hists.append(o)
    # hs.Draw("nostack")

    c1 = R.TCanvas()
    c1.cd()

    if ratio:
        if ReweightPlot:
            pads = TwoPadSplit(0.39, 0.01, 0.01)
        else:
            pads = TwoPadSplit(0.29, 0.01, 0.01)
    else:
        pads = OnePad()
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if log_x:
        pads[0].SetLogx(1)
    if custom_x_range:
        if x_axis_max > hists[0].GetXaxis().GetXmax():
            x_axis_max = hists[0].GetXaxis().GetXmax()
    if ratio:
        if log_x:
            pads[1].SetLogx(1)
        axish = createAxisHists(
            2,
            hists[0],
            hists[0].GetXaxis().GetXmin(),
            hists[0].GetXaxis().GetXmax() - 0.01,
        )
        axish[1].GetXaxis().SetTitle(x_title)
        axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetYaxis().SetNdivisions(4)
        axish[1].GetYaxis().SetTitle("Ratio")
        # if ReweightPlot:
        #  axish[1].GetYaxis().SetTitle("Correction")
        axish[1].GetYaxis().SetTitleOffset(1.6)
        axish[1].GetYaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetLabelSize(0.03)

        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
            axish[1].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    else:
        axish = createAxisHists(
            1,
            hists[0],
            hists[0].GetXaxis().GetXmin(),
            hists[0].GetXaxis().GetXmax() - 0.005,
        )
        axish[0].GetXaxis().SetLabelSize(0.03)
        axish[0].GetXaxis().SetTitle(x_title)
        axish[0].GetXaxis().SetTitleSize(0.04)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)

    hs.Draw("nostack same")

    uncert_hs = R.THStack()
    if uncert_hist is not None:
        if isinstance(uncert_hist, (list,)):
            # col_list = [12,6,4,2,3,4]
            col_list = colourlist
            count = 0
            for i in uncert_hist:
                if i is not None:
                    i.SetFillColor(CreateTransparentColor(col_list[count], 0.4))
                    i.SetLineColor(CreateTransparentColor(col_list[count], 0.4))
                    i.SetMarkerSize(0)
                    i.SetMarkerColor(CreateTransparentColor(col_list[count], 0.4))
                    i.SetFillStyle(1111)
                    uncert_hs.Add(i)
                count += 1
            uncert_hs.Draw("nostack e2same")
        else:
            uncert_hist.SetFillColor(CreateTransparentColor(12, 0.4))
            uncert_hist.SetLineColor(CreateTransparentColor(12, 0.4))
            uncert_hist.SetMarkerSize(0)
            uncert_hist.SetMarkerColor(CreateTransparentColor(12, 0.4))
            uncert_hist.SetFillStyle(1111)
            uncert_hs.Add(uncert_hist)
            uncert_hs.Draw("e2same")

        uncert_hs.Draw("nostack e2same")
    if not custom_y_range:
        if log_y:
            if hs.GetMinimum("nostack") > 0:
                axish[0].SetMinimum(hs.GetMinimum("nostack"))
            else:
                axish[0].SetMinimum(0.0009)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * hs.GetMaximum("nostack")
                            - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            maxi = (
                1.1
                * (1 + extra_pad)
                * max(hs.GetMaximum("nostack"), uncert_hs.GetMaximum("nostack"))
            )
            if not ReweightPlot:
                axish[0].SetMinimum(0)
            else:
                mini = None
                maxi = None
                for h in hists + uncert_hist:
                    if h is None:
                        continue
                    for i in range(1, h.GetNbinsX() + 1):
                        lo = h.GetBinContent(i) - h.GetBinError(i)
                        hi = h.GetBinContent(i) + h.GetBinError(i)
                        if mini is None:
                            mini = min(lo, hi)
                            maxi = max(lo, hi)
                        else:
                            mini = min(mini, lo, hi)
                            maxi = max(maxi, lo, hi)
                # mini = min(hs.GetMinimum("nostack"),uncert_hs.GetMinimum("nostack"))
                mini -= abs(mini) * extra_pad
                axish[0].SetMinimum(mini)
                maxi *= 1.0 + extra_pad
            axish[0].SetMaximum(maxi)
    axish[0].Draw()
    uncert_hs.Draw("nostack e2same")

    hs.Draw("nostack hist same")
    axish[0].Draw("axissame")

    # Setup legend
    tot = len(hists)
    if isinstance(uncert_hist, list):
        tot += len(uncert_hist)
    if tot > 4:
        legend = PositionedLegend(0.35, 0.3, 3, 0.03)
    else:
        legend = PositionedLegend(0.3, 0.15, 3, 0.03)
    legend.SetTextFont(42)
    legend.SetTextSize(0.025)
    legend.SetFillColor(0)

    for legi, hist in enumerate(legend_hists):
        legend.AddEntry(hist, legend_titles[legi], "l")
    if isinstance(uncert_hist, (list,)):
        count = 0
        for i in uncert_hist:
            if i is not None:
                legend.AddEntry(i, uncert_title[count], "f")
            count += 1
    else:
        if uncert_hist is not None and uncert_title:
            legend.AddEntry(uncert_hist, uncert_title, "f")
    legend.Draw("same")

    # CMS label and title
    # FixTopRange(pads[0], axish[0].GetMaximum(), extra_pad if extra_pad>0 else 0.30)
    # DrawCMSLogo(pads[0], 'CMS', 'Preliminary', 11, 0.045, 0.05, 1.0, '', 1.0)
    # DrawCMSLogo(pads[0], 'CMS', 'Simulation', 11, 0.045, 0.05, 1.0, '', 1.0)
    DrawTitle(pads[0], title, 3)

    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.145, 0.955, label)

    # Add ratio plot if required
    if ratio:
        ratio_hs = R.THStack("ratio_hs", "")
        hist_count = 0
        pads[1].cd()
        pads[1].SetGrid(0, 1)
        axish[1].Draw("axis")
        axish[1].SetMinimum(float(ratio_range.split(",")[0]))
        axish[1].SetMaximum(float(ratio_range.split(",")[1]))
        div_hist = hists[0].Clone()
        objects.append(div_hist)

        for i in range(0, div_hist.GetNbinsX() + 2):
            div_hist.SetBinError(i, 0)
        first_hist = True
        for hist in hists:
            h = hist.Clone()
            objects.append(h)

            h.SetFillColor(0)
            h.SetLineWidth(3)
            h.SetLineColor(colourlist[hist_count])
            h.SetMarkerColor(colourlist[hist_count])
            h.SetMarkerSize(0)

            h.Divide(div_hist)
            # if first_hist:
            #    for i in range(1,h.GetNbinsX()+1): h.SetBinError(i,0.00001)
            #    first_hist=False
            o = h.Clone()
            objects.append(o)
            ratio_hs.Add(o)
            hist_count += 1
        if uncert_hist is not None:
            if isinstance(uncert_hist, (list,)):
                ratio_err_hs = R.THStack("ratio_err_hs", "")
                count = 0
                for i in uncert_hist:
                    if i is not None:
                        h = i.Clone()
                        objects.append(h)
                        h.Divide(div_hist)
                        ratio_err_hs.Add(h)
                        h.Draw("e2same")
                    count += 1
                # ratio_err_hs.Draw("nostack e2same")
            else:
                h = uncert_hist.Clone()
                objects.append(h)
                h.Divide(div_hist)
                h.Draw("e2same")
        ratio_hs.Draw("nostack e same")
        pads[1].RedrawAxis("G")
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")
    if output_file is not None:
        output_file.WriteObject(c1, plot_name)
    for o in objects:
        o.IsA().Destructor(o)
    c1.Close()


def HTTPlotSignal(
    nodename,
    infile=None,
    signal_scale=1,
    signal_mass="",
    norm_bins=True,
    channel="mt",
    blind=False,
    x_blind_min=0,
    x_blind_max=4000,
    ratio=True,
    log_y=False,
    log_x=False,
    ratio_range="0.7,1.3",
    custom_x_range=False,
    x_axis_max=4000,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=4000,
    y_axis_min=0,
    x_title="",
    y_title="",
    extra_pad=0,
    signal_scheme="run2_mssm",
    do_custom_uncerts=False,
    add_stat_to_syst=False,
    add_flat_uncert=False,
    uncert_title="uncertainty",
    lumi="35.9",
    plot_name="htt_plot",
    custom_uncerts_up_name="total_bkg_custom_uncerts_up",
    custom_uncerts_down_name="total_bkg_custom_uncerts_down",
):
    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    # Define signal schemes here
    sig_schemes = {}
    sig_schemes["sm_default"] = (
        str(int(signal_scale))
        + "#times SM H("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH", "qqH"],
    )
    sig_schemes["smsummer16"] = (
        str(int(signal_scale))
        + "#times SM H("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH_htt", "qqH_htt"],
    )
    sig_schemes["run2_mssm"] = (
        str(int(signal_scale))
        + "#times gg#phi("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["ggH"],
    )
    sig_schemes["run2_mssm_bbH"] = (
        str(int(signal_scale))
        + "#times bb#phi("
        + signal_mass
        + " GeV)#rightarrow#tau#tau",
        ["bbH"],
    )

    ModTDRStyle(r=0.04, l=0.14)

    stack = R.THStack("hs", "")
    sighist = R.TH1F()
    sighist_copy = R.TH1F()
    if signal_mass != "":
        sig_scheme = sig_schemes[signal_scheme]
        for i in sig_scheme[1]:
            h = infile.Get(nodename + "/" + i + signal_mass).Clone()
            if sighist.GetEntries() == 0:
                sighist = h
            else:
                sighist.Add(h)
        sighist.SetLineColor(R.kBlue)
        sighist.SetLineWidth(3)
        sighist.SetMarkerSize(0)
        sighist.Scale(signal_scale)
        if norm_bins:
            sighist.Scale(1.0, "width")
        stack.Add(sighist.Clone())

        sighist_copy = sighist.Clone()
        stack.Draw("histsame")

    c1 = R.TCanvas()
    c1.cd()

    if ratio:
        pads = TwoPadSplit(0.29, 0.01, 0.01)
    else:
        pads = OnePad()
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if log_x:
        pads[0].SetLogx(1)
    if custom_x_range:
        if x_axis_max > sighist.GetXaxis().GetXmax():
            x_axis_max = sighist.GetXaxis().GetXmax()
    if ratio:
        if log_x:
            pads[1].SetLogx(1)
        axish = createAxisHists(
            2, sighist, sighist.GetXaxis().GetXmin(), sighist.GetXaxis().GetXmax()
        )
        axish[1].GetXaxis().SetTitle(x_title)
        axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetXaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetNdivisions(4)
        axish[1].GetYaxis().SetTitle("Ratio")
        axish[1].GetYaxis().SetTitleOffset(1.6)
        axish[1].GetYaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetLabelSize(0.03)

        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    else:
        axish = createAxisHists(
            1,
            sighist,
            sighist.GetXaxis().GetXmin(),
            sighist.GetXaxis().GetXmax() - 0.01,
        )
        axish[0].GetXaxis().SetTitle(x_title)
        axish[0].GetXaxis().SetTitleSize(0.04)
        axish[0].GetXaxis().SetLabelSize(0.03)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)
    if not ratio:
        axish[0].GetXaxis().SetLabelSize(0.03)
    if not custom_y_range:
        if log_y:
            axish[0].SetMinimum(0.0009)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * sighist.GetMaximum()
                            - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            axish[0].SetMinimum(0)
            axish[0].SetMaximum(1.1 * (1 + extra_pad) * sighist.GetMaximum())
    axish[0].Draw()

    # Draw uncertainty band
    sighist.SetFillColor(CreateTransparentColor(12, 0.4))
    sighist.SetLineColor(CreateTransparentColor(12, 0.4))
    sighist.SetMarkerSize(0)
    sighist.SetMarkerColor(CreateTransparentColor(12, 0.4))

    error_hist = sighist.Clone()
    if do_custom_uncerts:
        bkg_uncert_up = infile.Get(nodename + "/" + custom_uncerts_up_name).Clone()
        bkg_uncert_down = infile.Get(nodename + "/" + custom_uncerts_down_name).Clone()
        if norm_bins:
            bkg_uncert_up.Scale(1.0, "width")
            bkg_uncert_down.Scale(1.0, "width")

        for i in range(1, bkg_uncert_up.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            bin_up = bkg_uncert_up.GetBinContent(i)
            bin_down = bkg_uncert_down.GetBinContent(i)
            error = abs(bin_up - bin_down) / 2
            band_center = max(bin_up, bin_down) - error
            if add_stat_to_syst:
                error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinContent(i, band_center)
            error_hist.SetBinError(i, error)

    if add_flat_uncert > 0:
        for i in range(1, error_hist.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            error = add_flat_uncert * error_hist.GetBinContent(i)
            error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinError(i, error)

    error_hist.Draw("e2same")
    stack.Draw("hist same")
    axish[0].Draw("axissame")

    # Setup legend
    legend = PositionedLegend(0.30, 0.1, 3, 0.03)
    legend.SetTextFont(42)
    legend.SetTextSize(0.022)
    legend.SetFillColor(0)
    if signal_mass != "":
        legend.AddEntry(sighist_copy, sig_schemes[signal_scheme][0], "l")
    if do_custom_uncerts and uncert_title != "":
        legend.AddEntry(error_hist, uncert_title, "f")
    else:
        legend.AddEntry(error_hist, "Uncertainty", "f")
    legend.Draw("same")
    if channel == "em":
        channel_label = "e#mu"
    if channel == "et":
        channel_label = "e#tau_{h}"
    if channel == "mt":
        channel_label = "#mu#tau_{h}"
    if channel == "tt":
        channel_label = "#tau_{h}#tau_{h}"
    if channel == "zmm":
        channel_label = "Z#rightarrow#mu#mu"
    if channel == "zee":
        channel_label = "Z#rightarrow ee"
    channel_label = ""
    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.145, 0.955, channel_label)

    # CMS and lumi labels
    if not custom_y_range:
        FixTopRange(pads[0], GetPadYMax(pads[0]), extra_pad if extra_pad > 0 else 0.30)
    DrawCMSLogo(pads[0], "CMS", "Simulation", 11, 0.045, 0.05, 1.0, "", 1.0)
    DrawCMSLogo(pads[0], "CMS", "", 11, 0.045, 0.05, 1.0, "", 1.0)
    DrawTitle(pads[0], lumi, 3)

    # Add ratio plot if required
    if ratio:
        ratio_sighist = MakeRatioHist(error_hist.Clone(), sighist.Clone(), True, False)
        sighist = MakeRatioHist(
            sighist_copy.Clone(), sighist.Clone(), True, False, ZeroFix=True
        )
        pads[1].cd()
        pads[1].SetGrid(0, 1)
        axish[1].Draw("axis")
        axish[1].SetMinimum(float(ratio_range.split(",")[0]))
        axish[1].SetMaximum(float(ratio_range.split(",")[1]))
        ratio_sighist.SetMarkerSize(0)
        ratio_sighist.Draw("e2same")
        sighist.DrawCopy("hist same")
        pads[1].RedrawAxis("G")

    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")


def TagAndProbePlot(
    graphs=[],
    legend_titles=[],
    title="",
    ratio=True,
    log_y=False,
    log_x=False,
    ratio_range="0.1,2.0",
    custom_x_range=False,
    x_axis_max=100,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=1.5,
    y_axis_min=0,
    x_title="",
    y_title="",
    extra_pad=0,
    plot_name="plot",
    label="",
    fits=None,
    ratio_fits=None,
):
    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    colourlist = [
        R.kBlue,
        R.kRed,
        R.kGreen + 3,
        R.kBlack,
        R.kYellow + 2,
        R.kOrange,
        R.kCyan + 3,
        R.kMagenta + 2,
        R.kViolet - 5,
        R.kGray,
    ]
    mg = R.TMultiGraph("mg", "")
    graph_count = 0
    legend_graphs = []
    hs = R.THStack("hs", "")
    for graph in graphs:
        g = graph.Clone()
        g.SetFillColor(0)
        g.SetLineWidth(3)
        g.SetLineColor(colourlist[graph_count])
        g.SetMarkerSize(0)
        mg.Add(g)
        hs.Add(g.GetHistogram())
        graph_count += 1
        legend_graphs.append(g.Clone())

    c1 = R.TCanvas()
    c1.cd()

    if ratio:
        pads = TwoPadSplit(0.29, 0.01, 0.01)
    else:
        pads = OnePad()
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if custom_x_range:
        if x_axis_max > graphs[0].GetHistogram().GetXaxis().GetXmax():
            x_axis_max = graphs[0].GetHistogram().GetXaxis().GetXmax()
    if ratio:
        axish = createAxisHists(
            2,
            graphs[0].GetHistogram(),
            graphs[0].GetHistogram().GetXaxis().GetXmin(),
            graphs[0].GetHistogram().GetXaxis().GetXmax() - 0.01,
        )
        axish[1].GetXaxis().SetTitle(x_title)
        axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetXaxis().SetTitle(x_title)
        axish[1].GetXaxis().SetTitleSize(0.04)
        axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetYaxis().SetNdivisions(4)
        axish[1].GetYaxis().SetTitle("SF")
        # axish[1].GetYaxis().SetTitle("Embed/MC")
        axish[1].GetYaxis().CenterTitle()
        axish[1].GetYaxis().SetTitleOffset(1.6)
        axish[1].GetYaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetLabelSize(0.03)

        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
            axish[1].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    else:
        axish = createAxisHists(
            1,
            graphs[0].GetHistogram(),
            graphs[0].GetHistogram().GetXaxis().GetXmin(),
            graphs[0].GetHistogram().GetXaxis().GetXmax() - 0.01,
        )
        axish[0].GetXaxis().SetLabelSize(0.03)
        axish[0].GetXaxis().SetTitle(x_title)
        axish[0].GetXaxis().SetTitleSize(0.04)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)

    if not custom_y_range:
        if log_y:
            axish[0].SetMinimum(0.0009)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * hs.GetMaximum("nostack")
                            - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            axish[0].SetMinimum(0)
            axish[0].SetMaximum(1.1 * (1 + extra_pad) * hs.GetMaximum("nostack"))
    axish[0].Draw()

    mg.Draw("p")

    mg.Draw("p")
    if fits is not None:
        fit_count = 0
        for i in fits:
            i.SetLineWidth(2)
            i.SetLineColor(colourlist[fit_count])
            i.Draw("same")
            fit_count += 1

    axish[0].Draw("axissame")

    # Setup legend
    legend = PositionedLegend(0.20, 0.1, 3, 0.01)
    legend.SetTextFont(42)
    legend.SetTextSize(0.022)
    legend.SetFillColor(0)

    for legi, graph in enumerate(legend_graphs):
        legend.AddEntry(graph, legend_titles[legi], "l")
    legend.Draw("same")

    # CMS label and title
    FixTopRange(pads[0], axish[0].GetMaximum(), extra_pad if extra_pad > 0 else 0.13)
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.05, 1.0, "", 1.0)
    DrawTitle(pads[0], title, 3)

    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.145, 0.955, label)

    # Add ratio plot if required
    if ratio:
        # ratio_hs = R.THStack("ratio_hs","")
        # hist_count=0
        pads[1].cd()
        pads[1].SetGrid(0, 1)
        axish[1].Draw("axis")
        axish[1].SetMinimum(float(ratio_range.split(",")[0]))
        axish[1].SetMaximum(float(ratio_range.split(",")[1]))
        div_hist = graphs[0].GetHistogram().Clone()

        num = graphs[0].Clone()
        ratio_graphs = []
        for i in range(1, len(graphs)):
            sf_graph = GraphDivide(num, graphs[i].Clone())
            sf_graph.SetLineWidth(3)
            sf_graph.SetLineColor(colourlist[i])
            sf_graph.SetMarkerSize(0)
            # f = R.TF1('f','pol1')
            # f.SetParameter(0,1.)
            # sf_graph.Fit('f')
            # sf_graph.GetFunction('f').SetLineColor(R.kBlack)
            # sf_graph.GetFunction('f').SetLineWidth(2)
            # R.gStyle.SetStatY(0.38)
            # R.gStyle.SetStatY(0.355)
            # R.gStyle.SetStatX(0.96)
            ratio_graphs.append(sf_graph.Clone())
        for x in ratio_graphs:
            x.Draw("p")

        if ratio_fits is not None:
            fit_count = 0
            for i in ratio_fits:
                i.SetLineWidth(2)
                # i.SetLineColor(colourlist[fit_count+1])
                i.SetLineColor(R.kBlack)
                i.Draw("same")
                fit_count += 1

        pads[1].RedrawAxis("G")
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")
    # del c1
    c1.Close()


def HTTPlotUnrolled(
    nodename,
    infile=None,
    signal_scale=1,
    signal_mass="",
    FF=False,
    norm_bins=True,
    channel="mt",
    blind=False,
    x_blind_min=0,
    x_blind_max=4000,
    auto_blind=False,
    ratio=True,
    log_y=False,
    log_x=False,
    ratio_range="0.7,1.3",
    custom_x_range=False,
    x_axis_max=4000,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=4000,
    y_axis_min=0,
    x_title="",
    y_title="Events/bin",
    extra_pad=0,
    do_custom_uncerts=False,
    add_stat_to_syst=False,
    add_flat_uncert=False,
    uncert_title="background uncertainty",
    lumi="35.9",
    plot_name="htt_plot",
    custom_uncerts_up_name="total_bkg_custom_uncerts_up",
    custom_uncerts_down_name="total_bkg_custom_uncerts_down",
    scheme="mt",
    cat="",
    x_lines=None,
    y_labels_vec=None,
    embedding=False,
    vbf_background=False,
    signal_scheme="",
):
    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    # Define signal schemes here
    sig_schemes = {}

    # sig_schemes['sm_ggH'] = ( str(int(signal_scale))+"#times SM ggH("+signal_mass+" GeV)#rightarrow#tau#tau", ["ggHsm_htt"], False , R.kRed)
    # sig_schemes['sm_qqH'] = ( str(int(signal_scale))+"#times SM qqH("+signal_mass+" GeV)#rightarrow#tau#tau", ["qqH_htt"], False, R.kBlue)

    sig_schemes["sm_cp"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggH_sm_htt"],
        False,
        R.kRed,
    )
    sig_schemes["ps_cp"] = (
        str(int(signal_scale)) + "#times SM ggH#rightarrow#tau#tau",
        ["ggH_ps_htt"],
        False,
        R.kRed,
    )
    sig_schemes["sm_cp_decays"] = (
        str(int(signal_scale)) + "#times SM H#rightarrow#tau#tau",
        ["ggH_sm_htt", "qqH_sm_htt"],
        False,
        R.kRed,
    )
    sig_schemes["sm_cp_decays_ps"] = (
        str(int(signal_scale)) + "#times PS H#rightarrow#tau#tau",
        ["ggH_ps_htt", "qqH_ps_htt"],
        False,
        R.kGreen + 3,
    )
    # sig_schemes['sm_ps'] = ( str(int(signal_scale))+"#times PS ggH#rightarrow#tau#tau", ["ggHps_htt"], False, R.kGreen+3)
    # sig_schemes['sm_mm'] = ( str(int(signal_scale))+"#times MM ggH#rightarrow#tau#tau", ["ggHmm_htt"], False, R.kOrange-5)

    ModTDRStyle(width=1200, height=600, r=0.3, l=0.14, t=0.12, b=0.15)
    R.TGaxis.SetExponentOffset(-0.06, 0.01, "y")

    background_schemes = {
        "mt": [
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp(
                "Z#rightarrow#mu#mu", ["ZL", "ZJ"], R.TColor.GetColor(100, 192, 232)
            ),
            backgroundComp(
                "Z#rightarrow#tau#tau",
                ["ZTT", "EWKZ"],
                R.TColor.GetColor(248, 206, 104),
            ),
            backgroundComp(
                "qqH#rightarrow#tau#tau + VH#rightarrow#tau#tau",
                ["qqH_htt125", "ZH_htt125", "WplusH_htt125", "WminusH_htt125"],
                R.TColor.GetColor(51, 51, 255),
            ),
        ],
        "et": [
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp(
                "Z#rightarrowee", ["ZL", "ZJ"], R.TColor.GetColor(100, 192, 232)
            ),
            backgroundComp(
                "Z#rightarrow#tau#tau",
                ["ZTT", "EWKZ"],
                R.TColor.GetColor(248, 206, 104),
            ),
            # backgroundComp("SM EWK H#rightarrow#tau#tau",["qqH_htt125","ZH_htt125", "WplusH_htt125","WminusH_htt125"],R.TColor.GetColor(51,51,255)),
        ],
        "tt": [
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "Electroweak",
                ["VVT", "VVJ", "W", "ZL", "ZJ"],
                R.TColor.GetColor(222, 90, 106),
            ),
            backgroundComp(
                "Z#rightarrow#tau#tau",
                ["ZTT", "EWKZ"],
                R.TColor.GetColor(248, 206, 104),
            ),
            # backgroundComp("SM EWK H#rightarrow#tau#tau",["qqH_htt125","ZH_htt125", "WplusH_htt125","WminusH_htt125"],R.TColor.GetColor(51,51,255)),
        ],
        "em": [
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "Electroweak", ["VVJ", "VVT", "W"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp("Z#rightarrowll", ["ZLL"], R.TColor.GetColor(100, 192, 232)),
            backgroundComp(
                "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
            ),
        ],
        "zm": [
            backgroundComp(
                "Misidentified #mu", ["QCD"], R.TColor.GetColor(250, 202, 255)
            ),
            backgroundComp("t#bar{t}", ["TT"], R.TColor.GetColor(155, 152, 204)),
            backgroundComp(
                "Electroweak", ["VV", "W", "ZJ"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp(
                "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
            ),
            backgroundComp(
                "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
            ),
        ],
        "zmm": [
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp(
                "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp(
                "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
            ),
            backgroundComp(
                "Z#rightarrow#mu#mu", ["ZL", "ZJ"], R.TColor.GetColor(100, 192, 232)
            ),
        ],
        "zee": [
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "t#bar{t}", ["TTT", "TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp(
                "Electroweak", ["VVT", "VVJ", "W"], R.TColor.GetColor(222, 90, 106)
            ),
            backgroundComp(
                "Z#rightarrow ee", ["ZL", "ZJ", "ZTT"], R.TColor.GetColor(100, 192, 232)
            ),
        ],
        "dy": [
            backgroundComp("Z#rightarrowll", ["ZLL"], R.TColor.GetColor(100, 192, 232)),
            backgroundComp(
                "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
            ),
        ],
        "w": [backgroundComp("W", ["W"], R.TColor.GetColor(222, 90, 106))],
        "w_shape": [
            backgroundComp(
                "W loosened shape", ["W_shape"], R.TColor.GetColor(222, 90, 106)
            )
        ],
        "qcd": [backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255))],
        "qcd_shape": [
            backgroundComp(
                "QCD loosened shape", ["QCD_shape"], R.TColor.GetColor(250, 202, 255)
            )
        ],
        "ff_comp": [
            backgroundComp(
                "t#bar{t} j#rightarrow#tau", ["TTJ"], R.TColor.GetColor(155, 152, 204)
            ),
            backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
            backgroundComp(
                "Electroweak j#rightarrow#tau",
                ["VVJ", "W", "ZJ"],
                R.TColor.GetColor(222, 90, 106),
            ),
        ],
    }
    if channel == "zee" or channel == "zmm":
        background_schemes["dy"] = [
            backgroundComp("DY", ["ZLL"], R.TColor.GetColor(100, 192, 232))
        ]
    if FF:
        background_schemes = {
            "mt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrow#mu#mu", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
            ],
            "et": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp("Electroweak", ["VVT"], R.TColor.GetColor(222, 90, 106)),
                backgroundComp(
                    "Z#rightarrowee", ["ZL"], R.TColor.GetColor(100, 192, 232)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Jet#rightarrow#tau_{h}",
                    ["jetFakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
            ],
            "tt": [
                backgroundComp("t#bar{t}", ["TTT"], R.TColor.GetColor(155, 152, 204)),
                backgroundComp(
                    "Electroweak", ["VVT", "ZL"], R.TColor.GetColor(222, 90, 106)
                ),
                backgroundComp(
                    "Z#rightarrow#tau#tau", ["ZTT"], R.TColor.GetColor(248, 206, 104)
                ),
                backgroundComp(
                    "Jet#rightarrow#tau_{h}",
                    ["jetFakes", "Wfakes"],
                    R.TColor.GetColor(192, 232, 100),
                ),
            ],
            "ff_comp": [
                backgroundComp(
                    "t#bar{t} jet#rightarrow#tau_{h}",
                    ["TTJ"],
                    R.TColor.GetColor(155, 152, 204),
                ),
                backgroundComp("QCD", ["QCD"], R.TColor.GetColor(250, 202, 255)),
                backgroundComp(
                    "Electroweak jet#rightarrow#tau_{h}",
                    ["VVJ", "W"],
                    R.TColor.GetColor(222, 90, 106),
                ),
                backgroundComp(
                    "Z#rightarrow ll jet#rightarrow#tau_{h}",
                    ["ZJ"],
                    R.TColor.GetColor(100, 192, 232),
                ),
            ],
        }

    if vbf_background:
        for key in background_schemes:
            background_schemes[key].insert(
                0,
                backgroundComp(
                    "qqH#rightarrow#tau#tau + VH#rightarrow#tau#tau",
                    ["qqH_htt125", "ZH_htt125", "WplusH_htt125", "WminusH_htt125"],
                    R.TColor.GetColor(51, 51, 230),
                ),
            )

    if embedding:
        for chan in ["em", "et", "mt", "tt", "zmm"]:
            if chan not in background_schemes:
                continue
            schemes = background_schemes[chan]
            for bkg in schemes:
                if chan != "zmm" and bkg["leg_text"] == "Z#rightarrow#tau#tau":
                    bkg["plot_list"] = ["EmbedZTT"]
                    bkg["leg_text"] = "#mu#rightarrow#tau embedding"
                if chan == "zmm" and bkg["leg_text"] == "Z#rightarrow#mu#mu":
                    bkg["plot_list"] = ["EmbedZL", "ZJ"]

    total_datahist = infile.Get(nodename + "/data_obs").Clone()
    if scheme == "w_shape":
        total_datahist = infile.Get(nodename + "/W").Clone()
    if scheme == "qcd_shape":
        total_datahist = infile.Get(nodename + "/QCD").Clone()
    if scheme == "ff_comp":
        total_datahist = infile.Get(nodename + "/jetFakes").Clone()

    # Create stacked plot for the backgrounds
    bkg_histos = []
    for i, t in enumerate(background_schemes[scheme]):
        plots = t["plot_list"]
        h = R.TH1F()
        for j, k in enumerate(plots):
            if not infile.Get(nodename + "/" + k):
                continue
            if h.GetEntries() == 0:
                h = infile.Get(nodename + "/" + k).Clone()

                h.SetName(k)
            else:
                h.Add(infile.Get(nodename + "/" + k).Clone())
        h.SetFillColor(t["colour"])
        h.SetLineColor(R.kBlack)
        h.SetMarkerSize(0)

        bkg_histos.append(h)

    stack = R.THStack("hs", "")
    bkghist_blind = R.TH1F()
    for hists in bkg_histos:
        if bkghist_blind.GetEntries() == 0:
            bkghist_blind = hists.Clone()
        else:
            bkghist_blind.Add(hists.Clone())

    bkghist = R.TH1F()

    for hists in bkg_histos:
        if norm_bins:
            Norm2DBins(hists)
        stack.Add(hists.Clone())
        if bkghist.GetEntries() == 0:
            bkghist = hists.Clone()
        else:
            bkghist.Add(hists.Clone())

    c1 = R.TCanvas()
    c1.cd()

    if ratio:
        pads = TwoPadSplit(0.29, 0.01, 0.01)
    else:
        pads = OnePad()
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if log_x:
        pads[0].SetLogx(1)

    if custom_x_range:
        if x_axis_max > bkghist.GetXaxis().GetXmax():
            x_axis_max = bkghist.GetXaxis().GetXmax()
    if ratio:
        if log_x:
            pads[1].SetLogx(1)
        axish = createAxisHists(
            2,
            bkghist,
            bkghist.GetXaxis().GetXmin(),
            bkghist.GetXaxis().GetXmax() - 0.01,
        )
        axish[1].GetXaxis().SetTitle(x_title)
        axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetXaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetNdivisions(4)
        if scheme == "w_shape" or scheme == "qcd_shape" or scheme == "ff_comp":
            axish[1].GetYaxis().SetTitle("Ratio")
        else:
            axish[1].GetYaxis().SetTitle("Obs/Exp")
        axish[1].GetYaxis().SetTitleOffset(0.8)
        axish[1].GetYaxis().SetTitleSize(0.04)
        axish[1].GetYaxis().SetLabelSize(0.03)
        axish[1].GetXaxis().SetTitleOffset(1.7)
        axish[1].GetXaxis().SetLabelOffset(0.015)
        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
            axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    else:
        axish = createAxisHists(
            1,
            bkghist,
            bkghist.GetXaxis().GetXmin(),
            bkghist.GetXaxis().GetXmax() - 0.01,
        )
        axish[0].GetXaxis().SetTitle(x_title)
        axish[0].GetXaxis().SetTitleSize(0.04)
        axish[0].GetXaxis().SetLabelSize(0.03)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(0.8)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)
    if not ratio:
        axish[0].GetXaxis().SetLabelSize(0.03)
    if not custom_y_range:
        if log_y:
            axish[0].SetMinimum(0.01)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * bkghist.GetMaximum()
                            - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            axish[0].SetMinimum(0)
            axish[0].SetMaximum(1.1 * (1 + extra_pad) * bkghist.GetMaximum())
    axish[0].Draw()

    # Draw uncertainty band
    bkghist.SetFillColor(CreateTransparentColor(12, 0.4))
    bkghist.SetLineColor(CreateTransparentColor(12, 0.4))
    bkghist.SetMarkerSize(0)
    bkghist.SetMarkerColor(CreateTransparentColor(12, 0.4))

    totsighist = R.TH1F()
    sighist = R.TH1F()
    sighists = []
    sighist_blind = []

    signal_split_schemes = []
    if signal_scheme == "cpprod_split":
        signal_split_schemes = ["sm_cp", "ps_sm"]

        print(sig_schemes["sm_cp"][1])
        print(sig_schemes["sm_cp"][1][0])
        h = infile.Get(nodename + "/" + sig_schemes["sm_cp"][1][0] + "125").Clone()
        h2 = infile.Get(nodename + "/" + sig_schemes["ps_cp"][1][0] + "125").Clone()

        h.SetLineColor(R.kRed)
        h.SetLineWidth(2)
        h.Scale(signal_scale)

        h2.SetLineColor(R.kGreen - 2)
        h2.SetLineWidth(2)
        h2.Scale(signal_scale)
        stack.Draw("histsame")
        h.Draw("histsame")
        h2.Draw("histsame")

    elif signal_mass != "":
        for signal_scheme in sig_schemes:
            sighist = R.TH1F()
            sig_scheme = sig_schemes[signal_scheme]
            for i in sig_scheme[1]:
                h = infile.Get(nodename + "/" + i + signal_mass).Clone()
                if sighist.GetEntries() == 0:
                    sighist = h
                else:
                    sighist.Add(h)
            sighist.SetLineColor(sig_scheme[3])
            sighist.SetLineWidth(2)
            sighist.Scale(signal_scale)
            if sig_scheme[2]:
                s = sighist.Clone()
                if norm_bins:
                    Norm2DBins(s)
                stack.Add(s)
                if not custom_y_range:
                    axish[0].SetMaximum(1.1 * (1 + extra_pad) * stack.GetMaximum())
            sighist_blind.append(sighist.Clone())
            if norm_bins:
                Norm2DBins(sighist)
            stack.Draw("histsame")
            sighist.SetName(signal_scheme)
            sighists.append(sighist.Clone())
            # if not sig_scheme[2]: sighist.Draw("histsame")
        for sig in sighists:
            sig.Draw("histsame")

    else:
        stack.Draw("histsame")

    totsighist = R.TH1F()
    for hists in sighist_blind:
        if totsighist.GetEntries() == 0:
            totsighist = hists.Clone()
        else:
            totsighist.Add(hists.Clone())

    # blinding histograms, autoblind option overwrites others

    blind_datahist = total_datahist.Clone()
    total_datahist.SetMarkerStyle(20)
    blind_datahist.SetMarkerStyle(20)
    blind_datahist.SetLineColor(1)

    # Auto blinding option
    if auto_blind:
        for i in range(1, total_datahist.GetNbinsX() + 1):
            b = bkghist_blind.GetBinContent(i)
            s = totsighist.GetBinContent(i)
            if PassAutoBlindMetric(s, b, metric=0.00001) or True:
                blind_datahist.SetBinContent(i, 0)
                blind_datahist.SetBinError(i, 0)
    # Blinding by hand using requested range, set to 200-4000 by default:
    elif blind:
        for i in range(0, total_datahist.GetNbinsX()):
            low_edge = total_datahist.GetBinLowEdge(i + 1)
            high_edge = low_edge + total_datahist.GetBinWidth(i + 1)
            if (low_edge > float(x_blind_min) and low_edge < float(x_blind_max)) or (
                high_edge > float(x_blind_min) and high_edge < float(x_blind_max)
            ):
                blind_datahist.SetBinContent(i + 1, 0)
                blind_datahist.SetBinError(i + 1, 0)
    if norm_bins:
        Norm2DBins(blind_datahist)
        Norm2DBins(total_datahist)

    error_hist = bkghist.Clone()
    if do_custom_uncerts:
        bkg_uncert_up = infile.Get(nodename + "/" + custom_uncerts_up_name).Clone()
        bkg_uncert_down = infile.Get(nodename + "/" + custom_uncerts_down_name).Clone()
        if norm_bins:
            Norm2DBins(bkg_uncert_up)
            Norm2DBins(bkg_uncert_down)

        for i in range(1, bkg_uncert_up.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            bin_up = bkg_uncert_up.GetBinContent(i)
            bin_down = bkg_uncert_down.GetBinContent(i)
            error = abs(bin_up - bin_down) / 2
            band_center = max(bin_up, bin_down) - error
            if add_stat_to_syst:
                error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinContent(i, band_center)
            error_hist.SetBinError(i, error)

    if add_flat_uncert > 0:
        for i in range(1, error_hist.GetNbinsX() + 1):
            stat_error = error_hist.GetBinError(i)
            error = add_flat_uncert * error_hist.GetBinContent(i)
            error = math.sqrt(error**2 + stat_error**2)
            error_hist.SetBinError(i, error)

    # if norm_bins:
    # Norm2DBins(error_hist)
    # Norm2DBins(blind_datahist)

    error_hist.Draw("e2same")
    blind_datahist.Draw("E same")
    axish[0].Draw("axissame")

    # Setup legend
    legend = PositionedLegend(0.13, 0.45, 7, 0.02)
    legend.SetTextFont(42)
    legend.SetTextSize(0.022)
    legend.SetFillColor(0)
    if scheme == "w_shape" or scheme == "qcd_shape":
        legend.AddEntry(blind_datahist, "un-loosened shape", "PE")
    elif scheme == "ff_comp":
        legend.AddEntry(blind_datahist, "FF jet#rightarrow#tau_{h}", "PE")
    else:
        legend.AddEntry(blind_datahist, "Observation", "PE")
    # Drawn on legend in reverse order looks better
    bkg_histos.reverse()
    background_schemes[scheme].reverse()
    for legi, hists in enumerate(bkg_histos):
        legend.AddEntry(hists, background_schemes[scheme][legi]["leg_text"], "f")
    if do_custom_uncerts and uncert_title != "":
        legend.AddEntry(error_hist, uncert_title, "f")
    else:
        legend.AddEntry(error_hist, "Background uncertainty", "f")
    if signal_mass != "":
        for sig in sighists:
            legend.AddEntry(sig, sig_schemes[sig.GetName()][0], "l")
    if scheme == "qcd_shape" or scheme == "w_shape":
        ks_score = error_hist.KolmogorovTest(total_datahist)
        legend.AddEntry(R.TObject(), "K-S probability = %.3f" % ks_score, "")
    legend.Draw("same")
    if channel == "em":
        channel_label = "e#mu"
    if channel == "et":
        channel_label = "e#tau_{h}"
    if channel == "mt":
        channel_label = "#mu#tau_{h}"
    if channel == "tt":
        channel_label = "#tau_{h}#tau_{h}"
    if channel == "zmm":
        channel_label = "Z#rightarrow#mu#mu"
    if channel == "zee":
        channel_label = "Z#rightarrow ee"
    if cat != "":
        channel_label += " " + cat
    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextAlign(23)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.46, 0.92, channel_label)

    # CMS and lumi labels
    if not custom_y_range:
        FixTopRange(pads[0], GetPadYMax(pads[0]), extra_pad if extra_pad > 0 else 0.30)
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.01, -0.16, 1.0, "", 0.4)
    DrawTitle(pads[0], lumi, 3, scale=0.5)

    # Add ratio plot if required
    # replace commented with uncommented/vice versa to plot ratio of signal hists
    # for comparison of ggH signal shapes
    if ratio:
        ratio_bkghist = MakeRatioHist(error_hist.Clone(), bkghist.Clone(), True, False)
        blind_ratio = MakeRatioHist(
            blind_datahist.Clone(), bkghist.Clone(), True, False
        )

        sighist_ratios = []
        # ks_scores = []
        # for sighist in sighists:
        #     print(sighist)
        #     sighist_ratio = MakeRatioHist(sighist.Clone(),sighists[0].Clone(),True,False)
        #     sighist_ratio.SetMarkerSize(0)
        # if sighist != sighists[0]:
        #     ks_score =  sighist.KolmogorovTest(sighists[0],"DNOU")
        #     legend.AddEntry(R.TObject(), "K-S probability = %.3f" % ks_score, "")
        # sighist_ratios.append(sighist_ratio)

        pads[1].cd()
        pads[1].SetGrid(0, 1)
        axish[1].Draw("axis")
        axish[1].SetMinimum(float(ratio_range.split(",")[0]))
        axish[1].SetMaximum(float(ratio_range.split(",")[1]))

        # pads[0].cd()
        # if len(sighist_ratios) >= 1:
        #     sighist_ratios[0].SetLineColor(R.kRed)
        #     sighist_ratios[0].DrawCopy("e0same")
        # if len(sighist_ratios) >= 2:
        #     sighist_ratios[1].SetLineColor(R.kGreen+3)
        #     sighist_ratios[1].DrawCopy("e0same")

        ratio_bkghist.SetMarkerSize(0)
        ratio_bkghist.Draw("e2same")
        blind_ratio.DrawCopy("e0same")
        pads[1].RedrawAxis("G")

    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    if x_lines is not None:
        line = R.TLine()
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.SetLineColor(R.kBlack)
        for x in x_lines:
            pads[0].cd()
            ymax = axish[0].GetMaximum()
            ymin = axish[0].GetMinimum()
            line.DrawLine(x, ymin, x, ymax)
            if ratio:
                pads[1].cd()
                ymax = axish[1].GetMaximum()
                ymin = axish[1].GetMinimum()
                line.DrawLine(x, ymin, x, ymax)

    if y_labels_vec is not None:
        pads[0].cd()
        unit = y_labels_vec[1][2]
        var = y_labels_vec[1][0]
        if "(" in var:
            var = var.split(" ")[0]
        y_bins = y_labels_vec[0]
        latex = R.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(R.kBlack)
        latex.SetTextSize(0.028)

        Nybins = len(y_bins)
        if Nybins > 4:
            latex.SetTextSize(0.023)
        if Nybins > 5:
            latex.SetTextSize(0.02)
        for i in range(0, Nybins):
            ymin = y_labels_vec[0][i][0]
            ymax = y_labels_vec[0][i][1]
            if var not in ["jeta_1", "MVA Score", "NN Score"]:
                if ymax == -1:
                    y_bin_label = "%s #geq %0.f %s" % (var, ymin, unit)
                else:
                    y_bin_label = "%0.f #leq %s < %0.f %s" % (ymin, var, ymax, unit)
            else:
                if ymax == -1:
                    y_bin_label = "%s #geq %.1f %s" % (var, ymin, unit)
                else:
                    y_bin_label = "%.1f #leq %s < %.1f %s" % (ymin, var, ymax, unit)
            if "tau_decay_mode" in var and Nybins == 3:
                if i == 0:
                    y_bin_label = "1 prong"
                if i == 1:
                    y_bin_label = "1 prong + #pi^{0}"
                if i == 2:
                    y_bin_label = "3 prong"
            xshift = (
                0.78 / Nybins * i
            )  # bit annoying but will have to change the 0.78 if the plot proportions are changed
            if Nybins > 5:
                xshift = (
                    0.76 / Nybins * i
                )  # bit annoying but will have to change the 0.78 if the plot proportions are changed
            latex.DrawLatex(0.095 + xshift, 0.82, y_bin_label)

    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")


def SoverBPlot(
    nodename="",
    infile=None,
    channel="",
    log_y=False,
    log_x=False,
    custom_x_range=False,
    x_axis_max=4000,
    x_axis_min=0,
    custom_y_range=False,
    y_axis_max=4000,
    y_axis_min=0,
    x_title="",
    extra_pad=0,
    plot_name="plot",
):

    R.gROOT.SetBatch(R.kTRUE)
    R.TH1.AddDirectory(False)
    ModTDRStyle(r=0.04, l=0.14)

    bkg_hist = infile.Get(nodename + "/total_bkg").Clone()
    signal_hist = infile.Get(nodename + "/ggHsm_htt125").Clone()
    signal_hist_2j = infile.Get(nodename + "/ggH2jsm_htt125").Clone()
    signal_hist.Add(signal_hist_2j)
    vbf_hist = infile.Get(nodename + "/qqH_htt125").Clone()

    def SumHist(h, sroot=False, bbb=False):
        for i in range(1, h.GetNbinsX() + 2):
            if bbb:
                last_bin = i
            else:
                last_bin = -1
            if sroot:
                h.SetBinContent(i, math.sqrt(h.Integral(i, last_bin)))
            else:
                h.SetBinContent(i, h.Integral(i, last_bin))

    SumHist(signal_hist, False, False)
    SumHist(bkg_hist, True, False)
    SumHist(vbf_hist, False, False)

    r_1 = signal_hist.Clone()
    r_2 = signal_hist.Clone()

    r_1.Divide(bkg_hist)
    r_2.Divide(vbf_hist)

    r_1.SetFillColor(0)
    r_1.SetLineWidth(3)
    r_1.SetLineColor(R.kRed)
    r_1.SetMarkerSize(0)

    r_2.SetFillColor(0)
    r_2.SetLineWidth(3)
    r_2.SetLineColor(R.kBlue)
    r_2.SetMarkerSize(0)

    c1 = R.TCanvas()
    c1.cd()

    pads = TwoPadSplit(0.29, 0.01, 0.01)
    pads[0].cd()

    if log_y:
        pads[0].SetLogy(1)
    if log_x:
        pads[0].SetLogx(1)
    if custom_x_range:
        if x_axis_max > r_2.GetXaxis().GetXmax():
            x_axis_max = r_2.GetXaxis().GetXmax()
    if log_x:
        pads[1].SetLogx(1)
    axish = createAxisHists(
        2, r_2, r_2.GetXaxis().GetXmin(), r_2.GetXaxis().GetXmax() - 0.01
    )
    axish[1].GetXaxis().SetTitle(x_title)
    axish[1].GetXaxis().SetLabelSize(0.03)
    axish[1].GetXaxis().SetTitleSize(0.04)
    axish[1].GetYaxis().SetNdivisions(4)
    axish[1].GetYaxis().SetTitle("ggH/qqH")
    axish[1].GetYaxis().SetTitleOffset(1.6)
    axish[1].GetYaxis().SetTitleSize(0.04)
    axish[1].GetYaxis().SetLabelSize(0.03)

    axish[0].GetXaxis().SetTitleSize(0)
    axish[0].GetXaxis().SetLabelSize(0)
    if custom_x_range:
        axish[0].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
        axish[1].GetXaxis().SetRangeUser(x_axis_min, x_axis_max - 0.01)
    if custom_y_range:
        axish[0].GetYaxis().SetRangeUser(y_axis_min, y_axis_max)
    axish[1].GetYaxis().SetRangeUser(r_2.GetMinimum() * 0.9, r_2.GetMaximum() * 1.1)

    axish[0].GetYaxis().SetTitle("S/#sqrt{B}")
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)

    if not custom_y_range:
        if log_y:
            axish[0].SetMinimum(0.0009)
            axish[0].SetMaximum(
                10
                ** (
                    (1 + extra_pad)
                    * (
                        math.log10(
                            1.1 * r_1.GetMaximum() - math.log10(axish[0].GetMinimum())
                        )
                    )
                )
            )
        else:
            axish[0].SetMinimum(0)
            axish[0].SetMaximum(1.1 * (1 + extra_pad) * r_1.GetMaximum())
    axish[0].SetLineWidth(0)
    axish[0].Draw()

    r_1.Draw("hist same")
    axish[0].Draw("axissame")

    # CMS label and title
    FixTopRange(pads[0], axish[0].GetMaximum(), extra_pad if extra_pad > 0 else 0.30)
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.05, 1.0, "", 1.0)

    if channel == "em":
        channel_label = "e#mu"
    if channel == "et":
        channel_label = "e#tau_{h}"
    if channel == "mt":
        channel_label = "#mu#tau_{h}"
    if channel == "tt":
        channel_label = "#tau_{h}#tau_{h}"
    if channel == "zmm":
        channel_label = "Z#rightarrow#mu#mu"
    if channel == "zee":
        channel_label = "Z#rightarrow ee"

    latex2 = R.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(R.kBlack)
    latex2.SetTextSize(0.028)
    latex2.DrawLatex(0.145, 0.955, channel_label)

    pads[1].cd()
    axish[1].Draw("axis")
    r_2.Draw("hist same")
    pads[1].RedrawAxis("G")

    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")
