from uncertainties import ufloat
import ROOT
from array import array

ROOT.TH1.SetDefaultSumw2(True)


def PrintSummary(
    analysis,
    nodename="",
    data_strings=["data_obs"],
    add_names="",
    channel="",
    samples_dict={},
):
    print("")
    print("################### Summary ###################")
    nodes = analysis.nodes[nodename].SubNodes()
    bkg_total = ufloat(0.000000001, 0.000000001)
    sig_total = ufloat(0.000000001, 0.000000001)
    for node in nodes:
        if channel == "em" and node.name == "W":
            continue
        if node.shape.rate.n == 0:
            per_err = 0
        else:
            per_err = node.shape.rate.s / node.shape.rate.n
        print(
            node.name.ljust(10),
            ("%.2f" % node.shape.rate.n).ljust(10),
            "+/-".ljust(5),
            ("%.2f" % node.shape.rate.s).ljust(7),
            "(%.4f)" % per_err,
        )
        if True in [
            node.name.find(add_name) != -1 and add_name != "" for add_name in add_names
        ]:
            continue
        if len(samples_dict["signal_samples"]) != 0:
            sig_samp_cond = [
                node.name.find(sig) != -1
                for sig in samples_dict["signal_samples"].keys()
            ]
        else:
            sig_samp_cond = []
        if True in sig_samp_cond:
            sig_total += node.shape.rate
        elif node.name not in data_strings:
            bkg_total += node.shape.rate
    if bkg_total.n == 0:
        per_err = 0
    else:
        per_err = bkg_total.s / bkg_total.n
    print(
        "Total bkg".ljust(10),
        ("%.2f" % bkg_total.n).ljust(10),
        "+/-".ljust(5),
        ("%.2f" % bkg_total.s).ljust(7),
        "(%.4f)" % per_err,
    )
    if sig_total.n == 0:
        per_err = 0
    else:
        per_err = sig_total.s / sig_total.n
    print(
        "Total sig".ljust(10),
        ("%.2f" % sig_total.n).ljust(10),
        "+/-".ljust(5),
        ("%.2f" % sig_total.s).ljust(7),
        "(%.4f)" % per_err,
    )
    print("###############################################")
    print("")


def GetTotals(ana, nodename="", add_name="", samples_dict={}, outfile="outfile.root"):
    # add histograms to get totals for backgrounds split into real/fake taus and make a total backgrounds histogram
    from itertools import chain

    outfile.cd(nodename)
    nodes = ana.nodes[nodename].SubNodes()
    first_hist = True
    for node in nodes:
        if add_name not in node.name:
            continue
        if node.name in list(chain.from_iterable(samples_dict["signal_samples"])):
            continue
        # if node.name == "data_obs": continue
        if node.name.endswith("Up"):
            continue
        if node.name.endswith("Down"):
            continue
        if node.name.endswith("extrap"):
            continue
        if first_hist:
            total_bkg = ana.nodes[nodename].nodes[node.name].shape.hist.Clone()
            first_hist = False
        else:
            total_bkg.Add(ana.nodes[nodename].nodes[node.name].shape.hist.Clone())
    if not first_hist:
        total_bkg.SetName("total_bkg" + add_name)
        total_bkg.Write()
    outfile.cd()


def FixBins(ana, nodename="", outfile="output.root"):
    # Fix empty histograms
    nodes = ana.nodes[nodename].SubNodes()
    for node in nodes:
        if "data_obs" in node.name:
            continue
        hist = outfile.Get(nodename + "/" + node.name)
        outfile.cd(nodename)
        # Fix empty histogram
        if hist.Integral() == 0.0:
            hist.SetBinContent(hist.GetNbinsX() // 2, 0.00001)
            hist.SetBinError(hist.GetNbinsX() // 2, 0.00001)
            hist.Write(hist.GetName(), ROOT.TObject.kWriteDelete)
        outfile.cd()


def Get1DBinNumFrom2D(h2d, xbin, ybin):
    Nxbins = h2d.GetNbinsX()
    return (ybin - 1) * Nxbins + xbin - 1


def Get1DBinNumFrom3D(h3d, xbin, ybin, zbin):
    Nxbins = h3d.GetNbinsX()
    Nybins = h3d.GetNbinsY()
    return (zbin - 1) * Nxbins * Nybins + (ybin - 1) * Nxbins + xbin - 1


def UnrollHist2D(h2d, inc_y_of=True):
    # inc_y_of = True includes the y over-flow bins
    if inc_y_of:
        n = 1
    else:
        n = 0
    Nbins = (h2d.GetNbinsY() + n) * (h2d.GetNbinsX())
    h1d = ROOT.TH1D(h2d.GetName(), "", Nbins, 0, Nbins)
    for i in range(1, h2d.GetNbinsX() + 1):
        for j in range(1, h2d.GetNbinsY() + 1 + n):
            glob_bin = Get1DBinNumFrom2D(h2d, i, j)
            content = h2d.GetBinContent(i, j)
            error = h2d.GetBinError(i, j)
            h1d.SetBinContent(glob_bin + 1, content)
            h1d.SetBinError(glob_bin + 1, error)
    return h1d


def UnrollHist3D(h3d, inc_y_of=False, inc_z_of=True):
    if inc_y_of:
        ny = 1
    else:
        ny = 0
    if inc_z_of:
        nz = 1
    else:
        nz = 0

    Nbins = (h3d.GetNbinsZ() + nz) * (h3d.GetNbinsY() + ny) * (h3d.GetNbinsX())
    h1d = ROOT.TH1D(h3d.GetName(), "", Nbins, 0, Nbins)
    for i in range(1, h3d.GetNbinsX() + 1):
        for j in range(1, h3d.GetNbinsY() + 1 + ny):
            for k in range(1, h3d.GetNbinsZ() + 1 + nz):
                glob_bin = Get1DBinNumFrom3D(h3d, i, j, k)
                content = h3d.GetBinContent(i, j, k)
                error = h3d.GetBinError(i, j, k)
                h1d.SetBinContent(glob_bin + 1, content)
                h1d.SetBinError(glob_bin + 1, error)
    return h1d


def FindRebinning(hist, BinThreshold=100, BinUncertFraction=0.5):
    # hist.Print("all")
    # getting binning
    binning = []
    for i in range(1, hist.GetNbinsX() + 2):
        binning.append(hist.GetBinLowEdge(i))

    # remove outer 0 bins
    still_zero = True
    for i in range(1, hist.GetNbinsX()):
        if not (hist.GetBinContent(i) == 0 and hist.GetBinError(i) == 0):
            still_zero = False
        if still_zero:
            binning.remove(min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i))))
        else:
            break
    hist = RebinHist(hist, binning)

    still_zero = True
    for i in reversed(range(2, hist.GetNbinsX() + 1)):
        if not (hist.GetBinContent(i) == 0 and hist.GetBinError(i) == 0):
            still_zero = False
        if still_zero:
            binning.remove(
                min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i + 1)))
            )
        else:
            break
    hist = RebinHist(hist, binning)

    # left to right
    finished = False
    k = 0
    while finished == False and k < 1000:
        k += 1
        for i in range(1, hist.GetNbinsX()):
            if hist.GetBinContent(i) > 0:
                uncert_frac = hist.GetBinError(i) / hist.GetBinContent(i)
            else:
                uncert_frac = BinUncertFraction + 1
            # print hist.GetBinLowEdge(i), uncert_frac, BinUncertFraction,  hist.GetBinContent(i), BinThreshold
            if uncert_frac > BinUncertFraction and hist.GetBinContent(i) < BinThreshold:
                # binning.remove(hist.GetBinLowEdge(i+1))
                binning.remove(
                    min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i + 1)))
                )
                hist = RebinHist(hist, binning)
                break
            elif i + 1 == hist.GetNbinsX():
                finished = True

    # right to left
    finished = False
    k = 0
    while finished == False and k < 1000:
        k += 1
        for i in reversed(range(2, hist.GetNbinsX() + 1)):
            if hist.GetBinContent(i) > 0:
                uncert_frac = hist.GetBinError(i) / hist.GetBinContent(i)
            else:
                uncert_frac = BinUncertFraction + 1
            # print hist.GetBinLowEdge(i), uncert_frac, BinUncertFraction,  hist.GetBinContent(i), BinThreshold
            if uncert_frac > BinUncertFraction and hist.GetBinContent(i) < BinThreshold:
                #        binning.remove(hist.GetBinLowEdge(i))
                binning.remove(
                    min(binning, key=lambda x: abs(x - hist.GetBinLowEdge(i)))
                )
                hist = RebinHist(hist, binning)
                break
            elif i == 2:
                finished = True

    #  hist.Print("all")
    return binning


def RebinHist(hist, binning):
    # getting initial binning
    initial_binning = []
    for i in range(1, hist.GetNbinsX() + 2):
        initial_binning.append(hist.GetBinLowEdge(i))

    new_binning = array("f", map(float, binning))
    hout = ROOT.TH1D(hist.GetName(), "", len(new_binning) - 1, new_binning)
    for i in range(1, hout.GetNbinsX() + 1):
        for j in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinCenter(j) > hout.GetBinLowEdge(i) and hist.GetBinCenter(
                j
            ) < hout.GetBinLowEdge(i + 1):
                new_content = hout.GetBinContent(i) + hist.GetBinContent(j)
                new_error = (hout.GetBinError(i) ** 2 + hist.GetBinError(j) ** 2) ** 0.5
                hout.SetBinContent(i, new_content)
                hout.SetBinError(i, new_error)
    # hout.Print("all")
    return hout
