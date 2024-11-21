import ROOT
import math
from scipy.stats import f as f_pval_func
import time
ROOT.gROOT.SetBatch(0)

f = ROOT.TFile("smoothing_2022_v2.root")
chan = 'pirho'
h_mc = f.Get('tt_higgs_%(chan)s/Higgs_flat_htt125' % vars())

h_mc.Rebin(4)
h_mc_clone = h_mc.Clone()

def GetFitStats(fit_result):
    chi2 = fit_result.Chi2()
    ndf = fit_result.Ndf()
    p_value = ROOT.TMath.Prob(chi2, ndf)

    return (chi2,ndf,p_value)

def FTest(chi2_1, ndf_1, chi2_2, ndf_2):
    # need to check this test is implemented correctly
    delta_chi2 = chi2_1 - chi2_2
    delta_ndf = ndf_1 - ndf_2
    F = (delta_chi2/delta_ndf) / (chi2_2/ndf_2)
    p_val = 1 - f_pval_func.cdf(F, delta_ndf, ndf_2)

    return(p_val)

def FindBestFuncs(data, N=5):
    funcs = []
    fits = []

    func0 = ROOT.TF1('func0','[0]',0,2*math.pi) # always start from a streight line fit
    fit = data.Fit(func0,'RS')
    funcs.append(func0)
    fits.append(fit)

    N_vec = list(range(1,N))

    func = func0
    while len(N_vec) > 0:
        min_chi2 = None
        for N in N_vec:
            func_str = str(func.GetExpFormula()).replace('p','')
            func_temp = ROOT.TF1('func%i' % func.GetNpar(), func_str+'+[%i]*cos(%i*x)' % (func.GetNpar(),N), 0, 2*math.pi)
            fit = data.Fit(func_temp,'RS')
            chi2, _, _ = GetFitStats(fit)
            if not min_chi2 or chi2 < min_chi2:
                min_chi2 = chi2
                N_best = N
                func_best = func_temp
                fit_best = fit
        funcs.append(func_best)
        fits.append(fit_best)
        func = func_best
        N_vec.remove(N_best)

    f_tests = [1.]
    for i in range(1,len(fits)):
        chi2_2, ndf_2, _ = GetFitStats(fits[i])
        chi2_1, ndf_1, _ = GetFitStats(fits[i-1])
        p_val = FTest(chi2_1, ndf_1, chi2_2, ndf_2)

    return (funcs)

#h_mc.Scale(1./h_mc.Integral())
funcs = FindBestFuncs(h_mc)

fout = ROOT.TFile('smoothing_func_output_%(chan)s_v2.root' % vars(),'RECREATE')
fout.cd()
h_mc_clone.Write('hist_%(chan)s' % vars())

colours = [2,4,6,1,3]

for i, f in enumerate(funcs):
    f.SetLineColor(colours[i])
    f.Write()

fout.Close()
