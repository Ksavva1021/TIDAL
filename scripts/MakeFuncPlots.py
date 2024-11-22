import ROOT
ROOT.gROOT.SetBatch(1)


chan = 'pirho'

file_name = 'smoothing_func_output_%(chan)s_v2.root' % vars()

f = ROOT.TFile(file_name)

h = f.Get('Higgs_flat_htt125_total' % vars())

func0 = f.Get('Higgs_flat_htt125_total_func0')
func1 = f.Get('Higgs_flat_htt125_total_func1')
func2 = f.Get('Higgs_flat_htt125_total_func2')
func3 = f.Get('Higgs_flat_htt125_total_func3')
func4 = f.Get('Higgs_flat_htt125_total_func4')


c1 = ROOT.TCanvas('c1','',800,600)
c1.SetRightMargin(0.3)

h.SetTitle('')
h.SetStats(0)
h.GetYaxis().SetTitle('Entries')
h.GetXaxis().SetTitle('#phi_{CP}')
h.Draw('')
func3.Draw('same')
func2.Draw('same')
func1.Draw('same')
func0.Draw('same')

legend = ROOT.TLegend(0.72, 0.5, 0.99, 0.95)

legend.AddEntry(func0, str(func0.GetExpFormula()).replace('p',''), 'l') 
legend.AddEntry(func1, str(func1.GetExpFormula()).replace('p',''), 'l') 
legend.AddEntry(func2, str(func2.GetExpFormula()).replace('p',''), 'l') 
legend.AddEntry(func3, str(func3.GetExpFormula()).replace('p',''), 'l') 
legend.SetBorderSize(0)
legend.Draw()

c1.Print('func_plot_%(chan)s.pdf' % vars())
h.SetMinimum(0)
c1.Print('func_plot_%(chan)s_full_y.pdf' % vars())
