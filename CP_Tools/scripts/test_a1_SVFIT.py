import ROOT
from CP_Tools.python.Utilities import RotateToGJMax, PV_Calculator
from CP_Tools.python.Particle import Particle

#170
tau = Particle(ROOT.TLorentzVector(-0.7832594513893127, -0.6216949224472046, 12.521081924438477, 12.68602480440152), -1)
vis1 = Particle(ROOT.TLorentzVector(4.279638809216447, 14.840565012333753, -8.351602088649077, 17.5592259355213), 1)
vis2 = Particle(ROOT.TLorentzVector(7.159955256868739, 23.969074551602173, -14.481874237859758, 28.905460649619663),-1)
vis3 = Particle(ROOT.TLorentzVector(4.8515708037304055, 14.721150115951904, -8.532049285788943, 17.69365266981444),-1)

tau_vis = vis1.Fourvector + vis2.Fourvector + vis3.Fourvector

x = RotateToGJMax(tau_vis, tau.Fourvector)

tau_new = Particle(x,-1)

y = PV_Calculator(tau_new,[vis1,vis2,vis3],"a1")
