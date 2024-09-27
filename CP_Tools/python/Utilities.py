import ROOT
import math
import vector
import numpy as np
from CP_Tools.python.SCalculator import SCalculator
import sys

def Boost(p,boost):
   '''
   Description: Function to boost a particle to a reference frame
   Arguments:
   1) p: TLorentzVector object [Four-vector of the particle to be boosted]
   2) boost: TLorentzVector.BoostVector() [Reference frame to which the particle is to be boosted]
   '''
   boost = vector.obj(px=boost.Px(),py=boost.Py(),pz=boost.Pz())
   p = vector.obj(px=p.Px(),py=p.Py(),pz=p.Pz(),E=p.E())
   p = p.boost(boost)
   return ROOT.TLorentzVector(p.px,p.py,p.pz,p.E)

def AcoplanarityAngle(p1,p2,p3,p4,type1,type2):
   '''
   Description: Function to calculate the Acoplanarity Angle distribution
   Arguments:
   1) p1: TLorentzVector object
   2) p2: TLorentzVector object
   3) p3: TLorentzVector object
   4) p4: TLorentzVector object

   For pi-pi: p1 = IP1, p2 = IP2, p3 = pi1, p4 = pi2
   For pi-rho: pi1 = IP1, pi2 = pi0_2, pi3 = pi1, pi4 = pi2
   For pi-a1: pi1 = IP1, pi2 = *, pi3 = pi1, pi4 = *
   For rho-a1: pi1 = pi0_1, pi2 = *, pi3 = pi1, pi4 = *

   * To replicate what is done for a1 resonances in the Run2 CP Analysis, use the sortA1() function to sort the hadrons. pi2 = first element of sortA1() function, pi4 = second element of sortA1() function

   5) type1: string [Type of the first resonance]
   6) type2: string [Type of the second resonance]
   '''
   boost = (p3+p4).BoostVector()

   # If the type contains "IP" (impact parameter), convert to unit vector
   if "IP" in type1: p1 = ROOT.TLorentzVector(p1.Vect().Unit(),0.)
   if "IP" in type2: p2 = ROOT.TLorentzVector(p2.Vect().Unit(),0.)

   p1 = Boost(p1,-boost)
   p2 = Boost(p2,-boost)
   p3 = Boost(p3,-boost)
   p4 = Boost(p4,-boost)

   # Calculate the plane normals
   n1 = p1.Vect() - p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit()
   n2 = p2.Vect() - p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit()

   n1 = n1.Unit()
   n2 = n2.Unit()

   # Calculate the acoplanarity angle
   angle = math.acos(n1.Dot(n2))
   sign = p4.Vect().Unit().Dot(n1.Cross(n2))

   if(sign<0): angle = 2*math.pi - angle

   cp_sign = 1
   # Calculating the phase shift
   # Note I think for A1 it must be p3,p1 and p4,p2 but, I did this to match what is done in the IC FW
   if "NP" in type1:
      cp_sign *= YRho(p3,p1,ROOT.TVector3())
   elif "A1" in type1:
      cp_sign *= YA1(p1,p3,ROOT.TVector3())
   if "NP" in type2:
      cp_sign *= YRho(p4,p2,ROOT.TVector3())
   elif "A1" in type2:
      cp_sign *= YA1(p2,p4,ROOT.TVector3())

   if "NP" in type1 or "NP" in type2 or "A1" in type1 or "A1" in type2:
      if(cp_sign<0):
         if(angle < math.pi):
            angle += math.pi
         else:
            angle -= math.pi
   return angle

def YRho(charged_pion,neutral_pion,boost):
  '''
  Description: Function to calculate phase shift for rho resonances
  charged_pion: TLorentzVector object [Four-vector of the charged pion]
  neutral_pion: TLorentzVector object [Four-vector of the neutral pion]
  boost: TLorentzVector.BoostVector() [Reference frame to which the particle is to be boosted]
  '''
  pi = Boost(charged_pion,-boost)
  pi0 = Boost(neutral_pion,-boost)

  E_pi = pi.E()
  E_pi0 = pi0.E()
  y = (E_pi - E_pi0)/(E_pi + E_pi0)
  return y

def YA1(pion,pion_OS,boost):
  '''
  Description: Function to calculate phase shift for a1 resonances
  pion: TLorentzVector object [Four-vector of the charged pion]
  pion_OS: TLorentzVector object [Four-vector of the oppositely charged pion]
  boost: TLorentzVector.BoostVector() [Reference frame to which the particle is to be boosted]
  '''
  pi = Boost(pion,-boost)
  pi_OS = Boost(pion_OS,-boost)

  E_pi = pi.E()
  E_pi_OS = pi_OS.E()
  y = (E_pi - E_pi_OS)/(E_pi + E_pi_OS)
  return y

def AcoplanarityAngle_PV(p1, p2, p3, p4):
   '''
   Description: Function to calculate the acoplanarity angle distribution using Polarimetric Vectors
   Arguments:
   1) p1: TLorentzVector object [Four-vector of PV1]
   2) p2: TLorentzVector object [Four-vector of PV2]
   3) p3: TLorentzVector object [Four-vector of the leading tau]
   4) p4: TLorentzVector object [Four-vector of the sub-leading tau]
   '''
   boost = (p3+p4).BoostVector()

   p1 = Boost(p1,-boost)
   p2 = Boost(p2,-boost)
   p3 = Boost(p3,-boost)
   p4 = Boost(p4,-boost)

   h1 = p1.Vect().Unit()
   h2 = p2.Vect().Unit()

   n1 = p3.Vect().Unit()
   n2 = p4.Vect().Unit()

   k1 = h1.Cross(n1).Unit()
   k2 = h2.Cross(n2).Unit()

   angle = math.acos(k1.Dot(k2))
   sign = h1.Cross(h2).Dot(n1)

   if(sign<0): angle = 2*math.pi - angle
   return angle

def sortA1(pi,pi1,pi2):
   '''
   Description: Function to re-arrange/sort the hadrons in the A1 resonance such that the 
   oppositely charged pion is located at the first element of the container.
   pi: TLorentzVector object [Four-vector of the first pion]
   pi1: TLorentzVector object [Four-vector of the second pion]
   pi2: TLorentzVector object [Four-vector of the third pion]
   '''
   hads = []
   # sort the hadrons so, the oppositely charged pion is at the first element
   if pi1.charge != pi.charge and pi1.charge != pi2.charge:
      hads.append(pi1)
      hads.append(pi)
      hads.append(pi2)

   elif pi2.charge != pi.charge and pi2.charge != pi1.charge:
      hads.append(pi2)
      hads.append(pi1)
      hads.append(pi)

   else:
      hads.append(pi)
      hads.append(pi1)
      hads.append(pi2)

   # from the two SS hadrons, pick the one closest to the rho mass as the second element
   rho_mass = 0.7755
   dM1 = abs((hads[0].Fourvector + hads[1].Fourvector).M() - rho_mass)
   dM2 = abs((hads[0].Fourvector + hads[2].Fourvector).M() - rho_mass)

   if (dM2 < dM1):
      hads[0], hads[1], hads[2] = hads[0], hads[2], hads[1]
   return hads

def AcoplanarityAngle_PV_A1A1(tau1, tau2, hads1, hads2):
   '''
   Description: Function to calculate the Acoplanarity Angle distribution using Polarimetric Vectors for A1 A1 (Full PV method)
   Arguments:
   1) tau1: TLorentzVector object [Four-vector of the leading tau]
   2) tau2: TLorentzVector object [Four-vector of the sub-leading tau]
   3) hads1: list of Particle objects [List of the hadrons in the leading tau]
   4) hads2: list of Particle objects [List of the hadrons in the sub-leading tau]
   '''
   angle = -9999.0
   a1_1 = ROOT.TLorentzVector(0,0,0,0)
   a1_2 = ROOT.TLorentzVector(0,0,0,0)
   for index in range(len(hads1)):
      a1_1 += hads1[index].Fourvector
   for index in range(len(hads2)):
      a1_2 += hads2[index].Fourvector

   charges_1 = [hads1[0].charge, hads1[1].charge, hads1[2].charge]
   charges_2 = [hads2[0].charge, hads2[1].charge, hads2[2].charge]

   pis_1 = [hads1[0].Fourvector, hads1[1].Fourvector, hads1[2].Fourvector]
   pis_2 = [hads2[0].Fourvector, hads2[1].Fourvector, hads2[2].Fourvector]

   Tauminus = tau1
   Tauplus = tau2

   Scalc = SCalculator("a1")
   if (Scalc.isOK("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2)):
      angle = Scalc.AcopAngle("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2)
   return angle

def PV_Calculator(tau1,hads1,type):
   '''
   Description: Function to calculate the Polarimetric Vectors
   Arguments:
   1) tau1: Particle object [Tau]
   2) hads1: list of Particle objects [List of the hadrons in the leading tau]
   3) type: string [Type of the resonance]
   '''
   vis = ROOT.TLorentzVector(0,0,0,0)

   for index in range(len(hads1)):
      vis += hads1[index].Fourvector

   if type == "a1":
      hads1 = sortA1(hads1[0],hads1[1],hads1[2])

   pions = [x.Fourvector for x in hads1]

   tauandprod = pions
   tauandprod.insert(0,tau1.Fourvector)

   SCalc = SCalculator(type)
   SCalc.Configure(tauandprod,(tau1.Fourvector),tau1.charge)
   pv1 = SCalc.pv()
   pv1 = ROOT.TLorentzVector(pv1.Px(),pv1.Py(),pv1.Pz(),0)
   pv_final = Boost(pv1,(tau1.Fourvector).BoostVector())
   return pv_final

# ------------------------------------------------------------------------------------------------
# Run 2 CP Tools for PV A1 A1

def AcoplanarityAngle_A1A1(tau1, tau2, hads1, hads2):
   '''
   tau1: TVector3 [Leading Tau Direction]
   tau2: TVector3 [Sub-leading Tau Direction]
   hads1: list of Particle objects [List of the hadrons in the leading tau]
   hads2: list of Particle objects [List of the hadrons in the sub-leading tau]

   In the Run 2 CP-Analysis, tau1 and tau2 were set to SV-PV as a proxy of the tau direction.
   '''
   angle = -9999.0
   a1_1 = ROOT.TLorentzVector(0,0,0,0)
   a1_2 = ROOT.TLorentzVector(0,0,0,0)
   for index in range(len(hads1)):
      a1_1 += hads1[index].Fourvector
   for index in range(len(hads2)):
      a1_2 += hads2[index].Fourvector

   charges_1 = [hads1[0].charge, hads1[1].charge, hads1[2].charge]
   charges_2 = [hads2[0].charge, hads2[1].charge, hads2[2].charge]
   pis_1 = [hads1[0].Fourvector, hads1[1].Fourvector, hads1[2].Fourvector]
   pis_2 = [hads2[0].Fourvector, hads2[1].Fourvector, hads2[2].Fourvector]

   dummy1 = ROOT.TLorentzVector()
   dummy2 = ROOT.TLorentzVector()
   solutions = tauPairMomentumSolutions(tau1, a1_1, dummy1, True, tau2, a1_2, dummy2, True, True)

   Tauminus = solutions[3]
   Tauplus = solutions[7]

   Scalc = SCalculator("a1")
   if (Scalc.isOK("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2)):
      angle = Scalc.AcopAngle("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2)
   return angle

def tauPairMomentumSolutions(tau1Dir, a1RefitLV1, a1LV1, is1Real, tau2Dir, a1RefitLV2, a1LV2, is2Real, isRefit):

   if (isRefit):
      tau1PairConstraintLV, tau2PairConstraintLV = tauPairConstraint(tau1Dir, a1RefitLV1, is1Real, tau2Dir, a1RefitLV2, is2Real)
      tau1Solutions = tauMomentumSolutions(tau1Dir, a1RefitLV1, is1Real)
      tau2Solutions = tauMomentumSolutions(tau2Dir, a1RefitLV2, is2Real)
   else:
      tau1PairConstraintLV, tau2PairConstraintLV = tauPairConstraint(tau1Dir, a1LV1, is1Real, tau2Dir, a1LV2, is2Real)
      tau1Solutions = tauMomentumSolutions(tau1Dir, a1LV1, is1Real)
      tau2Solutions = tauMomentumSolutions(tau2Dir, a1LV2, is2Real)
   solutions = tau1Solutions
   solutions.append(tau1PairConstraintLV)
   solutions += tau2Solutions
   solutions.append(tau2PairConstraintLV)
   return solutions

def tauPairConstraint(tau1Dir, a1LV1, is1Real,tau2Dir, a1LV2, is2Real):
   Hmass = 125.10

   tau1Solutions = tauMomentumSolutions(tau1Dir, a1LV1, is1Real)
   tau2Solutions = tauMomentumSolutions(tau2Dir, a1LV2, is2Real)

   massPairs = []
   massPairs.append((tau1Solutions[0] + tau2Solutions[0]).M()) # tau1:small solution -- tau2:small solution
   massPairs.append((tau1Solutions[0] + tau2Solutions[1]).M()) # tau1:small solution -- tau2:large solution
   massPairs.append((tau1Solutions[1] + tau2Solutions[0]).M()) # tau1:large solution -- tau2:small solution
   massPairs.append((tau1Solutions[1] + tau2Solutions[1]).M()) # tau1:large solution -- tau2:large solution

   massConstraint = []
   massConstraint.append(abs(massPairs[0] - Hmass))
   massConstraint.append(abs(massPairs[1] - Hmass))
   massConstraint.append(abs(massPairs[2] - Hmass))
   massConstraint.append(abs(massPairs[3] - Hmass))

   bestMassConstraint = massConstraint[0]
   bestCouple = 0

   for index in range(4):
      if (massConstraint[index] < bestMassConstraint):
         bestMassConstraint = massConstraint[index]
         bestCouple = index

   if bestCouple == 0:
      tau1PairConstraintLV = tau1Solutions[0]
      tau2PairConstraintLV = tau2Solutions[0]
   elif bestCouple == 1:
      tau1PairConstraintLV = tau1Solutions[0]
      tau2PairConstraintLV = tau2Solutions[1]
   elif bestCouple == 2:
      tau1PairConstraintLV = tau1Solutions[1]
      tau2PairConstraintLV = tau2Solutions[0]
   elif bestCouple == 3:
      tau1PairConstraintLV = tau1Solutions[1]
      tau2PairConstraintLV = tau2Solutions[1]
   return tau1PairConstraintLV, tau2PairConstraintLV

def tauMomentumSolutions(tauDir, a1LV, isReal):

   tauMass = 1.77682
   thetaGJ = a1LV.Angle(tauDir)
   tauDirUnit = tauDir.Unit()

   a = 4.0 * (a1LV.M2() + a1LV.P()**2 * np.sin(thetaGJ)**2)
   b = -4.0 * (a1LV.M2() + tauMass**2) * a1LV.P() * np.cos(thetaGJ)
   c = 4.0 * tauMass**2 * (a1LV.M2() + a1LV.P()**2) - (a1LV.M2() + tauMass**2)**2

   tauMomentumSmall, tauMomentumLarge = quadratic_alternate(a,b,c,isReal)
   tauMomentumMean = (tauMomentumSmall + tauMomentumLarge)/2.0

   tauSmall = ROOT.TLorentzVector(tauMomentumSmall * tauDirUnit.x(), tauMomentumSmall * tauDirUnit.y(), tauMomentumSmall * tauDirUnit.z(), np.sqrt(tauMomentumSmall**2 + tauMass**2))
   tauLarge = ROOT.TLorentzVector(tauMomentumLarge * tauDirUnit.x(), tauMomentumLarge * tauDirUnit.y(), tauMomentumLarge * tauDirUnit.z(), np.sqrt(tauMomentumLarge**2 + tauMass**2))
   tauMean = ROOT.TLorentzVector(tauMomentumMean * tauDirUnit.x(), tauMomentumMean * tauDirUnit.y(), tauMomentumMean * tauDirUnit.z(), np.sqrt(tauMomentumMean**2 + tauMass**2))

   solutions = [tauSmall, tauLarge, tauMean]
   return solutions

def quadratic_alternate(a,b,c,isReal):
   D = b**2 - 4.0*a*c
   if D >= 0:
      isReal = True
      q = -0.5 * (b + math.copysign(np.sqrt(D),b))
      x_0 = c/q
      x_1 = q/a
   else:
      isReal = False
      q = -0.5 * b
      x_0 = c/q
      x_1 = q/a
   return x_0, x_1


def RotateToGJMax(vis_tau_vec,  tau_vec):
   tau_dir = tau_vec.Vect().Unit()
   vis_dir = vis_tau_vec.Vect().Unit()
   tau_vis = vis_tau_vec

   m_tau = tau_vec.M()
   m_vis = vis_tau_vec.M()
   theta_GJ = math.acos(np.clip(tau_dir.Dot(vis_dir.Unit()), -1, 1))
   theta_GJ_max = math.asin(np.clip((m_tau**2 - m_vis**2)/(2*m_tau*tau_vis.P()), -1, 1))

   if (theta_GJ > theta_GJ_max):
      theta_GJ = theta_GJ_max

      n_1_x = 1 / math.sqrt(1+(tau_vis.X()/tau_vis.Y())**2)
      n_1_y = -n_1_x * tau_vis.X()/tau_vis.Y()
      n_1 = ROOT.TVector3(n_1_x, n_1_y, 0)
      n_2 = n_1.Cross(tau_vis.Vect().Unit())

      phi_opt_1 = math.atan(tau_dir.Dot(n_2)/tau_dir.Dot(n_1))
      new_dir_1 = math.cos(theta_GJ)*vis_dir + math.sin(theta_GJ)*(math.cos(phi_opt_1)*n_1 + math.sin(phi_opt_1)*n_2)

      phi_opt_2 = phi_opt_1 + math.pi
      new_dir_2 = math.cos(theta_GJ)*vis_dir + math.sin(theta_GJ)*(math.cos(phi_opt_2)*n_1 + math.sin(phi_opt_2)*n_2)

      if (new_dir_1.Dot(tau_dir) > new_dir_2.Dot(tau_dir)):
         new_dir = new_dir_1
      else:
         new_dir = new_dir_2

      new_pt = tau_vec.P()*math.sin(math.atan(math.exp(-new_dir.Eta()))*2)
      new_tau_vec = ROOT.Math.PtEtaPhiEVector(new_pt, new_dir.Eta(), new_dir.Phi(), tau_vec.E())
      new_tau_vec = ROOT.TLorentzVector(new_tau_vec.Px(), new_tau_vec.Py(), new_tau_vec.Pz(), new_tau_vec.E())
      return new_tau_vec
   else:
      return ROOT.TLorentzVector(tau_vec.Px(), tau_vec.Py(), tau_vec.Pz(), tau_vec.E())