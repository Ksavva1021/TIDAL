import numpy as np
import ROOT
import math
from CP_Tools.python.PolarimetricA1 import PolarimetricA1

class SCalculator:

    def __init__(self,type):
        self.TauAndProd_HRF_ = []
        self.charge_ = 0
        self.type_ = type

    def Configure(self, TauAndProd, Frame, charge):
        for index in range(len(TauAndProd)):
            self.TauAndProd_HRF_.append(self.Boost(TauAndProd[index],Frame))
        self.charge_ = charge

    def convertToMatrix(self, V):
        M = np.empty((len(V),1), dtype=float)
        for index in range(len(V)):
            M[index][0] = V[index]
        return M

    def Boost(self, pB, frame):
        if frame.Vect().Mag() == 0:
            print("RH Boost is not set, perfrom calculation in the Lab Frame")
            return pB
        if frame.E() == 0:
            print("Caution: Please check that you perform boost correctly!")
            return pB
        else:
            b = frame.Vect() * (1/frame.E())

        vec = [pB.E(), pB.Px(), pB.Py(), pB.Pz()]
        gamma = 1 / np.sqrt(1 - b.Mag2())

        transform = np.empty((4, 4), dtype=float)

        transform[0][0] = gamma
        transform[0][1] = -gamma * b.X()
        transform[0][2] = -gamma * b.Y()
        transform[0][3] = -gamma * b.Z()

        transform[1][0] = -gamma * b.X()
        transform[1][1] = 1 + (gamma - 1) * b.X() * b.X() / b.Mag2()
        transform[1][2] = (gamma - 1) * b.X() * b.Y() / b.Mag2()
        transform[1][3] = (gamma - 1) * b.X() * b.Z() / b.Mag2()

        transform[2][0] = -gamma * b.Y()
        transform[2][1] = (gamma - 1) * b.Y() * b.X() / b.Mag2()
        transform[2][2] = 1 + (gamma - 1) * b.Y() * b.Y() / b.Mag2()
        transform[2][3] = (gamma - 1) * b.Y() * b.Z() / b.Mag2()

        transform[3][0] = -gamma * b.Z()
        transform[3][1] = (gamma - 1) * b.Z() * b.X() / b.Mag2()
        transform[3][2] = (gamma - 1) * b.Z() * b.Y() / b.Mag2()
        transform[3][3] = 1 + (gamma - 1) * b.Z() * b.Z() / b.Mag2()

        result = transform @ self.convertToMatrix(vec)
        return ROOT.TLorentzVector(result[1][0],result[2][0],result[3][0],result[0][0])

    def SortPions(self, pionsvec, charges):
        npim = 0
        npip = 0
        npin = 0

        OSCharge = -99
        OSMCPionIndex = -99

        SS1Charge = -99
        SSMCPion1Index = -99

        SS2Charge = -99
        SSMCPion2Index = -99

        MCNeutralPionIndex = -99
        MCChargedPionIndex = -99
        NeutralPionCharge = -99
        ChargedPionCharge = -99

        if len(charges) == 3:
            for index in range(3):
                if charges[index] == 1:
                    npip += 1
                if charges[index] == -1:
                    npim += 1

            if npip == 1 and npim == 2:
                nss = 0
                for index in range(3):
                    if charges[index] == 1:
                        OSCharge = 1
                        OSMCPionIndex = index
                    if charges[index] == -1 and nss == 0:
                        nss += 1
                        SS1Charge = -1
                        SSMCPion1Index = index
                    if charges[index] == -1 and nss == 1:
                        SS2Charge = -1
                        SSMCPion2Index = index

            elif npip == 2 and npim == 1:
                nss = 0
                for index in range(3):
                    if charges[index] == -1:
                        OSCharge = -1
                        OSMCPionIndex = index
                    if charges[index] == 1 and nss == 0:
                        nss += 1
                        SS1Charge = 1
                        SSMCPion1Index = index
                    if charges[index] == 1 and nss == 1:
                        SS2Charge = 1
                        SSMCPion2Index = index

            if OSMCPionIndex != -99 and SSMCPion1Index != -99 and SSMCPion2Index != -99:
                os = pionsvec[OSMCPionIndex]
                ss1 = pionsvec[SSMCPion1Index]
                ss2 = pionsvec[SSMCPion2Index]

                charges.clear()
                charges.append(OSCharge)
                charges.append(SS1Charge)
                charges.append(SS2Charge)

                pionsvec.clear()
                pionsvec.append(os)
                pionsvec.append(ss1)
                pionsvec.append(ss2)

            else:
                charges.clear()
                charges.append(-99)
                charges.append(-99)
                charges.append(-99)

        if len(charges) == 2:
            for index in range(2):
                if charges[index] == 1: npip += 1
                if charges[index] == -1: npim += 1
                if charges[index] == 0: npin += 1

            if npip == 1 and npin == 1:
                for index in range(2):
                    if charges[index] == 1:
                        ChargedPionCharge = 1
                        MCChargedPionIndex = index
                    if charges[index] == 0:
                        NeutralPionCharge = 0
                        MCNeutralPionIndex = index

            if npim == 1 and npin == 1:
                for index in range(2):
                    if charges[index] == -1:
                        ChargedPionCharge = -1
                        MCChargedPionIndex = index
                    if charges[index] == 0:
                        NeutralPionCharge = 0
                        MCNeutralPionIndex = index

            ChargedPion = pionsvec[MCChargedPionIndex]
            NeutralPion = pionsvec[MCNeutralPionIndex]

            charges.clear()
            charges.append(ChargedPionCharge)
            charges.append(NeutralPionCharge)

            pionsvec.clear()
            pionsvec.append(ChargedPion)
            pionsvec.append(NeutralPion)

    def pv(self):
        if self.type_ == "pion":
            out = self.TauAndProd_HRF_[1].Vect()

        elif self.type_ == "rho":
            pi = self.TauAndProd_HRF_[1]
            pi0 = self.TauAndProd_HRF_[2]
            Tau = self.TauAndProd_HRF_[0]

            q = pi - pi0
            P = Tau
            N = Tau - pi - pi0

            out = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)))

        elif self.type_ == "a1":
            a1pol = PolarimetricA1()
            a1pol.Configure(self.TauAndProd_HRF_, self.charge_)
            out = -a1pol.PVC().Vect()

        return out

    def isOK(self,type1,type2,tauMinus,sumPionsMinus,sumPionsChargeMinus,tauPlus,sumPionsPlus,sumPionsChargePlus):

        zeroLV = ROOT.TLorentzVector(0,0,0,0)

        Scalc1 = SCalculator(type1)
        Scalc2 = SCalculator(type2)

        if type1 != "pion": Scalc1.SortPions(sumPionsMinus,sumPionsChargeMinus)
        if type2 != "pion": Scalc2.SortPions(sumPionsPlus,sumPionsChargePlus)

        isgoodcharge = True
        if (sumPionsChargeMinus[0] == -99 or sumPionsChargeMinus[1] == -99 or sumPionsChargeMinus[2] == -99 or sumPionsChargePlus[0] == -99 or sumPionsChargePlus[1] == -99 or sumPionsChargePlus[2] == -99):
            isgoodcharge = False

        tauandprodminus = []
        tauandprodplus = []
        pionszero = False

        tauandprodminus.append(tauMinus)
        for index in range(len(sumPionsMinus)):
            tauandprodminus.append(sumPionsMinus[index])
            if sumPionsMinus[index] == zeroLV:
                pionszero = True

        tauandprodplus.append(tauPlus)
        for index in range(len(sumPionsPlus)):
            tauandprodplus.append(sumPionsPlus[index])
            if sumPionsPlus[index] == zeroLV:
                pionszero = True

        Scalc1.Configure(tauandprodminus, tauandprodminus[0] + tauandprodplus[0], -1)
        h1 = Scalc1.pv()

        Scalc2.Configure(tauandprodplus, tauandprodminus[0] + tauandprodplus[0], +1) 
        h2 = Scalc2.pv()

        if math.isnan(h1.Mag()) or math.isnan(h2.Mag()) or tauMinus == zeroLV or tauPlus == zeroLV or tauMinus == tauPlus or pionszero or sumPionsPlus == sumPionsMinus or not isgoodcharge:
            return False
        else:
            return True

    def AcopAngle(self, type1, type2, tauMinus, sumPionsMinus, sumPionsChargeMinus, tauPlus, sumPionsPlus, sumPionsChargePlus):

        Scalc1 = SCalculator(type1)
        Scalc2 = SCalculator(type2)

        if type1 != "pion": Scalc1.SortPions(sumPionsMinus,sumPionsChargeMinus)
        if type2 != "pion": Scalc2.SortPions(sumPionsPlus,sumPionsChargePlus)

        tauandprodminus = []
        tauandprodplus = []

        tauandprodminus.append(tauMinus)
        for index in range(len(sumPionsMinus)):
            tauandprodminus.append(sumPionsMinus[index])

        tauandprodplus.append(tauPlus)
        for index in range(len(sumPionsPlus)):
            tauandprodplus.append(sumPionsPlus[index])

        Scalc1.Configure(tauandprodminus, tauandprodminus[0] + tauandprodplus[0], -1)
        h1 = Scalc1.pv()

        Scalc2.Configure(tauandprodplus, tauandprodminus[0] + tauandprodplus[0], +1)
        h2 = Scalc2.pv()

        tauMinus_HRF = Scalc1.Boost(tauandprodminus[0], tauandprodminus[0] + tauandprodplus[0])
        tauPlus_HRF = Scalc2.Boost(tauandprodplus[0], tauandprodminus[0] + tauandprodplus[0])

        h1Norm = 1.0 / h1.Mag()
        h2Norm = 1.0 / h2.Mag()
        h1 = h1 * h1Norm
        h2 = h2 * h2Norm
        k1Norm=1./((h1.Cross(tauMinus_HRF.Vect().Unit())).Mag())
        k2Norm=1./((h2.Cross(tauPlus_HRF.Vect().Unit())).Mag())
        k1 = (h1.Cross(tauMinus_HRF.Vect().Unit()))*k1Norm
        k2 = (h2.Cross(tauPlus_HRF.Vect().Unit()))*k2Norm

        if(((h1.Cross(h2))*(tauMinus_HRF.Vect().Unit()))<=0):
            return math.atan2((k1.Cross(k2)).Mag(),k1*k2)
        else:
            return (2.*math.pi-math.atan2((k1.Cross(k2)).Mag(),k1*k2))