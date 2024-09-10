from CP_Tools.python.Particle import Particle
import ROOT
import uproot3
import pandas as pd
from CP_Tools.python.Utilities import AcoplanarityAngle, sortA1, AcoplanarityAngle_A1A1
from alive_progress import alive_bar

filepath = "/vols/cms/ks1021/TIDAL/CP_Tools/samples/merged.root"
file = uproot3.open(filepath)

variables = [
    'pi_pt_1','pi_eta_1','pi_phi_1','pi_mass_1',
    'pi_pt_2','pi_eta_2','pi_phi_2','pi_mass_2',
    'pi2_pt_1','pi2_eta_1','pi2_phi_1','pi2_mass_1',
    'pi2_pt_2','pi2_eta_2','pi2_phi_2','pi2_mass_2',
    'pi3_pt_1','pi3_eta_1','pi3_phi_1','pi3_mass_1',
    'pi3_pt_2','pi3_eta_2','pi3_phi_2','pi3_mass_2',
    'pi0_pt_1','pi0_eta_1','pi0_phi_1','pi0_mass_1',
    'pi0_pt_2','pi0_eta_2','pi0_phi_2','pi0_mass_2',
    'ip_x_1','ip_y_1','ip_z_1',
    'ip_x_2','ip_y_2','ip_z_2',
    'sv_x_1','sv_y_1','sv_z_1',
    'sv_x_2','sv_y_2','sv_z_2',
    'pv_x','pv_y','pv_z',
    'decayMode_1','decayMode_2',
    'pi_charge_1','pi_charge_2',
    'pi2_charge_1','pi2_charge_2',
    'pi3_charge_1','pi3_charge_2',
    'pi0_charge_1','pi0_charge_2',
    'charge_1','charge_2',
    'event','os',
]

df = file['ntuple'].pandas.df(variables)

dict = {
    'aco_angle_0_0': [],
    'aco_angle_0_1': [],
    'aco_angle_1_0': [],
    'aco_angle_1_1': [],
    'aco_angle_0_10': [],
    'aco_angle_10_0': [],
    'aco_angle_1_10': [],
    'aco_angle_10_1': [],
    'aco_angle_10_10': []
}

# set up alive bar tracking
with alive_bar(df.shape[0]) as bar:
    for index, row in df.iterrows():
        ip1 = ROOT.TLorentzVector(row['ip_x_1'], row['ip_y_1'], row['ip_z_1'], 0)

        ip2 = ROOT.TLorentzVector(row['ip_x_2'], row['ip_y_2'], row['ip_z_2'], 0)

        sv1 = ROOT.TLorentzVector(row['sv_x_1'], row['sv_y_1'], row['sv_z_1'], 0)

        sv2 = ROOT.TLorentzVector(row['sv_x_2'], row['sv_y_2'], row['sv_z_2'], 0)

        pv = ROOT.TLorentzVector(row['pv_x'], row['pv_y'], row['pv_z'], 0)

        sv1 = ROOT.TLorentzVector(sv1.Vect()-pv.Vect(),0)
        sv2 = ROOT.TLorentzVector(sv2.Vect()-pv.Vect(),0)

        vector_pi_1 = ROOT.TLorentzVector()
        vector_pi_1.SetPtEtaPhiM(row['pi_pt_1'], row['pi_eta_1'], row['pi_phi_1'], row['pi_mass_1'])
        cartesian_vector_pi_1 = ROOT.TLorentzVector(vector_pi_1.Px(), vector_pi_1.Py(), vector_pi_1.Pz(), vector_pi_1.E())
        pi_1 = Particle(cartesian_vector_pi_1, row['pi_charge_1'])

        vector_pi_2 = ROOT.TLorentzVector()
        vector_pi_2.SetPtEtaPhiM(row['pi_pt_2'], row['pi_eta_2'], row['pi_phi_2'], row['pi_mass_2'])
        cartesian_vector_pi_2 = ROOT.TLorentzVector(vector_pi_2.Px(), vector_pi_2.Py(), vector_pi_2.Pz(), vector_pi_2.E())
        pi_2 = Particle(cartesian_vector_pi_2, row['pi_charge_2'])

        vector_pi2_1 = ROOT.TLorentzVector()
        vector_pi2_1.SetPtEtaPhiM(row['pi2_pt_1'], row['pi2_eta_1'], row['pi2_phi_1'], row['pi2_mass_1'])
        cartesian_vector_pi2_1 = ROOT.TLorentzVector(vector_pi2_1.Px(), vector_pi2_1.Py(), vector_pi2_1.Pz(), vector_pi2_1.E())
        pi2_1 = Particle(cartesian_vector_pi2_1, row['pi2_charge_1'])

        vector_pi2_2 = ROOT.TLorentzVector()
        vector_pi2_2.SetPtEtaPhiM(row['pi2_pt_2'], row['pi2_eta_2'], row['pi2_phi_2'], row['pi2_mass_2'])
        cartesian_vector_pi2_2 = ROOT.TLorentzVector(vector_pi2_2.Px(), vector_pi2_2.Py(), vector_pi2_2.Pz(), vector_pi2_2.E())
        pi2_2 = Particle(cartesian_vector_pi2_2, row['pi2_charge_2'])

        vector_pi3_1 = ROOT.TLorentzVector()
        vector_pi3_1.SetPtEtaPhiM(row['pi3_pt_1'], row['pi3_eta_1'], row['pi3_phi_1'], row['pi3_mass_1'])
        cartesian_vector_pi3_1 = ROOT.TLorentzVector(vector_pi3_1.Px(), vector_pi3_1.Py(), vector_pi3_1.Pz(), vector_pi3_1.E())
        pi3_1 = Particle(cartesian_vector_pi3_1, row['pi3_charge_1'])

        vector_pi3_2 = ROOT.TLorentzVector()
        vector_pi3_2.SetPtEtaPhiM(row['pi3_pt_2'], row['pi3_eta_2'], row['pi3_phi_2'], row['pi3_mass_2'])
        cartesian_vector_pi3_2 = ROOT.TLorentzVector(vector_pi3_2.Px(), vector_pi3_2.Py(), vector_pi3_2.Pz(), vector_pi3_2.E())
        pi3_2 = Particle(cartesian_vector_pi3_2, row['pi3_charge_2'])

        vector_pi0_1 = ROOT.TLorentzVector()
        vector_pi0_1.SetPtEtaPhiM(row['pi0_pt_1'], row['pi0_eta_1'], row['pi0_phi_1'], row['pi0_mass_1'])
        cartesian_vector_pi0_1 = ROOT.TLorentzVector(vector_pi0_1.Px(), vector_pi0_1.Py(), vector_pi0_1.Pz(), vector_pi0_1.E())
        pi0_1 = Particle(cartesian_vector_pi0_1, row['pi0_charge_1'])

        vector_pi0_2 = ROOT.TLorentzVector()
        vector_pi0_2.SetPtEtaPhiM(row['pi0_pt_2'], row['pi0_eta_2'], row['pi0_phi_2'], row['pi0_mass_2'])
        cartesian_vector_pi0_2 = ROOT.TLorentzVector(vector_pi0_2.Px(), vector_pi0_2.Py(), vector_pi0_2.Pz(), vector_pi0_2.E())
        pi0_2 = Particle(cartesian_vector_pi0_2, row['pi0_charge_2'])

        if row['decayMode_1'] == 0 and row['decayMode_2'] == 0:
            aco_angle_0_0 = AcoplanarityAngle(ip1,ip2,pi_1.Fourvector,pi_2.Fourvector,type1="IP",type2="IP")
            dict['aco_angle_0_0'].append(aco_angle_0_0)
        else:
            dict['aco_angle_0_0'].append(-9999.0)

        if row['decayMode_1'] == 0 and row['decayMode_2'] == 1:
            aco_angle_0_1 = AcoplanarityAngle(ip1,pi0_2.Fourvector,pi_1.Fourvector,pi_2.Fourvector,type1="IP",type2="NP")
            dict['aco_angle_0_1'].append(aco_angle_0_1)
        else:
            dict['aco_angle_0_1'].append(-9999.0)

        if row['decayMode_1'] == 1 and row['decayMode_2'] == 0:
            aco_angle_1_0 = AcoplanarityAngle(pi0_1.Fourvector,ip2, pi_1.Fourvector, pi_2.Fourvector,type1="NP",type2="IP")
            dict['aco_angle_1_0'].append(aco_angle_1_0)
        else:
            dict['aco_angle_1_0'].append(-9999.0)

        if row['decayMode_1'] == 1 and row['decayMode_2'] == 1:
            aco_angle_1_1 = AcoplanarityAngle(pi0_1.Fourvector,pi0_2.Fourvector,pi_1.Fourvector,pi_2.Fourvector,type1="NP",type2="NP")
            dict['aco_angle_1_1'].append(aco_angle_1_1)
        else:
            dict['aco_angle_1_1'].append(-9999.0)

        if row['decayMode_1'] == 0 and row['decayMode_2'] == 10:
            hads = sortA1(pi_2,pi2_2,pi3_2)
            aco_angle_0_10 = AcoplanarityAngle(ip1,hads[0].Fourvector,pi_1.Fourvector,hads[1].Fourvector,type1="IP",type2="A1")
            dict['aco_angle_0_10'].append(aco_angle_0_10)
        else:
            dict['aco_angle_0_10'].append(-9999.0)

        if row['decayMode_1'] == 10 and row['decayMode_2'] == 0:
            hads = sortA1(pi_1,pi2_1,pi3_1)
            aco_angle_10_0 = AcoplanarityAngle(hads[0].Fourvector,ip2,hads[1].Fourvector,pi_2.Fourvector,type1="A1",type2="IP")
            dict['aco_angle_10_0'].append(aco_angle_10_0)
        else:
            dict['aco_angle_10_0'].append(-9999.0)

        if row['decayMode_1'] == 1 and row['decayMode_2'] == 10:
            hads = sortA1(pi_2,pi2_2,pi3_2)
            aco_angle_1_10 = AcoplanarityAngle(pi0_1.Fourvector, hads[0].Fourvector, pi_1.Fourvector, hads[1].Fourvector, type1="NP", type2="A1")
            dict['aco_angle_1_10'].append(aco_angle_1_10)
        else:
            dict['aco_angle_1_10'].append(-9999.0)

        if row['decayMode_1'] == 10 and row['decayMode_2'] == 1:
            hads = sortA1(pi_1,pi2_1,pi3_1)
            aco_angle_10_1 = AcoplanarityAngle(hads[0].Fourvector, pi0_2.Fourvector, hads[1].Fourvector, pi_2.Fourvector, type1="A1", type2="NP")
            dict['aco_angle_10_1'].append(aco_angle_10_1)
        else:
            dict['aco_angle_10_1'].append(-9999.0)

        if row['decayMode_1'] == 10 and row['decayMode_2'] == 10:

            hads1 = [pi_1,pi2_1,pi3_1]
            hads2 = [pi_2,pi2_2,pi3_2]

            sv_minus_pv_1 = sv1.Vect().Unit()
            sv_minus_pv_2 = sv2.Vect().Unit()

            if row['charge_1'] > 0:
                hads2 = [pi_1,pi2_1,pi3_1]
                hads1 = [pi_2,pi2_2,pi3_2]

                sv_minus_pv_2 = sv1.Vect().Unit()
                sv_minus_pv_1 = sv2.Vect().Unit()

            aco_angle_11_11 = AcoplanarityAngle_A1A1(sv_minus_pv_1,sv_minus_pv_2,hads1,hads2)
            dict['aco_angle_10_10'].append(aco_angle_11_11)
        else:
            dict['aco_angle_10_10'].append(-9999.0)

        bar()

del df

df_new = pd.DataFrame(dict)

# read the original tree
tree = file['ntuple']
# convert tree to pandas dataframe
df_original = tree.pandas.df()

# check if any columns have a boolean dtype and convert it to int
#for col in df_original.columns:
#    if df_original[col].dtype == bool:
#        df_original[col] = df_original[col].astype(int)

# merge the two dataframes
df_final = pd.concat([df_original, df_new], axis=1)

data_dict = {}
type_dict = {}
for column in df_final.columns:
    data_dict[column] = df_final[column].to_numpy()
    dtype = data_dict[column].dtype
    type_dict[column] = dtype

root_file_path = filepath.replace(".root", "_cp.root")
with uproot3.recreate(root_file_path) as root_file:
    root_file["ntuple"] = uproot3.newtree(type_dict)
    root_file["ntuple"].extend(data_dict)

