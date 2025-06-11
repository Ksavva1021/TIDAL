# script to generate systematics dictionary
from collections import OrderedDict


def generate_systematics_dict(specific_era='Run3_2022', specific_channel='mt', specific_systematic='MuonID', specific_name=''):
    systematics = OrderedDict()

    # Muon ID/Isolation systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Muon_ID' or specific_systematic == 'Muon_Isolation':
        for kind in ['ID', 'Isolation']:
            up_var = f'w_Muon_{kind}_Up'
            down_var = f'w_Muon_{kind}_Down'

            for updown in [up_var,down_var]:
                formula_leading = (
                    f"(({updown}) * (genPartFlav_1 == 1 || genPartFlav_1 == 15) + "
                    f"(!(genPartFlav_1 == 1 || genPartFlav_1 == 15)))"
                )

                formula_subleading = (
                    f"(({updown}) * (genPartFlav_2 == 1 || genPartFlav_2 == 15) + "
                    f"(!(genPartFlav_2 == 1 || genPartFlav_2 == 15)))"
                )

                if kind == "ID":
                    extension = 'id'
                elif kind == "Isolation":
                    extension = 'iso'

                systematic_name = 'Muon_' + kind + '_' + updown.split('_')[-1]
                if specific_name == '':
                    histogram_name = 'syst_muon_' + extension + updown.split('_')[-1]
                else:
                    histogram_name = '_' + specific_name + updown.split('_')[-1]

                if specific_channel == 'mm':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading}) * ({formula_subleading})", [], False)
                elif specific_channel == 'mt':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading})", [], False)

            del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # Electron ID systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Electron_ID':
        up_var = 'w_Electron_ID_Up'
        down_var = 'w_Electron_ID_Down'

        for updown in [up_var,down_var]:
            formula_leading = (
                f"(({updown}) * (genPartFlav_1 == 1 || genPartFlav_1 == 15) + "
                f"(!(genPartFlav_1 == 1 || genPartFlav_1 == 15)))"
            )

            formula_subleading = (
                f"(({updown}) * (genPartFlav_2 == 1 || genPartFlav_2 == 15) + "
                f"(!(genPartFlav_2 == 1 || genPartFlav_2 == 15)))"
            )

            systematic_name = 'Electron_ID' + updown.split('_')[-1]
            if specific_name == '':
                histogram_name = 'syst_electron_id' + updown.split('_')[-1]
            else:
                histogram_name = '_' + specific_name + updown.split('_')[-1]

            if specific_channel == 'ee':
                systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading}) * ({formula_subleading})", [], False)
            elif specific_channel == 'et':
                systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading})", [], False)

        del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # Tau ID systematics
    # ----------------------------------------------------------------------------------------------------

    # First Kind: stat1, stat2, syst_TES_era_dm
    # should be uncorrelated across DMs and eras
    if specific_systematic == 'Tau_ID_PNet':
        kinds = ['stat1','stat2']
        eras = ["Run3_2023", "Run3_2022EE", "Run3_2022", "Run3_2023BPix"]
        decay_modes = ["0", "1", "2", "10", "11"]

        for kind in kinds:
            for era in eras:
                for dm in decay_modes:
                    up_weights = []
                    down_weights = []
                    for obj_index, obj_type in enumerate(specific_channel):
                        if obj_type == 't':
                            up_var = f'w_Tau_ID_PNet_{obj_index+1}_{kind}_Up'
                            down_var = f'w_Tau_ID_PNet_{obj_index+1}_{kind}_Down'

                            if specific_era == era:
                                formula = (
                                    f"((variation_to_replace) * (decayMode_{obj_index+1} == {dm}) + "
                                    f"(!(decayMode_{obj_index+1} == {dm})))"
                                )
                            else:
                                formula = (
                                    f"((1) * (decayMode_{obj_index+1} == {dm}) + "
                                    f"(!(decayMode_{obj_index+1} == {dm})))"
                                )

                            up_weights.append(formula.replace('variation_to_replace', up_var))
                            down_weights.append(formula.replace('variation_to_replace', down_var))

                            del up_var, down_var

                    systematic_name = f'Tau_ID_{kind.replace("_era_", "_")}_DM{dm}_{era}'
                    histogram_name = f'syst_tau_id_{kind.replace("_era_", "_")}_DM{dm}_{era}'

                    if specific_channel in ["et","mt","tt"]:
                        systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), [], False)
                        systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), [], False)

        del up_weights, down_weights

        # Second Kind: syst_era
        # should be correlated across DMs but uncorrelated across eras

        for era in eras:
            up_weights = []
            down_weights = []
            for obj_index, obj_type in enumerate(specific_channel):
                if obj_type == 't':
                    up_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_era_Up'
                    down_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_era_Down'

                    if specific_era == era:
                        up_weights.append(f"({up_var})")
                        down_weights.append(f"({down_var})")
                    else:
                        up_weights.append("(1)")
                        down_weights.append("(1)")

                    del up_var, down_var

            systematic_name = f'Tau_ID_syst_era_{era}'
            histogram_name = f'syst_tau_id_{era}'
            if specific_channel in ["et","mt","tt"]:
                systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), [], False)
                systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), [], False)

        del up_weights, down_weights

        # Third Kind: syst_all_eras
        # should be correlated across DMs and eras

        up_weights = []
        down_weights = []
        for obj_index, obj_type in enumerate(specific_channel):
            if obj_type == 't':
                up_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_all_eras_Up'
                down_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_all_eras_Down'

                up_weights.append(f"({up_var})")
                down_weights.append(f"({down_var})")

                del up_var, down_var

        systematic_name = 'Tau_ID_syst_all_eras'
        histogram_name = 'syst_tau_id_all_eras'
        if specific_channel in ["et","mt","tt"]:
            systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), [], False)
            systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), [], False)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Fakerate (e) systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Tau_FakeRate_e':
        eta_bins = ["0.0", "1.5", "2.5"]
        eras = ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]
        for era in eras:
            for i, eta in enumerate(eta_bins):
                up_weights = []
                down_weights = []
                for obj_index, obj_type in enumerate(specific_channel):
                    if obj_type == 't':
                        up_var = f"w_Tau_e_FakeRate_{obj_index+1}_Up"
                        down_var = f"w_Tau_e_FakeRate_{obj_index+1}_Down"

                        if specific_era == era:
                            if i == len(eta_bins) - 1:
                                continue
                            else:
                                formula = (
                                    f"((variation_to_replace * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                                    f"(fabs(eta_{obj_index+1}) >= {eta_bins[i+1]}))"
                                )
                        else:
                            if i == len(eta_bins) - 1:
                                continue
                            else:
                                formula = (
                                    f"((1 * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                                    f"(fabs(eta_{obj_index+1}) >= {eta_bins[i+1]}))"
                                )

                        up_weights.append(formula.replace('variation_to_replace', up_var))
                        down_weights.append(formula.replace('variation_to_replace', down_var))

                        del up_var, down_var

                systematic_name = f'Tau_e_FakeRate_{era}_eta_{eta}'
                histogram_name = f'syst_etau_fakerate_{era}_eta_{eta}'

                if specific_channel in ["et","mt","tt"]:
                    systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), [], False)
                    systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), [], False)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Fakerate (mu) systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Tau_FakeRate_mu':
        eta_bins = ["0.0", "0.4", "0.8", "1.2", "1.7", "2.4"]
        eras = [ "Run3_2023", "Run3_2022", "Run3_2022EE", "Run3_2023BPix"]

        for era in eras:
            for i, eta in enumerate(eta_bins):
                up_weights = []
                down_weights = []
                for obj_index, obj_type in enumerate(specific_channel):
                    if obj_type == 't':
                        up_var = f"w_Tau_mu_FakeRate_{obj_index+1}_Up"
                        down_var = f"w_Tau_mu_FakeRate_{obj_index+1}_Down"

                        if specific_era == era:
                            if i == len(eta_bins) - 1:
                                formula = (
                                    f"((variation_to_replace * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}))) + "
                                    f"(fabs(eta_{obj_index+1}) < {eta_bins[i-1]}))"
                                )
                            else:
                                formula = (
                                    f"((variation_to_replace * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                                    f"(fabs(eta_{obj_index+1}) >= {eta_bins[i+1]}))"
                                )
                        else:
                            if i == len(eta_bins) - 1:
                                formula = (
                                    f"((1 * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}))) + "
                                    f"(fabs(eta_{obj_index+1}) < {eta_bins[i-1]}))"
                                )
                            else:
                                formula = (
                                    f"((1 * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                                    f"(fabs(eta_{obj_index+1}) >= {eta_bins[i+1]}))"
                                )

                        up_weights.append(formula.replace('variation_to_replace', up_var))
                        down_weights.append(formula.replace('variation_to_replace', down_var))

                        del up_var, down_var

                systematic_name = f'Tau_mu_FakeRate_{era}_eta_{eta}'
                histogram_name = f'syst_mutau_fakerate_{era}_eta_{eta}'
                if specific_channel in ["et","mt","tt"]:
                    systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), [], False)
                    systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), [], False)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Energy Scale systematics (This are recommended to be uncorrelated across eras but we will leave it for now)
    # ----------------------------------------------------------------------------------------------------

    if specific_systematic in ['Tau_EnergyScale_TSCALE', 'Tau_EnergyScale_ESCALE', 'Tau_EnergyScale_MUSCALE']:
        kinds = {
            '1prong': '1PRONG',
            '1prong1pizero': '1PRONG_1PI0',
            '1prong2pizero': '1PRONG_2PI0',
            '3prong': '3PRONG',
            '3prong1pizero': '3PRONG_1PI0'
        }

        # Genuine Taus
        if specific_systematic == 'Tau_EnergyScale_TSCALE':
            prefixes = ['Tau_EnergyScale_TSCALE_']
            samples_to_skip = ["QCD"]
        # Genuine electrons misidentified as taus
        elif specific_systematic == 'Tau_EnergyScale_ESCALE':
            prefixes = ['Tau_EnergyScale_ESCALE_']
            samples_to_skip = ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W']
        # Genuine muons misidentified as taus
        elif specific_systematic == 'Tau_EnergyScale_MUSCALE':
            prefixes = ['Tau_EnergyScale_MUSCALE_']
            samples_to_skip = ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W']

        # temporary patch for derivation of tau id scale factor systematics
        # Tau_EnergyScale_forTauIDSFs_TSCALE, Tau_EnergyScale_forTauIDSFs_ESCALE, Tau_EnergyScale_forTauIDSFs_MUSCALE,
        for index in range(len(prefixes)):
            # TODO: CHANGE BACK TO no PNet
            prefixes[index] = prefixes[index].replace("Tau_EnergyScale", "Tau_EnergyScale_forTauIDSFs_PNet")

        for name, folder_suffix in kinds.items():
            for prefix in prefixes:
                for updown in ['up', 'down']:
                    scale_type = prefix.removesuffix('_').split('_')[-1]
                    systematic_name = 'syst_tau_escale_' + scale_type + '_' + name + '_' + updown
                    folder_name = prefix + folder_suffix + '_' + updown
                    if specific_name == '':
                        histogram_name = 'syst_tau_escale_' + scale_type + '_' + name + updown.capitalize()
                    else:
                        histogram_name = '_' + specific_name.replace("*group", name).replace("*year", specific_era) + updown.capitalize()

                    if specific_channel in ["et","mt","tt"]:
                        systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', samples_to_skip, False)

    if specific_systematic == 'Tau_EnergyScale_JSCALE':

        samples_to_skip = [
            'ZJ','ZL','ZLL','ZTT'
            'VVT','VVJ',
            'TTT','TTJ',
            'QCD','signal',
        ]

        prefix = 'Tau_EnergyScale_JSCALE_'
        # temporary patch for derivation of tau id scale factor systematics
        # Tau_EnergyScale_forTauIDSFs_JSCALE
        # TODO: CHANGE BACK TO no PNet
        prefix = prefix.replace("Tau_EnergyScale", "Tau_EnergyScale_forTauIDSFs_PNet")
        for updown in ['up', 'down']:
            systematic_name = 'syst_tau_escale_jscale_' + updown
            folder_name = prefix + updown
            if specific_name == '':
                histogram_name = 'syst_tau_escale_jscale' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name.replace("*year", specific_era) + updown.capitalize()

            if specific_channel in ["et","mt","tt"]:
                systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', samples_to_skip, False)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Scale systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Jet_EnergyScale_Total':
        prefix = 'jec_syst_Total'
        name = 'syst_jet_scale_Total'
        for updown in ['up', 'down']:
            systematic_name = name + '_' + updown
            folder_name = prefix + '_' + updown
            if specific_name == '':
                histogram_name = '_' + prefix + updown.capitalize()
            else:
                histogram_name = '_' + specific_name.replace("*year",specific_era) + updown.capitalize()
            systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Resolution systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Jet_EnergyResolution':
        prefix = 'jer_syst'
        name = 'syst_jet_resolution'
        for updown in ['up', 'down']:
            systematic_name = name + '_' + updown
            folder_name = prefix + '_' + updown
            if specific_name == '':
                histogram_name = '_' + prefix + updown.capitalize()
            else:
                histogram_name = '_' + specific_name.replace("*year",specific_era) + updown.capitalize()
            systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # Electron Energy Scale + Smearing systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Electron_Scale' or specific_systematic == 'Electron_Smearing':
        for kind in ['Scale', 'Smearing']:
            for updown in ['up', 'down']:
                systematic_name = 'syst_electron_' + kind.lower() + '_' + updown
                folder_name = 'Electron_' + kind + '_' + updown
                if specific_name == '':
                    histogram_name = '_' +'syst_electron_' + kind.lower() + updown.capitalize()
                else:
                    histogram_name = '_' + specific_name + updown.capitalize()

                if specific_channel in ["ee","et"]:
                    systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # TTbar pT reweighting systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'TTbar_Shape':
        samples_to_skip = [
            "ZTT", "ZLL", "ZL", "ZJ",
            "VV", "VVT", "VVJ",
            "W","signal", 'QCD'
        ]
        up_var = 'w_Top_pt_Reweighting'
        down_var = '(1/w_Top_pt_Reweighting)'

        for updown in ["up", "down"]:
            systematic_name = 'syst_ttbar_shape_' + updown
            if specific_name == '':
                histogram_name = '_syst_ttbar_shape' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            weight_updown = up_var if updown == "up" else down_var
            systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", samples_to_skip, False)

        del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # DY pT reweighting systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'DY_Shape':
        samples_to_skip = [
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "W","signal", 'QCD'
        ]
        up_var = 'w_Zpt_Reweighting'
        down_var = '(1/w_Zpt_Reweighting)'

        for updown in ["up", "down"]:
            systematic_name = 'syst_dy_shape_' + updown
            if specific_name == '':
                histogram_name = '_syst_dy_shape' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            weight_updown = up_var if updown == "up" else down_var
            systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", samples_to_skip, False)

        del up_var, down_var

    if specific_systematic == 'Fake_Flat_Uncertainty':
        samples_to_skip = [
            "TT", "TTT", "TTJ",
            "ZTT", "ZLL", "ZL", "ZJ",
            "VV", "VVT", "VVJ",
            "W","signal"
        ]

        for updown in ["up", "down"]:
            systematic_name = 'flat_fake_sub_' + updown
            if specific_name == '':
                histogram_name = '_flat_fake_sub' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", samples_to_skip, False)


    # ----------------------------------------------------------------------------------------------------

    # QCD Background systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "QCD_Background":
        samples_to_skip = [
            "ZTT", "ZLL", "ZL", "ZJ",
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "W","signal"
        ]
        for updown in ['up','down']:
            systematic_name = 'qcd_sub_' + updown
            if specific_name == '':
                histogram_name = '_qcd_sub' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            systematics[systematic_name] = ('nominal', histogram_name, 'weight_to_replace', samples_to_skip, False)

    return systematics
