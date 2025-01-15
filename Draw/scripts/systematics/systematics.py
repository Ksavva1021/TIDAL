# script to generate systematics dictionary
from collections import OrderedDict


def generate_systematics_dict(specific_era='Run3_2022', specific_channel='mt'):
    systematics = OrderedDict()

    # Muon ID/Isolation systematics
    # ----------------------------------------------------------------------------------------------------
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
            histogram_name = 'syst_muon_' + extension + updown.split('_')[-1]

            if specific_channel == 'mm':
                systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading}) * ({formula_subleading})", [], False)
            elif specific_channel == 'mt':
                systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading})", [], False)

        del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # Electron ID systematics
    # ----------------------------------------------------------------------------------------------------
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
        histogram_name = 'syst_electron_id' + updown.split('_')[-1]

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

    kinds = ['stat1','stat2','syst_TES_era_dm']
    eras = ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]
    decay_modes = ["0", "1", "10", "11"]

    for kind in kinds:
        for era in eras:
            for dm in decay_modes:
                up_weights = []
                down_weights = []
                for obj_index, obj_type in enumerate(specific_channel):
                    if obj_type == 't':
                        up_var = f'w_Tau_ID_{obj_index+1}_{kind}_Up'
                        down_var = f'w_Tau_ID_{obj_index+1}_{kind}_Down'

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

                if kind != "syst_TES_era_dm":
                    systematic_name = f'Tau_ID_{kind}_DM{dm}_{era}'
                    histogram_name = f'syst_tau_id_{kind}_DM{dm}_{era}'
                else:
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
                up_var = f'w_Tau_ID_{obj_index+1}_syst_era_Up'
                down_var = f'w_Tau_ID_{obj_index+1}_syst_era_Down'

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
            up_var = f'w_Tau_ID_{obj_index+1}_syst_all_eras_Up'
            down_var = f'w_Tau_ID_{obj_index+1}_syst_all_eras_Down'

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
    eta_bins = ["0.0", "0.4", "0.8", "1.2", "1.7", "2.4"]
    eras = ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]

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

    kinds = {
        '1prong': '0PI',
        '1prong1pizero': '1PI',
        '3prong': '3PRONG',
        '3prong1pizero': '3PRONG1PI0'
    }

    # Genuine Taus, Genuine electrons/muons misidentified as taus
    prefixes = ['Tau_EnergyScale_TSCALE_', 'Tau_EnergyScale_ESCALE_', 'Tau_EnergyScale_MUSCALE_']

    for name, folder_suffix in kinds.items():
        for prefix in prefixes:
            for updown in ['up', 'down']:
                systematic_name = 'syst_tau_escale_' + name + '_' + updown
                folder_name = prefix + folder_suffix + '_' + updown.upper()
                histogram_name = 'syst_tau_escale_' + name + updown.capitalize()

                if specific_channel in ["et","mt","tt"]:
                    systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Scale systematics
    # ----------------------------------------------------------------------------------------------------

    prefix = 'jec_syst_Total'
    name = 'syst_jet_scale_Total'
    for updown in ['up', 'down']:
        systematic_name = name + '_' + updown
        folder_name = prefix + '_' + updown
        histogram_name = '_' + prefix + updown.capitalize()
        systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Resolution systematics
    # ----------------------------------------------------------------------------------------------------

    prefix = 'jer_syst'
    name = 'syst_jet_resolution'
    for updown in ['up', 'down']:
        systematic_name = name + '_' + updown
        folder_name = prefix + '_' + updown
        histogram_name = '_' + prefix + updown.capitalize()
        systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    # Electron Energy Scale + Smearing systematics
    # ----------------------------------------------------------------------------------------------------
    for kind in ['Scale', 'Smearing']:
        for updown in ['up', 'down']:
            systematic_name = 'syst_electron_' + kind.lower() + '_' + updown
            folder_name = 'Electron_' + kind + '_' + updown
            histogram_name = '_' +'syst_electron_' + kind.lower() + updown.capitalize()

            if specific_channel in ["ee","et"]:
                systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', [], False)

    # ----------------------------------------------------------------------------------------------------

    return systematics
