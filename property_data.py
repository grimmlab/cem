
"""
This file mainly contains classes with property data (either parameters or experimental lle data) for
various systems and sources.
"""

import itertools
import numpy as np
import os
import pandas as pd
import importlib

import thermo_models
import analysis_class


def convert_mass_fr_to_molar_fr(mass_fractions, molar_masses):
    molar_fractions = np.zeros(len(mass_fractions))
    denominator = 0
    for i in range(len(mass_fractions)):
        denominator = denominator + (mass_fractions[i] / molar_masses[i])

    for i in range(len(mass_fractions)):
        molar_fractions[i] = (mass_fractions[i] / molar_masses[i]) / denominator

    return molar_fractions


class InitClass:
    """
    Takes names and some other inputs and sets all the parameters, which are needed to set
    up a lle analysis (for example the interactions parameters are immediately processed to the
    relevant ge_model here). The lle analysis is then executed.

    NRTL parameters are in classes NRTL_ + name_for_data
    UNIQUAC parameters are in classes UNIQUAC_ + name_for_data
    Experimental data is in classes Experimental_ + name_for_data
    """
    def __init__(self, names_list, point_discretization_rec_steps, name_for_data,
                 temperature_kelvin, store_lle, deciding_index, actors_parallelized,
                 ge_model_name, tau_given, experimental_data_given):
        self.names_list = names_list
        self.num_comp = len(self.names_list)
        self.point_discretization_rec_steps = point_discretization_rec_steps
        self.name_for_data = name_for_data
        self.store_lle = store_lle
        self.temperature_kelvin = temperature_kelvin
        self.deciding_index = deciding_index
        self.actors_para = actors_parallelized

        # to be able to get instances of classes by strings later on
        my_module = importlib.import_module("property_data")

        if ge_model_name == "NRTL":
            parameter_name = "NRTL_" + name_for_data
            temp_class = getattr(my_module, parameter_name)
            nrtl_collection = temp_class()
            specific_paras = nrtl_collection.combine_interactions(names_list)

            alphas = specific_paras[0]
            a = specific_paras[1]
            b = specific_paras[2]
            e = specific_paras[3]
            f = specific_paras[4]
            tau = specific_paras[5]

            self.ge_model = thermo_models.NRTL(alphas, a, b, e, f, tau, tau_given=tau_given)

        elif ge_model_name == "UNIQUAC":
            parameter_name = "UNIQUAC_" + name_for_data
            temp_class = getattr(my_module, parameter_name)
            uniquac_collection = temp_class()
            specific_paras = uniquac_collection.combine_interactions(names_list)

            a = specific_paras[3]
            b = specific_paras[4]
            c = specific_paras[5]
            d = specific_paras[6]
            e = specific_paras[7]

            q = specific_paras[0]
            q_tilde = specific_paras[1]
            r = specific_paras[2]

            self.ge_model = thermo_models.UNIQUAC(a, b, c, d, e, r, q, q_tilde)

        if experimental_data_given:
            experimental_data_name = "Experimental_" + name_for_data
            temp_class2 = getattr(my_module, experimental_data_name)
            self.experimental_data = temp_class2()

        else:
            self.experimental_data = None

        self.analysis = analysis_class.Analysis(num_comp=self.num_comp, names_list=self.names_list,
                                                point_discretization_rec_steps=self.point_discretization_rec_steps,
                                                ge_model=self.ge_model, store_lle=self.store_lle,
                                                experimental_data=self.experimental_data,
                                                temperature_kelvin=self.temperature_kelvin,
                                                deciding_index=self.deciding_index,
                                                actors_para=self.actors_para,
                                                general_name=self.name_for_data)


class UNIQUAC_Pure_Ryll2009_1:
    """
    Parameters for pure components listed in dissertation of Ryll. Not available in digital form,
    therefore listed by hand.
    """
    def __init__(self):
        self.names = ["water", "benzene", "nitromethane", "1_hexanol", "ethanol",
                      "ethylene_glycol", "ethyl_formate", "furfural", "isobutanol",
                      "phenol", "acetone", "n_heptane", "acetonitrile", "1_butanol",
                      "acetic_acid", "1_butyl_acetate"]

        self.q = np.array([1.4, 2.4, 1.876, 4.132, 1.972, 2.248, 2.576, 2.5, 3.048,
                           2.716, 2.336, 4.396, 1.724, 3.048, 2.072, 4.196])
        self.q_tilde = np.array([1.4, 2.4, 1.876, 4.132, 1.972, 2.248, 2.576, 2.5, 3.048,
                                2.716, 2.336, 4.396, 1.724, 3.048, 2.072, 4.196])
        self.r = np.array([0.92, 3.19051, 2.00857, 4.8029, 2.10547, 2.4087, 2.8174,
                           3.17205, 3.45353, 3.54647, 2.5735, 5.17403, 1.87014,
                           3.45419, 2.19512, 4.82729])

    def get_pure_para(self, name):
        for ind, name_2 in enumerate(self.names):
            if str(name) == str(name_2):
                return [self.q[ind], self.q_tilde[ind], self.r[ind]]

        return None


class UNIQUAC_Pure_Ryll2009_2:
    """
    Parameters for pure components listed in dissertation of Ryll. Not available in digital form,
    therefore listed by hand. This class is just needed as those components are also partly present
    with other binary interactions in the first class.
    """
    def __init__(self):
        self.names = ["water", "nitromethane", "nonanol"]

        self.q = np.array([1.4, 1.876, 5.752])
        self.q_tilde = np.array([1.4, 1.876, 5.752])
        self.r = np.array([0.92, 2.00857, 6.82927])

    def get_pure_para(self, name):
        for ind, name_2 in enumerate(self.names):
            if str(name) == str(name_2):
                return [self.q[ind], self.q_tilde[ind], self.r[ind]]

        return None


class UNIQUAC_Ryll2009_1:
    def __init__(self):
        self.binary_names = []

        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.e = []

        self.add_pair("nitromethane", "1_hexanol", np.zeros(2), np.array([-59.312, -342.15]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("benzene", "ethanol", np.zeros(2), np.array([-177.34, 171.7]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("benzene", "water", np.zeros(2), np.array([-860.81, -369.01]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("ethanol", "water", np.zeros(2), np.array([185.3, 167.38]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("water", "nitromethane", np.zeros(2), np.array([-229.72, -421.7]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("water", "1_hexanol", np.zeros(2), np.array([-276.54, -109.38]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("water", "ethylene_glycol", np.zeros(2), np.array([221.13, 39.881]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("water", "ethyl_formate", np.zeros(2), np.array([-170.57, -363.06]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("ethylene_glycol", "ethyl_formate", np.zeros(2), np.array([-211.54, -104.7]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("furfural", "water", np.zeros(2), np.array([-185.52, -103.14]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("furfural", "isobutanol", np.zeros(2), np.array([67.704, -401.61]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("isobutanol", "water", np.zeros(2), np.array([-52.254, -207.5]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("phenol", "acetone", np.zeros(2), np.array([278.76, 248.26]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("phenol", "water", np.zeros(2), np.array([187.88, -290.44]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("acetone", "water", np.zeros(2), np.array([-30.5, -48.44]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("n_heptane", "benzene", np.zeros(2), np.array([108.5, 120.28]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("n_heptane", "water", np.zeros(2), np.array([-1283.20, -575.85]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("n_heptane", "acetonitrile", np.zeros(2), np.array([-579.72, -34.444]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("benzene", "acetonitrile", np.zeros(2), np.array([-29.417, 209.79]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("acetonitrile", "water", np.zeros(2), np.array([-156.47, -91.529]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("1_butanol", "acetic_acid", np.array([-0.651, 0.639]), np.array([-31.89, -18.5]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("1_butanol", "1_butyl_acetate", np.array([-0.073, 0.608]), np.array([60.29, -327.08]),
                      np.zeros(2), np.zeros(2), np.zeros(2))
        self.add_pair("1_butanol", "water", np.array([-2.690, -5.284]), np.array([168.67, 1016.79]), np.zeros(2),
                      np.array([0.00709, 0.00356]), np.zeros(2))
        self.add_pair("acetic_acid", "1_butyl_acetate", np.array([-1.278, 1.967]), np.array([581.68, -950.36]),
                      np.zeros(2), np.zeros(2), np.zeros(2))
        self.add_pair("acetic_acid", "water", np.array([1.316, 0.770]), np.array([-1074.60, 33.816]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("1_butyl_acetate", "water", np.array([1.093, 0.252]), np.array([-872.21, -198.78]), np.zeros(2),
                      np.zeros(2), np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        q_list = np.zeros(len(names))
        q_tilde_list = np.zeros(len(names))
        r_list = np.zeros(len(names))

        pure_storage = UNIQUAC_Pure_Ryll2009_1()

        for i, name in enumerate(names):
            paras = pure_storage.get_pure_para(name)
            if paras is None:
                return None

            else:
                q_list[i] = paras[0]
                q_tilde_list[i] = paras[1]
                r_list[i] = paras[2]

        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        c = np.zeros((len(names), len(names)))
        d = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))

        matrices = [a, b, c, d, e]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return [q_list, q_tilde_list, r_list] + matrices

    def add_pair(self, name_1, name_2, a, b, c, d, e):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.a.append(a)
        self.b.append(b)
        self.c.append(c)
        self.d.append(d)
        self.e.append(e)


class UNIQUAC_Ryll2009_2:
    def __init__(self):
        self.binary_names = []

        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.e = []

        self.add_pair("nonanol", "nitromethane", np.zeros(2), np.array([-375.55, -45.711]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("nonanol", "water", np.zeros(2), np.array([-120.3, -366.53]), np.zeros(2),
                      np.zeros(2), np.zeros(2))
        self.add_pair("nitromethane", "water", np.zeros(2), np.array([-421.7, -229.72]), np.zeros(2),
                      np.zeros(2), np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        q_list = np.zeros(len(names))
        q_tilde_list = np.zeros(len(names))
        r_list = np.zeros(len(names))

        pure_storage = UNIQUAC_Pure_Ryll2009_2()

        for i, name in enumerate(names):
            paras = pure_storage.get_pure_para(name)
            if paras is None:
                return None

            else:
                q_list[i] = paras[0]
                q_tilde_list[i] = paras[1]
                r_list[i] = paras[2]

        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        c = np.zeros((len(names), len(names)))
        d = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))

        matrices = [a, b, c, d, e]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return [q_list, q_tilde_list, r_list] + matrices

    def add_pair(self, name_1, name_2, a, b, c, d, e):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.a.append(a)
        self.b.append(b)
        self.c.append(c)
        self.d.append(d)
        self.e.append(e)


class NRTL_Chen2000_1:
    def __init__(self):
        self.binary_names = []

        self.alphas = []
        self.a = []
        self.b = []
        self.e = []
        self.f = []
        self.tau = []

        self.add_pair("n_hexane", "benzene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([1.523, 90.89]))
        self.add_pair("n_hexane", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([2045, 870.6]))
        self.add_pair("benzene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([403.2, -103.8]))
        self.add_pair("n_hexane", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([347, -267.6]))
        self.add_pair("toluene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([562.3, -53.59]))
        self.add_pair("n_hexane", "xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([115.5, -98.56]))
        self.add_pair("xylene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([413.7, 145.1]))
        self.add_pair("n_octane", "benzene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([-163.5, 283.7]))
        self.add_pair("n_octane", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([2235, 1231]))
        self.add_pair("n_octane", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([115.5, -84.22]))
        self.add_pair("n_octane", "xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([-331.8, 364.5]))
        self.add_pair("n_octane", "n_hexane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("benzene", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("benzene", "xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("toluene", "xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        alphas = np.zeros((len(names), len(names)))
        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))
        f = np.zeros((len(names), len(names)))
        tau = np.zeros((len(names), len(names)))

        matrices = [alphas, a, b, e, f, tau]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return matrices

    def add_pair(self, name_1, name_2, alphas, a, b, e, f, tau):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.alphas.append(alphas)
        self.a.append(a)
        self.b.append(b)
        self.e.append(e)
        self.f.append(f)
        self.tau.append(tau)


class Experimental_Chen2000_1:
    def __init__(self):
        self.names = []
        self.feeds = []
        self.corresponding_phases = []
        self.temperatures = []

        temp = 298.15
        names = [["n_hexane", "benzene", "sulfolane"], ["n_hexane", "toluene", "sulfolane"],
                 ["n_hexane", "xylene", "sulfolane"], ["n_octane", "benzene", "sulfolane"],
                 ["n_octane", "toluene", "sulfolane"], ["n_octane", "xylene", "sulfolane"],
                 ["n_hexane", "n_octane", "benzene", "sulfolane"],
                 ["n_hexane", "benzene", "xylene", "sulfolane"],
                 ["n_octane", "toluene", "xylene", "sulfolane"],
                 ["n_hexane", "n_octane", "benzene", "toluene", "sulfolane"]]
        molar_fractions_given = True
        molar_masses = []
        sheet_names = ["Table_1_1", "Table_1_2", "Table_1_3",
                       "Table_2_1", "Table_2_2", "Table_2_3",
                       "Table_5_1", "Table_5_2", "Table_5_3",
                       "Table_7"]
        excel_path = os.path.join(os.getcwd(), "source_property_data", "chen_2000_1", "to_read.xlsx")

        for j in range(len(names)):
            data = pd.read_excel(excel_path, sheet_name=sheet_names[j], header=None, index_col=None)
            self.add_system(temp, names[j], molar_fractions_given, molar_masses, data)

    def add_system(self, temp, names, molar_fractions_given, molar_masses, data):
        matrix_to_read = data.to_numpy()
        num_comp = len(names)
        phase_1 = np.zeros((len(matrix_to_read), num_comp))
        phase_2 = np.zeros((len(matrix_to_read), num_comp))

        for i in range(len(matrix_to_read)):
            phase_1[i] = matrix_to_read[i][:num_comp]
            phase_2[i] = matrix_to_read[i][num_comp:2 * num_comp]

        if not molar_fractions_given:
            for i in range(len(phase_1)):
                phase_1[i] = convert_mass_fr_to_molar_fr(phase_1[i], molar_masses)
                phase_2[i] = convert_mass_fr_to_molar_fr(phase_2[i], molar_masses)

        corresponding_phases = []
        for i in range(len(phase_1)):
            corresponding_phases.append([phase_1[i], phase_2[i]])

        feeds = []
        for i in range(len(corresponding_phases)):
            feed = np.zeros(len(names))
            for phase in corresponding_phases[i]:
                feed = feed + phase

            feed = feed / sum(feed)
            feeds.append(feed)

        self.names.append(names)
        self.feeds.append(feeds)
        self.corresponding_phases.append(corresponding_phases)
        self.temperatures.append(temp)

    @staticmethod
    def compare_names(names_1, names_2):
        if len(names_1) != len(names_2):
            return False

        else:
            equal = [False] * len(names_1)
            for i, name_1 in enumerate(names_1):
                for j, name_2 in enumerate(names_2):
                    if name_1 == name_2:
                        equal[i] = True
                        break

            if sum(equal) == len(names_1):
                return True

            else:
                return False


class NRTL_Chen2001:
    def __init__(self):
        self.binary_names = []

        self.alphas = []
        self.a = []
        self.b = []
        self.e = []
        self.f = []
        self.tau = []

        self.add_pair("hexane", "heptane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("hexane", "octane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("hexane", "benzene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([1.523, 90.89]))
        self.add_pair("hexane", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([347, -267.6]))
        self.add_pair("hexane", "m_xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([115.5, -98.56]))
        self.add_pair("hexane", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([2045, 870.6]))
        self.add_pair("heptane", "octane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("heptane", "benzene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([-18.13, 98.93]))
        self.add_pair("heptane", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([359.8, -279.8]))
        self.add_pair("heptane", "m_xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([139.8, -58.49]))
        self.add_pair("heptane", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([1960, 870.6]))
        self.add_pair("octane", "benzene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([-163.5, 283.7]))
        self.add_pair("octane", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([115.5, -84.22]))
        self.add_pair("octane", "m_xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([-331.8, 364.5]))
        self.add_pair("octane", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([2235, 1231]))
        self.add_pair("benzene", "toluene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("benzene", "m_xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("benzene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([403.2, -103.8]))
        self.add_pair("toluene", "m_xylene", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([0, 0]))
        self.add_pair("toluene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([562.3, -53.59]))
        self.add_pair("m_xylene", "sulfolane", alphas=np.array([0.2, 0.2]), a=np.zeros(2),
                      b=np.zeros(2), e=np.zeros(2), f=np.zeros(2), tau=np.array([413.7, 145.1]))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        alphas = np.zeros((len(names), len(names)))
        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))
        f = np.zeros((len(names), len(names)))
        tau = np.zeros((len(names), len(names)))

        matrices = [alphas, a, b, e, f, tau]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return matrices

    def add_pair(self, name_1, name_2, alphas, a, b, e, f, tau):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.alphas.append(alphas)
        self.a.append(a)
        self.b.append(b)
        self.e.append(e)
        self.f.append(f)
        self.tau.append(tau)


class Experimental_Chen2001:
    def __init__(self):
        self.names = []
        self.feeds = []
        self.corresponding_phases = []
        self.temperatures = []

        temp = 298.15
        names = [["hexane", "heptane", "toluene", "sulfolane"],
                 ["heptane", "octane", "m_xylene", "sulfolane"],
                 ["heptane", "benzene", "toluene", "sulfolane"],
                 ["hexane", "heptane", "toluene", "m_xylene", "sulfolane"],
                 ["heptane", "octane", "benzene", "m_xylene", "sulfolane"]]
        molar_fractions_given = True
        molar_masses = []
        sheet_names = ["Table_1_1", "Table_1_2", "Table_1_3",
                       "Table_2_1", "Table_2_2"]
        excel_path = os.path.join(os.getcwd(), "source_property_data", "chen_2001", "to_read.xlsx")

        for j in range(len(names)):
            data = pd.read_excel(excel_path, sheet_name=sheet_names[j], header=None, index_col=None)
            self.add_system(temp, names[j], molar_fractions_given, molar_masses, data)

    def add_system(self, temp, names, molar_fractions_given, molar_masses, data):
        matrix_to_read = data.to_numpy()
        num_comp = len(names)
        phase_1 = np.zeros((len(matrix_to_read), num_comp))
        phase_2 = np.zeros((len(matrix_to_read), num_comp))

        for i in range(len(matrix_to_read)):
            phase_1[i] = matrix_to_read[i][:num_comp]
            phase_2[i] = matrix_to_read[i][num_comp:2 * num_comp]

        if not molar_fractions_given:
            for i in range(len(phase_1)):
                phase_1[i] = convert_mass_fr_to_molar_fr(phase_1[i], molar_masses)
                phase_2[i] = convert_mass_fr_to_molar_fr(phase_2[i], molar_masses)

        corresponding_phases = []
        for i in range(len(phase_1)):
            corresponding_phases.append([phase_1[i], phase_2[i]])

        feeds = []
        for i in range(len(corresponding_phases)):
            feed = np.zeros(len(names))
            for phase in corresponding_phases[i]:
                feed = feed + phase

            feed = feed / sum(feed)
            feeds.append(feed)

        self.names.append(names)
        self.feeds.append(feeds)
        self.corresponding_phases.append(corresponding_phases)
        self.temperatures.append(temp)

    @staticmethod
    def compare_names(names_1, names_2):
        if len(names_1) != len(names_2):
            return False

        else:
            equal = [False] * len(names_1)
            for i, name_1 in enumerate(names_1):
                for j, name_2 in enumerate(names_2):
                    if name_1 == name_2:
                        equal[i] = True
                        break

            if sum(equal) == len(names_1):
                return True

            else:
                return False


class NRTL_Yuan2019:
    def __init__(self):
        self.binary_names = []

        self.alphas = []
        self.a = []
        self.b = []
        self.e = []
        self.f = []
        self.tau = []

        self.add_pair("octanol", "nonanol", alphas=np.array([0.3, 0.3]), a=np.array([5.46, -0.66]),
                      b=np.array([-1329.55, -158.88]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("octanol", "water", alphas=np.array([0.26, 0.26]), a=np.array([-1.45, 5.92]),
                      b=np.array([-3011.98, -1581.72]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("nonanol", "water", alphas=np.array([0.2, 0.2]), a=np.array([-0.27, 10.30]),
                      b=np.array([-2326.25, -1670.40]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("dodecane", "water", alphas=np.array([0.2, 0.2]), a=np.array([-6.09, 23.43]),
                      b=np.array([-248.17, -4478.56]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("tridecane", "water", alphas=np.array([0.2, 0.2]), a=np.array([-6.91, 24.64]),
                      b=np.array([346.42, -4906.39]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("ethanol", "water", alphas=np.array([0.3, 0.3]), a=np.array([-0.80, 3.46]),
                      b=np.array([-3837.76, -884.89]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("octanol", "dodecane", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-306.76, -428.98]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("octanol", "ethanol", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-804.11, -342.30]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("dodecane", "ethanol", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([439.59, 1047.59]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("octanol", "tridecane", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-202.10, -614.853]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("tridecane", "ethanol", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([478.69, 1133.72]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("nonanol", "dodecane", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-444.49, -408.667]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("nonanol", "ethanol", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-685.72, -707.19]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("nonanol", "tridecane", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-410.03, -561.45]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))
        self.add_pair("dodecane", "tridecane", alphas=np.array([0.3, 0.3]), a=np.array([0, 0]),
                      b=np.array([-48.69, 50.33]), e=np.zeros(2), f=np.zeros(2), tau=np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.alphas[ind], self.a[ind], self.b[ind], self.e[ind], self.f[ind], self.tau[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        alphas = np.zeros((len(names), len(names)))
        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))
        f = np.zeros((len(names), len(names)))
        tau = np.zeros((len(names), len(names)))

        matrices = [alphas, a, b, e, f, tau]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return matrices

    def add_pair(self, name_1, name_2, alphas, a, b, e, f, tau):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.alphas.append(alphas)
        self.a.append(a)
        self.b.append(b)
        self.e.append(e)
        self.f.append(f)
        self.tau.append(tau)


class Experimental_Yuan2019:
    def __init__(self):
        self.names = []
        self.feeds = []
        self.corresponding_phases = []
        self.temperatures = []

        temp = [293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15]
        names = [
            ["ethanol", "octanol", "dodecane", "water"],
            ["ethanol", "octanol", "dodecane", "water"],
            ["ethanol", "octanol", "dodecane", "water"],
            ["ethanol", "octanol", "tridecane", "water"],
            ["ethanol", "octanol", "tridecane", "water"],
            ["ethanol", "octanol", "tridecane", "water"],
            ["ethanol", "nonanol", "dodecane", "water"],
            ["ethanol", "nonanol", "dodecane", "water"],
            ["ethanol", "nonanol", "dodecane", "water"],
            ["ethanol", "nonanol", "tridecane", "water"],
            ["ethanol", "nonanol", "tridecane", "water"],
            ["ethanol", "nonanol", "tridecane", "water"],
            ["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"],
            ["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"],
            ["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"]
        ]
        molar_fractions_given = False
        molar_masses = [
            np.array([46.07, 130.23, 170.33, 18.01528]),
            np.array([46.07, 130.23, 170.33, 18.01528]),
            np.array([46.07, 130.23, 170.33, 18.01528]),
            np.array([46.07, 130.23, 184.37, 18.01528]),
            np.array([46.07, 130.23, 184.37, 18.01528]),
            np.array([46.07, 130.23, 184.37, 18.01528]),
            np.array([46.07, 144.26, 170.33, 18.01528]),
            np.array([46.07, 144.26, 170.33, 18.01528]),
            np.array([46.07, 144.26, 170.33, 18.01528]),
            np.array([46.07, 144.26, 184.37, 18.01528]),
            np.array([46.07, 144.26, 184.37, 18.01528]),
            np.array([46.07, 144.26, 184.37, 18.01528]),
            np.array([46.07, 130.23, 144.26, 170.33, 184.37, 18.01528]),
            np.array([46.07, 130.23, 144.26, 170.33, 184.37, 18.01528]),
            np.array([46.07, 130.23, 144.26, 170.33, 184.37, 18.01528])
        ]
        sheet_names = ["Table_2_1", "Table_2_2", "Table_2_3",
                       "Table_3_1", "Table_3_2", "Table_3_3",
                       "Table_4_1", "Table_4_2", "Table_4_3",
                       "Table_5_1", "Table_5_2", "Table_5_3",
                       "Table_7_1", "Table_7_2", "Table_7_3"]
        excel_path = os.path.join(os.getcwd(), "source_property_data", "yuan_2019", "to_read.xlsx")

        for j in range(len(names)):
            data = pd.read_excel(excel_path, sheet_name=sheet_names[j], header=None, index_col=None)
            self.add_system(temp[j], names[j], molar_fractions_given, molar_masses[j], data)

    def add_system(self, temp, names, molar_fractions_given, molar_masses, data):
        matrix_to_read = data.to_numpy()
        num_comp = len(names)
        phase_1 = np.zeros((len(matrix_to_read), num_comp))
        phase_2 = np.zeros((len(matrix_to_read), num_comp))

        for i in range(len(matrix_to_read)):
            phase_1[i] = matrix_to_read[i][:num_comp]
            phase_2[i] = matrix_to_read[i][num_comp:2 * num_comp]

        if not molar_fractions_given:
            for i in range(len(phase_1)):
                phase_1[i] = convert_mass_fr_to_molar_fr(phase_1[i], molar_masses)
                phase_2[i] = convert_mass_fr_to_molar_fr(phase_2[i], molar_masses)

        corresponding_phases = []
        for i in range(len(phase_1)):
            corresponding_phases.append([phase_1[i], phase_2[i]])

        feeds = []
        for i in range(len(corresponding_phases)):
            feed = np.zeros(len(names))
            for phase in corresponding_phases[i]:
                feed = feed + phase

            feed = feed / sum(feed)
            feeds.append(feed)

        self.names.append(names)
        self.feeds.append(feeds)
        self.corresponding_phases.append(corresponding_phases)
        self.temperatures.append(temp)

    @staticmethod
    def compare_names(names_1, names_2):
        if len(names_1) != len(names_2):
            return False

        else:
            equal = [False] * len(names_1)
            for i, name_1 in enumerate(names_1):
                for j, name_2 in enumerate(names_2):
                    if name_1 == name_2:
                        equal[i] = True
                        break

            if sum(equal) == len(names_1):
                return True

            else:
                return False


class UNIQUAC_Pure_Yuan2018:
    def __init__(self):
        self.names = ["ethanol", "hexanol", "heptanol", "decane", "undecane", "water"]

        self.q = np.array([1.97, 4.13, 4.67, 6.02, 6.56, 1.40])
        self.q_tilde = np.array([1.97, 4.13, 4.67, 6.02, 6.56, 1.40])
        self.r = np.array([2.11, 4.80, 5.48, 7.20, 7.87, 0.92])

    def get_pure_para(self, name):
        for ind, name_2 in enumerate(self.names):
            if str(name) == str(name_2):
                return [self.q[ind], self.q_tilde[ind], self.r[ind]]

        return None


class UNIQUAC_Yuan2018:
    def __init__(self):
        self.binary_names = []

        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.e = []

        self.add_pair("ethanol", "hexanol", a=np.array([0.774, -1.173]), b=np.array([-325.241, 816.280]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "heptanol", a=np.array([0, 0]), b=np.array([-121.942, 385.105]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "decane", a=np.array([0.133, -0.645]), b=np.array([-26.0029, -188.648]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "undecane", a=np.array([0, 0]), b=np.array([-57.4879, -305.506]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "water", a=np.array([2.005, -2.494]), b=np.array([-684.446, 915.981]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("hexanol", "heptanol", a=np.array([0, 0]), b=np.array([-378.437, 227.650]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("hexanol", "decane", a=np.array([0, 0]), b=np.array([209.283, -211.834]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("hexanol", "undecane", a=np.array([0, 0]), b=np.array([268.333, -330.156]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("hexanol", "water", a=np.array([-0.283, -0.184]), b=np.array([27.309, 29.930]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("heptanol", "decane", a=np.array([0, 0]), b=np.array([269.119, -364.908]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("heptanol", "undecane", a=np.array([0, 0]), b=np.array([302.893, -467.431]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("heptanol", "water", a=np.array([0, 0]), b=np.array([-70.5651, -18.8341]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decane", "undecane", a=np.array([0, 0]), b=np.array([54.5282, -60.2206]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decane", "water", a=np.array([0, 0]), b=np.array([-1216.844, -149.114]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("undecane", "water", a=np.array([4.965, -2.579]), b=np.array([-2698.77, 663.267]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        q_list = np.zeros(len(names))
        q_tilde_list = np.zeros(len(names))
        r_list = np.zeros(len(names))

        pure_storage = UNIQUAC_Pure_Yuan2018()

        for i, name in enumerate(names):
            paras = pure_storage.get_pure_para(name)
            if paras is None:
                return None

            else:
                q_list[i] = paras[0]
                q_tilde_list[i] = paras[1]
                r_list[i] = paras[2]

        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        c = np.zeros((len(names), len(names)))
        d = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))

        matrices = [a, b, c, d, e]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return [q_list, q_tilde_list, r_list] + matrices

    def add_pair(self, name_1, name_2, a, b, c, d, e):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.a.append(a)
        self.b.append(b)
        self.c.append(c)
        self.d.append(d)
        self.e.append(e)


class Experimental_Yuan2018:
    def __init__(self):
        self.names = []
        self.feeds = []
        self.corresponding_phases = []
        self.temperatures = []

        temp = 298.15
        names = [
            ["ethanol", "hexanol", "decane", "water"],
            ["ethanol", "hexanol", "undecane", "water"],
            ["ethanol", "heptanol", "decane", "water"],
            ["ethanol", "heptanol", "undecane", "water"],
            ["ethanol", "hexanol", "heptanol", "decane", "undecane", "water"]
        ]
        molar_fractions_given = False
        molar_masses = [
            np.array([46.07, 102.162, 142.29, 18.01528]),
            np.array([46.07, 102.162, 156.31, 18.01528]),
            np.array([46.07, 116.2, 142.29, 18.01528]),
            np.array([46.07, 116.2, 156.31, 18.01528]),
            np.array([46.07, 102.162, 116.2, 142.29, 156.31, 18.01528])
        ]
        sheet_names = ["Table_3", "Table_4", "Table_5",
                       "Table_6", "Table_11"]
        excel_path = os.path.join(os.getcwd(), "source_property_data", "yuan_2018", "to_read.xlsx")

        for j in range(len(names)):
            data = pd.read_excel(excel_path, sheet_name=sheet_names[j], header=None, index_col=None)
            self.add_system(temp, names[j], molar_fractions_given, molar_masses[j], data)

    def add_system(self, temp, names, molar_fractions_given, molar_masses, data):
        matrix_to_read = data.to_numpy()
        num_comp = len(names)
        phase_1 = np.zeros((len(matrix_to_read), num_comp))
        phase_2 = np.zeros((len(matrix_to_read), num_comp))

        for i in range(len(matrix_to_read)):
            phase_1[i] = matrix_to_read[i][:num_comp]
            phase_2[i] = matrix_to_read[i][num_comp:2 * num_comp]

        if not molar_fractions_given:
            for i in range(len(phase_1)):
                phase_1[i] = convert_mass_fr_to_molar_fr(phase_1[i], molar_masses)
                phase_2[i] = convert_mass_fr_to_molar_fr(phase_2[i], molar_masses)

        corresponding_phases = []
        for i in range(len(phase_1)):
            corresponding_phases.append([phase_1[i], phase_2[i]])

        feeds = []
        for i in range(len(corresponding_phases)):
            feed = np.zeros(len(names))
            for phase in corresponding_phases[i]:
                feed = feed + phase

            feed = feed / sum(feed)
            feeds.append(feed)

        self.names.append(names)
        self.feeds.append(feeds)
        self.corresponding_phases.append(corresponding_phases)
        self.temperatures.append(temp)

    @staticmethod
    def compare_names(names_1, names_2):
        if len(names_1) != len(names_2):
            return False

        else:
            equal = [False] * len(names_1)
            for i, name_1 in enumerate(names_1):
                for j, name_2 in enumerate(names_2):
                    if name_1 == name_2:
                        equal[i] = True
                        break

            if sum(equal) == len(names_1):
                return True

            else:
                return False


class UNIQUAC_Pure_Yuan2020:
    def __init__(self):
        self.names = ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water", "octadecane"]

        self.q = np.array([1.97, 6.29, 6.83, 8.18, 8.72, 1.40, 15.26])
        self.q_tilde = np.array([1.97, 6.29, 6.83, 8.18, 8.72, 1.40, 15.26])
        self.r = np.array([2.10, 7.50, 8.17, 9.89, 10.57, 0.92, 16.22])

    def get_pure_para(self, name):
        for ind, name_2 in enumerate(self.names):
            if str(name) == str(name_2):
                return [self.q[ind], self.q_tilde[ind], self.r[ind]]

        return None


class UNIQUAC_Yuan2020:
    def __init__(self):
        self.binary_names = []

        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.e = []

        self.add_pair("ethanol", "decanol", a=np.array([0, 0]), b=np.array([-357.84, 250.89]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "undecanol", a=np.array([0, 0]), b=np.array([-376.10, 256.25]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "tetradecane", a=np.array([0, 0]), b=np.array([-9.12, -439.29]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "pentadecane", a=np.array([0, 0]), b=np.array([-2.17, -461.43]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "water", a=np.array([2.00, -2.49]), b=np.array([-313.05, -5730.88]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decanol", "tetradecane", a=np.array([0, 0]), b=np.array([156.86, -216.72]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decanol", "pentadecane", a=np.array([0, 0]), b=np.array([191.48, -280.23]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decanol", "water", a=np.array([0, 0]), b=np.array([-246.41, -564.17]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("undecanol", "tetradecane", a=np.array([0, 0]), b=np.array([171.09, -229.20]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("undecanol", "pentadecane", a=np.array([0, 0]), b=np.array([180.48, -249.00]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("undecanol", "water", a=np.array([-1.06, -0.01]), b=np.array([16.22, -608.69]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("tetradecane", "water", a=np.array([5.59, -2.37]), b=np.array([-4338.48, 480.59]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("pentadecane", "water", a=np.array([0, 0]), b=np.array([-11367.90, -240.30]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decanol", "undecanol", a=np.array([0, 0]), b=np.array([53.12, -58.75]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("tetradecane", "pentadecane", a=np.array([0, 0]), b=np.array([54.53, -60.22]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("ethanol", "octadecane", a=np.array([0.97, 3.82]), b=np.array([-218.61, -1812.48]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("decanol", "octadecane", a=np.array([0, 0]), b=np.array([277.56, -474.97]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("undecanol", "octadecane", a=np.array([0, 0]), b=np.array([274.68, -465.18]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("tetradecane", "octadecane", a=np.array([0, 0]), b=np.array([130.53, -149.61]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("pentadecane", "octadecane", a=np.array([0, 0]), b=np.array([127.10, -142.93]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))
        self.add_pair("water", "octadecane", a=np.array([0, 0]), b=np.array([-422.24, -1668.58]), c=np.zeros(2),
                      d=np.zeros(2), e=np.zeros(2))

    def get_binary_interaction_para(self, name_1, name_2):
        for ind, name in enumerate(self.binary_names):
            if str(name_1 + "--" + name_2) == str(name):
                return [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]

            elif str(name_2 + "--" + name_1) == str(name):
                ret = [self.a[ind], self.b[ind], self.c[ind], self.d[ind], self.e[ind]]
                return [np.flip(el) for el in ret]

        return None

    def combine_interactions(self, names):
        q_list = np.zeros(len(names))
        q_tilde_list = np.zeros(len(names))
        r_list = np.zeros(len(names))

        pure_storage = UNIQUAC_Pure_Yuan2020()

        for i, name in enumerate(names):
            paras = pure_storage.get_pure_para(name)
            if paras is None:
                return None

            else:
                q_list[i] = paras[0]
                q_tilde_list[i] = paras[1]
                r_list[i] = paras[2]

        a = np.zeros((len(names), len(names)))
        b = np.zeros((len(names), len(names)))
        c = np.zeros((len(names), len(names)))
        d = np.zeros((len(names), len(names)))
        e = np.zeros((len(names), len(names)))

        matrices = [a, b, c, d, e]

        todo = itertools.combinations(list(range(len(names))), 2)
        for combination in todo:
            paras = self.get_binary_interaction_para(names[combination[0]], names[combination[1]])
            if paras is None:
                return None

            else:
                for j in range(len(matrices)):
                    matrices[j][combination[0]][combination[1]] = paras[j][0]
                    matrices[j][combination[1]][combination[0]] = paras[j][1]

        return [q_list, q_tilde_list, r_list] + matrices

    def add_pair(self, name_1, name_2, a, b, c, d, e):
        self.binary_names.append(str(name_1 + "--" + name_2))
        self.a.append(a)
        self.b.append(b)
        self.c.append(c)
        self.d.append(d)
        self.e.append(e)


class Experimental_Yuan2020:
    def __init__(self):
        self.names = []
        self.feeds = []
        self.corresponding_phases = []
        self.temperatures = []

        temp = [293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15,
                293.15, 298.15, 303.15]
        names = [
            ["ethanol", "decanol", "tetradecane", "water"],
            ["ethanol", "decanol", "tetradecane", "water"],
            ["ethanol", "decanol", "tetradecane", "water"],
            ["ethanol", "decanol", "pentadecane", "water"],
            ["ethanol", "decanol", "pentadecane", "water"],
            ["ethanol", "decanol", "pentadecane", "water"],
            ["ethanol", "undecanol", "tetradecane", "water"],
            ["ethanol", "undecanol", "tetradecane", "water"],
            ["ethanol", "undecanol", "tetradecane", "water"],
            ["ethanol", "undecanol", "pentadecane", "water"],
            ["ethanol", "undecanol", "pentadecane", "water"],
            ["ethanol", "undecanol", "pentadecane", "water"],
            ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"],
            ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"],
            ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"]
        ]
        molar_fractions_given = False
        molar_masses = [
            np.array([46.07, 158.28, 198.39, 18.01528]),
            np.array([46.07, 158.28, 198.39, 18.01528]),
            np.array([46.07, 158.28, 198.39, 18.01528]),
            np.array([46.07, 158.28, 212.42, 18.01528]),
            np.array([46.07, 158.28, 212.42, 18.01528]),
            np.array([46.07, 158.28, 212.42, 18.01528]),
            np.array([46.07, 172.31, 198.39, 18.01528]),
            np.array([46.07, 172.31, 198.39, 18.01528]),
            np.array([46.07, 172.31, 198.39, 18.01528]),
            np.array([46.07, 172.31, 212.42, 18.01528]),
            np.array([46.07, 172.31, 212.42, 18.01528]),
            np.array([46.07, 172.31, 212.42, 18.01528]),
            np.array([46.07, 158.28, 172.31, 198.39, 212.42, 18.01528]),
            np.array([46.07, 158.28, 172.31, 198.39, 212.42, 18.01528]),
            np.array([46.07, 158.28, 172.31, 198.39, 212.42, 18.01528])
        ]
        sheet_names = ["Table_2_1", "Table_2_2", "Table_2_3",
                       "Table_3_1", "Table_3_2", "Table_3_3",
                       "Table_4_1", "Table_4_2", "Table_4_3",
                       "Table_5_1", "Table_5_2", "Table_5_3",
                       "Table_10_1", "Table_10_2", "Table_10_3"]
        excel_path = os.path.join(os.getcwd(), "source_property_data", "yuan_2020", "to_read.xlsx")

        for j in range(len(names)):
            data = pd.read_excel(excel_path, sheet_name=sheet_names[j], header=None, index_col=None)
            self.add_system(temp[j], names[j], molar_fractions_given, molar_masses[j], data)

    def add_system(self, temp, names, molar_fractions_given, molar_masses, data):
        matrix_to_read = data.to_numpy()
        num_comp = len(names)
        phase_1 = np.zeros((len(matrix_to_read), num_comp))
        phase_2 = np.zeros((len(matrix_to_read), num_comp))

        for i in range(len(matrix_to_read)):
            phase_1[i] = matrix_to_read[i][:num_comp]
            phase_2[i] = matrix_to_read[i][num_comp:2 * num_comp]

        if not molar_fractions_given:
            for i in range(len(phase_1)):
                phase_1[i] = convert_mass_fr_to_molar_fr(phase_1[i], molar_masses)
                phase_2[i] = convert_mass_fr_to_molar_fr(phase_2[i], molar_masses)

        corresponding_phases = []
        for i in range(len(phase_1)):
            corresponding_phases.append([phase_1[i], phase_2[i]])

        feeds = []
        for i in range(len(corresponding_phases)):
            feed = np.zeros(len(names))
            for phase in corresponding_phases[i]:
                feed = feed + phase

            feed = feed / sum(feed)
            feeds.append(feed)

        self.names.append(names)
        self.feeds.append(feeds)
        self.corresponding_phases.append(corresponding_phases)
        self.temperatures.append(temp)

    @staticmethod
    def compare_names(names_1, names_2):
        if len(names_1) != len(names_2):
            return False

        else:
            equal = [False] * len(names_1)
            for i, name_1 in enumerate(names_1):
                for j, name_2 in enumerate(names_2):
                    if name_1 == name_2:
                        equal[i] = True
                        break

            if sum(equal) == len(names_1):
                return True

            else:
                return False
