
import os

import property_data
import point_discretization


if __name__ == "__main__":
    # create point discretizations if necessary, list of lists, first el is num_comp, second rec_steps
    discretizations_todo = [[2, 8], [3, 7], [4, 6], [5, 5], [6, 4]]

    # dir for disc data
    discretization_path = os.path.join(os.getcwd(), "data", "discretization")
    if not os.path.isdir(discretization_path):
        os.mkdir(discretization_path)

    for todo_el in discretizations_todo:
        filename = os.path.join(discretization_path, str(todo_el[0]) + "_" + str(todo_el[1]))

        # check if not already existing
        if not os.path.isdir(filename):
            point_discretization.PointDisc(num_comp=todo_el[0], recursion_steps=todo_el[1],
                                           load=False, store=True)

    # specify ray actors
    actors_parallelized = 0

    # dir for results
    results_path = os.path.join(os.getcwd(), "data", "lle_results")
    if not os.path.isdir(results_path):
        os.mkdir(results_path)

    """# ryll2009_1
    names_list = [["nitromethane", "1_hexanol"], ["benzene", "ethanol", "water"],
                  ["1_hexanol", "nitromethane", "water"], ["ethyl_formate", "ethylene_glycol", "water"],
                  ["water", "isobutanol", "furfural"], ["water", "acetone", "phenol"]]
    store_lle = True
    temperature_kelvin = [273.15 + 25, 273.15 + 25, 273.15 + 21, 273.15 + 30, 273.15 + 25, 273.15 + 56.5]
    deciding_index = None
    general_name = "Ryll2009_1"
    ge_model_name = "UNIQUAC"
    tau_given = None
    experimental_data_given = False

    property_data.InitClass(names_list=names_list[0], store_lle=store_lle,
                            point_discretization_rec_steps=8,
                            name_for_data=general_name, temperature_kelvin=temperature_kelvin[0],
                            deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                            ge_model_name=ge_model_name, tau_given=tau_given,
                            experimental_data_given=experimental_data_given)

    for j in [7]:
        for i in range(1, len(names_list)):
            property_data.InitClass(names_list=names_list[i], store_lle=store_lle,
                                    point_discretization_rec_steps=j,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin[i],
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # ryll2009_2
    names_list = ["water", "nitromethane", "nonanol"]
    store_lle = True
    temperature_kelvin = 273.15 + 45.1
    deciding_index = None
    general_name = "Ryll2009_2"
    ge_model_name = "UNIQUAC"
    tau_given = None
    experimental_data_given = False

    for j in [7]:
        property_data.InitClass(names_list=names_list, store_lle=store_lle,
                                point_discretization_rec_steps=j,
                                name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                ge_model_name=ge_model_name, tau_given=tau_given,
                                experimental_data_given=experimental_data_given)"""

    # chen2000
    names_list = [["n_hexane", "benzene", "sulfolane"], ["n_hexane", "toluene", "sulfolane"],
                  ["n_hexane", "xylene", "sulfolane"], ["n_octane", "benzene", "sulfolane"],
                  ["n_octane", "toluene", "sulfolane"], ["n_octane", "xylene", "sulfolane"]]
    store_lle = True
    temperature_kelvin = 298.15
    deciding_index = 0
    general_name = "Chen2000_1"
    ge_model_name = "NRTL"
    tau_given = True
    experimental_data_given = True

    for j in range(len(names_list)):
        for i in [7]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    """# 4 comp systems
    names_list = [["n_hexane", "n_octane", "benzene", "sulfolane"],
                  ["n_hexane", "benzene", "xylene", "sulfolane"],
                  ["n_octane", "toluene", "xylene", "sulfolane"]]

    for j in range(len(names_list)):
        for i in [6]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # 5 comp system
    names_list = ["n_hexane", "n_octane", "benzene", "toluene", "sulfolane"]
    property_data.InitClass(names_list=names_list, store_lle=store_lle,
                            point_discretization_rec_steps=5,
                            name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                            deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                            ge_model_name=ge_model_name, tau_given=tau_given,
                            experimental_data_given=experimental_data_given)

    # chen2001
    store_lle = True
    temperature_kelvin = 298.15
    deciding_index = 0
    general_name = "Chen2001"
    ge_model_name = "NRTL"
    tau_given = True
    experimental_data_given = True

    # 4 comp systems
    names_list = [["hexane", "heptane", "toluene", "sulfolane"],
                  ["heptane", "octane", "m_xylene", "sulfolane"],
                  ["heptane", "benzene", "toluene", "sulfolane"]]

    for j in range(len(names_list)):
        for i in [6]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # 5 comp systems
    names_list = [["hexane", "heptane", "toluene", "m_xylene", "sulfolane"],
                  ["heptane", "octane", "benzene", "m_xylene", "sulfolane"]]

    for j in range(len(names_list)):
        for i in [5]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # yuan2019
    store_lle = True
    temperature_kelvin = [293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15]
    deciding_index = 3
    general_name = "Yuan2019"
    ge_model_name = "NRTL"
    tau_given = False
    experimental_data_given = True

    # 4 comp systems
    names_list = [
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
            ["ethanol", "nonanol", "tridecane", "water"]]

    for j in range(len(names_list)):
        for i in [6]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin[j],
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # 6 comp systems
    names_list = [["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"],
                  ["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"],
                  ["ethanol", "octanol", "nonanol", "dodecane", "tridecane", "water"]]
    temperature_kelvin = [293.15, 298.15, 303.15]
    deciding_index = 5

    for j in range(len(names_list)):
        property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                point_discretization_rec_steps=4,
                                name_for_data=general_name, temperature_kelvin=temperature_kelvin[j],
                                deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                ge_model_name=ge_model_name, tau_given=tau_given,
                                experimental_data_given=experimental_data_given)

    # yuan2018
    store_lle = True
    temperature_kelvin = 298.15
    deciding_index = 3
    general_name = "Yuan2018"
    ge_model_name = "UNIQUAC"
    tau_given = False
    experimental_data_given = True

    # 4 comp systems
    names_list = [
        ["ethanol", "hexanol", "decane", "water"],
        ["ethanol", "hexanol", "undecane", "water"],
        ["ethanol", "heptanol", "decane", "water"],
        ["ethanol", "heptanol", "undecane", "water"]
        ]

    for j in range(len(names_list)):
        for i in [6]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # 6 comp system
    names_list = ["ethanol", "hexanol", "heptanol", "decane", "undecane", "water"]
    temperature_kelvin = 298.15
    deciding_index = 5

    property_data.InitClass(names_list=names_list, store_lle=store_lle,
                            point_discretization_rec_steps=4,
                            name_for_data=general_name, temperature_kelvin=temperature_kelvin,
                            deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                            ge_model_name=ge_model_name, tau_given=tau_given,
                            experimental_data_given=experimental_data_given)

    # yuan2020
    store_lle = True
    temperature_kelvin = [293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15,
                          293.15, 298.15, 303.15]

    deciding_index = 3
    general_name = "Yuan2020"
    ge_model_name = "UNIQUAC"
    tau_given = False
    experimental_data_given = True

    # 4 comp systems
    names_list = [
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
        ["ethanol", "undecanol", "pentadecane", "water"]
    ]

    for j in range(len(names_list)):
        for i in [6]:
            property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                    point_discretization_rec_steps=i,
                                    name_for_data=general_name, temperature_kelvin=temperature_kelvin[j],
                                    deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                    ge_model_name=ge_model_name, tau_given=tau_given,
                                    experimental_data_given=experimental_data_given)

    # 6 comp systems
    names_list = [
        ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"],
        ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"],
        ["ethanol", "decanol", "undecanol", "tetradecane", "pentadecane", "water"]
    ]
    temperature_kelvin = [293.15, 298.15, 303.15]
    deciding_index = 5

    for j in range(len(names_list)):
        property_data.InitClass(names_list=names_list[j], store_lle=store_lle,
                                point_discretization_rec_steps=4,
                                name_for_data=general_name, temperature_kelvin=temperature_kelvin[j],
                                deciding_index=deciding_index, actors_parallelized=actors_parallelized,
                                ge_model_name=ge_model_name, tau_given=tau_given,
                                experimental_data_given=experimental_data_given)"""
