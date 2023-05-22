
import time
import os
import numpy as np

import point_discretization
import plotter
import lle


class Analysis:
    """
    class which takes property data and several information as input and generates results for the
    paper
    """
    def __init__(self, num_comp, names_list, point_discretization_rec_steps, ge_model,
                 store_lle, experimental_data, temperature_kelvin,
                 deciding_index, actors_para, general_name):
        # create a general folder for this source if not existing
        self.general_path = os.path.join(os.getcwd(), "data", "lle_results", general_name)
        if not os.path.isdir(self.general_path):
            os.mkdir(self.general_path)

        self.num_comp = num_comp
        self.comp_names = names_list
        self.temperature = temperature_kelvin

        # simpl disc class, just loaded
        self.point_disc = point_discretization.PointDisc(num_comp=num_comp,
                                                         recursion_steps=point_discretization_rec_steps,
                                                         load=True, store=False)

        # init plotter, obviously, we can only plot something for 4 or less components:
        if num_comp <= 4:
            self.plotter_instance = plotter.Plotter(num_comp=num_comp)

        # initialize the ge model
        self.ge_model = ge_model
        print("lle analysis for", self.comp_names)

        name = ""
        for i in range(len(self.comp_names) - 1):
            name = name + self.comp_names[i] + "--"

        name = name + self.comp_names[-1]

        path = os.path.join(self.general_path, str(self.num_comp) + "_" + str(
            point_discretization_rec_steps) + "___" + name + "_temp_" +
                            str(self.temperature))

        # only create if not existing
        if store_lle and not os.path.isdir(path):
            # lle analysis
            start_time = time.time()
            lle_analysis = lle.miscibility_analysis(self.point_disc, self.ge_model, self.temperature,
                                                    construct=store_lle, path=None, actors_for_para=actors_para)
            construction_time = time.time() - start_time
            print("\ntime", construction_time)

            path = lle_analysis.store_phase_eq_liquid(name, path=self.general_path)

            file = open(os.path.join(path, "construction_time.txt"), "w+")
            file.write(str(construction_time))
            file.close()

            file = open(os.path.join(path, "len_points.txt"), "w+")
            file.write(str(len(self.point_disc.points_mfr)))
            file.close()

            if num_comp <= 4:
                self.plotter_instance.plot_outer_simplex("black")
                for simplex in lle_analysis.miscibility_gap_simplices:
                    self.plotter_instance.plot_misc_gap_simplex(simplex, "red", only_hetero_edges=False,
                                                                plot_surface=True)

                self.plotter_instance.save_plot(os.path.join(path, name + "_lle.png"))

        self.loaded_lle = lle.miscibility_analysis(self.point_disc, self.ge_model, self.temperature,
                                                   construct=False, path=path, actors_for_para=actors_para)

        if experimental_data is not None:
            # search for correct index
            index = None
            for i in range(len(experimental_data.names)):
                if experimental_data.compare_names(self.comp_names, experimental_data.names[i]) and \
                        np.abs(experimental_data.temperatures[i] - self.temperature) < 0.01:
                    index = i
                    break

            sum_squared_differences = 0
            sum_abs_differences = 0
            missing_indices = []
            for j, feed in enumerate(experimental_data.feeds[index]):
                phases_flowrates, relevant_simplex = self.loaded_lle.find_phase_split(feed, relevant_simplex=None,
                                                discretized_system=self.point_disc,
                                                miscibility_gap_simplices=self.loaded_lle.miscibility_gap_simplices,
                                                num_comp=num_comp)
                phases_mfr = [el / sum(el) for el in phases_flowrates]

                if len(phases_mfr) != len(experimental_data.corresponding_phases[index][j]):
                    missing_indices.append([j, len(phases_mfr), len(experimental_data.corresponding_phases[index][j])])
                    print("missing feed")
                    print("simu", len(phases_mfr))
                    print("exp", len(experimental_data.corresponding_phases[index][j]))
                    print("\n\n")

                else:
                    # depending on deciding index, enumerate the phases and compare reduced and experimental
                    first_index_reduced = np.argmax([phases_mfr[0][deciding_index],
                                                     phases_mfr[1][deciding_index]])
                    first_index_exp = np.argmax([experimental_data.corresponding_phases[index][j][0][deciding_index],
                                                experimental_data.corresponding_phases[index][j][1][deciding_index]])

                    sum_squared_differences = sum_squared_differences + np.sum(np.square(
                        phases_mfr[first_index_reduced] -
                        experimental_data.corresponding_phases[index][j][first_index_exp]))
                    sum_squared_differences = sum_squared_differences + np.sum(np.square(
                        phases_mfr[first_index_reduced - 1] -
                        experimental_data.corresponding_phases[index][j][first_index_exp - 1]))

                    sum_abs_differences = sum_abs_differences + np.sum(np.abs(
                        phases_mfr[first_index_reduced] -
                        experimental_data.corresponding_phases[index][j][first_index_exp]))
                    sum_abs_differences = sum_abs_differences + np.sum(np.abs(
                        phases_mfr[first_index_reduced - 1] -
                        experimental_data.corresponding_phases[index][j][first_index_exp - 1]))

            if len(experimental_data.feeds[index]) > len(missing_indices):
                # we divide by the number of phases, 2, and the number of feeds and the number of components as chen
                sum_squared_differences = np.sqrt(sum_squared_differences /
                    (self.num_comp * 2 * (len(experimental_data.feeds[index]) - len(missing_indices))))
                sum_abs_differences = sum_abs_differences / (
                    self.num_comp * 2 * (len(experimental_data.feeds[index]) - len(missing_indices)))

            print("\n\nfinal RMSD", sum_squared_differences)
            print("final abs diff:", sum_abs_differences)
            file = open(os.path.join(path, "comparison_simu_exp.txt"), "w+")
            file.write("root mean squared error: " + str(sum_squared_differences) + "\n")
            file.write("abs diff: " + str(sum_abs_differences) + "\n")
            file.write("number of feeds: " + str(len(experimental_data.feeds[index])) + "\n")
            file.write("missing feeds: " + str(len(missing_indices)) + "\n\n")
            for m_f in missing_indices:
                file.write("index " + str(m_f[0]) + " simu " + str(m_f[1]) + " exp " + str(m_f[2]) + "\n")

            file.close()
