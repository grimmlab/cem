
import numpy as np
import matplotlib.pyplot as plt

import point_discretization


class Plotter:
    def __init__(self, num_comp):
        self.num_comp = num_comp

        # to be able to plot the simplex for the whole component system, we just construct it
        # in the same way as in the simplex discretization
        self.vertices_outer_simplex = []
        for i in range(self.num_comp-1):
            basis_vector = np.zeros(self.num_comp-1)
            basis_vector[i] = 1 / np.sqrt(2)
            self.vertices_outer_simplex.append(basis_vector)

        # the last point
        self.vertices_outer_simplex.append(np.ones(self.num_comp-1) * (1 + np.sqrt(self.num_comp)) / (
                (self.num_comp - 1) * np.sqrt(2)))

        if self.num_comp == 3:
            # rotation with psi
            psi = 2 * np.pi * 285 / 360
            rotation_matrix = np.array([[np.cos(psi), -1 * np.sin(psi)], [np.sin(psi), np.cos(psi)]])
            for i in range(len(self.vertices_outer_simplex)):
                self.vertices_outer_simplex[i] = np.matmul(rotation_matrix, self.vertices_outer_simplex[i])

        # the transformation matrices may be useful at some point
        self.matrix_mfr_to_cart, self.matrix_cart_to_mfr = point_discretization.PointDisc.get_basis_change(
            self.vertices_outer_simplex)

        if self.num_comp == 4:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')

    def plot_outer_simplex(self, color):
        if self.num_comp == 3:
            plt.gca().set_aspect('equal', adjustable='box')
            points = [self.vertices_outer_simplex[j] for j in range(self.num_comp)]
            for j in range(self.num_comp):
                plt.plot([points[j - 1][0], points[j][0]], [points[j - 1][1], points[j][1]], color=color, linewidth=1)

        elif self.num_comp == 4:
            points = [self.vertices_outer_simplex[j] for j in range(self.num_comp)]
            for i in range(self.num_comp):
                for j in range(i+1, self.num_comp):
                    self.ax.plot([points[i][0], points[j][0]], [points[i][1], points[j][1]],
                                 zs=[points[i][2], points[j][2]], color=color)

    def plot_misc_gap_simplex(self, simplex, color, only_hetero_edges=None, size=None, plot_surface=False):
        if self.num_comp == 2:
            # we just plot the points
            for i in range(self.num_comp):
                plt.plot(simplex.points_molar_fractions[i][0], 0, marker="o", markersize=size, color=color)

        elif self.num_comp == 3:
            # plot also the surface between points, except for 3 phase simplices
            if plot_surface and np.sum(simplex.edge_classification) < 6:
                t1 = plt.Polygon(simplex.points_coordinates_cart, color=color)
                plt.gca().add_patch(t1)

            else:
                # plot three phase simplex with transparent color
                if plot_surface and np.sum(simplex.edge_classification) == 6:
                    t1 = plt.Polygon(simplex.points_coordinates_cart, color=color, alpha=0.5)
                    plt.gca().add_patch(t1)

                else:

                    if only_hetero_edges is None:
                        for i in range(self.num_comp):
                            for j in range(i+1, self.num_comp):
                                plt.plot([simplex.points_coordinates_cart[i][0], simplex.points_coordinates_cart[j][0]],
                                         [simplex.points_coordinates_cart[i][1], simplex.points_coordinates_cart[j][1]],
                                         color=color, linewidth=1)

                    elif only_hetero_edges:
                        for i in range(self.num_comp):
                            for j in range(i+1, self.num_comp):
                                if simplex.edge_classification[i][j] > 0.5:
                                    plt.plot([simplex.points_coordinates_cart[i][0], simplex.points_coordinates_cart[j][0]],
                                             [simplex.points_coordinates_cart[i][1], simplex.points_coordinates_cart[j][1]],
                                             color=color, linewidth=1)

                    else:
                        for i in range(self.num_comp):
                            for j in range(i+1, self.num_comp):
                                if simplex.edge_classification[i][j] < 0.5:
                                    plt.plot([simplex.points_coordinates_cart[i][0], simplex.points_coordinates_cart[j][0]],
                                             [simplex.points_coordinates_cart[i][1], simplex.points_coordinates_cart[j][1]],
                                             color=color, linewidth=1)

        elif self.num_comp == 4:
            for i in range(self.num_comp):
                for j in range(i+1, self.num_comp):
                    self.ax.plot([simplex.points_coordinates_cart[i][0], simplex.points_coordinates_cart[j][0]],
                                 [simplex.points_coordinates_cart[i][1], simplex.points_coordinates_cart[j][1]],
                                 zs=[simplex.points_coordinates_cart[i][2], simplex.points_coordinates_cart[j][2]],
                                 color=color)

    def transform_molar_fr_to_cartesian(self, molar_fractions):
        """
        A * lambda = (1, p), we cut off the first entry
        """

        return np.matmul(self.matrix_mfr_to_cart, molar_fractions)[1:]

    def transform_cartesian_to_molar_fr(self, cartesian_point):
        """
        lambda = A_inv * (1, p)
        """
        vector = np.empty(self.num_comp)
        vector[0] = 1
        vector[1:] = cartesian_point

        return np.matmul(self.matrix_cart_to_mfr, vector)

    @staticmethod
    def save_plot(path):
        plt.savefig(path, dpi=400)
        plt.close()

    @staticmethod
    def show_plot():
        plt.show()
        plt.close()

    @staticmethod
    def turn_axis_off():
        plt.axis('off')
