
import numpy as np


class NRTL:
    def __init__(self, alpha_matrix, a_matrix, b_matrix, e_matrix, f_matrix, tau_matrix, tau_given,
                 different_tau=False):
        """
        Class stores the binary parameters of all components.

        Computes the activity coefficient of component index as described in Molecular
        Thermodynamics of Fluid-Phase Equilibria 1998. There, we have two parameters tau_ij and
        tau_ji (non symmetric) and alpha (symmetric).

        As the interaction parameters are sometimes provided in different form, we have multiple
        options for this. Normally, we use the Aspen standard:
        tau_ij = a_ij + b_ij/T + e_ij*logT + f_ij*T (T in K).

        If the tau_matrix is directly given, we assume that it is rescaled by division by the
        temperature in Kelvin:
        tau_ij = tau_matrix_ij / T

        The bool different_tau indicates, if tau is given but not rescaled with temperature.
        """
        self.alpha_matrix = alpha_matrix
        self.tau_given = tau_given
        self.a_matrix = a_matrix
        self.b_matrix = b_matrix
        self.e_matrix = e_matrix
        self.f_matrix = f_matrix
        self.tau_without_temperature = tau_matrix
        self.tau_different = different_tau

        self.num_comp = len(self.alpha_matrix)

    def get_binary_parameters(self, temperature):
        """
        Binary interactions depend on the temperature (provided in Kelvin) and are calculated here.
        """
        tau_matrix = np.empty((self.num_comp, self.num_comp))
        upper_G_matrix = np.empty((self.num_comp, self.num_comp))
        if self.tau_given:
            if self.tau_different:
                for i in range(self.num_comp):
                    for j in range(self.num_comp):
                        tau_matrix[i][j] = self.tau_without_temperature[i][j]
                        upper_G_matrix[i][j] = np.exp(-1 * self.alpha_matrix[i][j] * tau_matrix[i][j])

            else:
                for i in range(self.num_comp):
                    for j in range(self.num_comp):
                        tau_matrix[i][j] = self.tau_without_temperature[i][j] / temperature
                        upper_G_matrix[i][j] = np.exp(-1 * self.alpha_matrix[i][j] * tau_matrix[i][j])

        else:
            for i in range(self.num_comp):
                for j in range(self.num_comp):
                    tau_matrix[i][j] = self.a_matrix[i][j] + (self.b_matrix[i][j] / temperature) + (
                            self.e_matrix[i][j] * np.log(temperature)) + (self.f_matrix[i][j] * temperature)

                    upper_G_matrix[i][j] = np.exp(-1 * self.alpha_matrix[i][j] * tau_matrix[i][j])

        return tau_matrix, upper_G_matrix

    def compute_activity_coefficient(self, molar_fractions, index, temperature):
        """
        Index refers to the component, where the activity coefficient should be calculated.
        """
        # get interaction parameters to given temperature (unit K)
        tau_matrix, upper_G_matrix = self.get_binary_parameters(temperature)

        # we use the same formula as in the reference from the description (6-169)
        # it is split in several terms for the sake of readability
        numerator_1 = 0
        denominator_1 = 0
        for j in range(self.num_comp):
            numerator_1 = numerator_1 + (tau_matrix[j][index] * upper_G_matrix[j][index] *
                                         molar_fractions[j])

            denominator_1 = denominator_1 + (upper_G_matrix[j][index] * molar_fractions[j])

        summand_1 = numerator_1 / denominator_1

        summand_2 = 0
        for j in range(self.num_comp):
            denom_first_part_sum_2 = 0
            for l in range(self.num_comp):
                denom_first_part_sum_2 = denom_first_part_sum_2 + (
                        upper_G_matrix[l][j] * molar_fractions[l])

            first_part_summand_2 = molar_fractions[j] * upper_G_matrix[index][j] / denom_first_part_sum_2

            numerator_2 = 0
            for r in range(self.num_comp):
                numerator_2 = numerator_2 + (molar_fractions[r] * tau_matrix[r][j] *
                                             upper_G_matrix[r][j])

            second_part_summand_2 = tau_matrix[index][j] - (numerator_2 / denom_first_part_sum_2)
            summand_2 = summand_2 + (first_part_summand_2 * second_part_summand_2)

        activity_coefficient = np.exp(summand_1 + summand_2)

        return activity_coefficient


class UNIQUAC:
    def __init__(self, a, b, c, d, e, r, q, q_2):
        """
        Similar as NRTL above.

        Scaling with temperature in Kelvin:
        tau_ij = exp(a_ij + b_ij/T + c_ij*ln(T) + d_ij*T + e_ij/T^2)

        If one wants to use the simplified approach with q == q' (=q_2), just insert the same values...
        """
        # parameter matrices
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e

        # parameter vectors
        self.r = r
        self.q = q
        self.q_2 = q_2

        self.num_comp = len(self.r)

    def get_binary_parameters(self, temperature):
        tau_matrix = np.empty((self.num_comp, self.num_comp))
        for i in range(self.num_comp):
            for j in range(self.num_comp):
                tau_matrix[i][j] = np.exp(self.a[i][j] + (self.b[i][j] / temperature) + (self.c[i][j] * np.log(
                    temperature)) + (self.d[i][j] * temperature) + (self.e[i][j] / np.square(temperature)))

        return tau_matrix

    def compute_activity_coefficient(self, molar_fractions, index, temperature):
        """
        ln gamma_i = ln gamma_i^comb + ln gamma_i_red
        """
        tau_matrix = self.get_binary_parameters(temperature)

        # some helper variables
        z = 10  # coordination number
        theta_index = self.q[index] * molar_fractions[index] / sum([
            self.q[i] * molar_fractions[i] for i in range(self.num_comp)])

        # referring to theta'
        theta_2_vector = np.array([self.q_2[i] * molar_fractions[i] for i in range(self.num_comp)])
        theta_2_vector = theta_2_vector / sum(theta_2_vector)

        t_vector = np.empty(self.num_comp)
        for i in range(self.num_comp):
            t_vector[i] = 0
            for k in range(self.num_comp):
                t_vector[i] = t_vector[i] + (theta_2_vector[k] * tau_matrix[k][i])

        phi_index = self.r[index] * molar_fractions[index] / sum([
            self.r[i] * molar_fractions[i] for i in range(self.num_comp)])

        l_vector = np.array([((z / 2) * (self.r[i] - self.q[i])) + 1 - self.r[i] for i in range(self.num_comp)])

        ln_gamma_index_comb = np.log(phi_index / molar_fractions[index]) + ((z / 2) * self.q[index] * np.log(
            theta_index / phi_index)) + l_vector[index] - ((phi_index / molar_fractions[index]) * sum([
                molar_fractions[i] * l_vector[i] for i in range(self.num_comp)]))

        ln_gamma_index_res = self.q_2[index] * (1 - np.log(t_vector[index]) - sum([
            theta_2_vector[j] * tau_matrix[index][j] / t_vector[j] for j in range(self.num_comp)]))

        ln_gamma_index = ln_gamma_index_comb + ln_gamma_index_res

        return np.exp(ln_gamma_index)
