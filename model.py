from math import sqrt, log, exp, pi
from copy import copy

class Model:

    def __init__(
            self, body_core_mass, body_mass, body_radius, timestep, max_time_core_formation, fO2, core_fraction_per_timestep,
            conc_bulk_182hf, conc_bulk_184w, alpha, beta, chi, delta, epsilon, surface_temperature, surface_pressure,
            silicate_thermal_expansivity, gravitational_acceleration, silicate_heat_capacity
    ):

        # general modeling parameters
        self.timestep = timestep  # the modeling timestep, set
        self.time = 0.0  # the model time, do not set
        self.max_time_core_formation = max_time_core_formation  # the maximum timing of core formation, set
        self.__molar_mass_w = 183.84 / 1000  # kg/mol
        self.__molar_mass_hf = 178.49 / 1000  # kg/mol
        self.__half_life_182hf = 8.9 * 10 ** 6  # the half life of 182Hf, do not set
        self.__decay_const_182hf = log(0.5) / self.__half_life_182hf  # the decay constant of 182Hf, do not set
        self.alpha = alpha
        self.beta = beta
        self.chi = chi
        self.delta = delta
        self.epislon = epsilon

        self.nbo_t = 2.6

        self.surface_temperature = surface_temperature
        self.surface_pressure = surface_pressure

        self.gravitational_acceleration = gravitational_acceleration
        self.silicate_thermal_expansivity = silicate_thermal_expansivity
        self.silicate_heat_capacity = silicate_heat_capacity

        self.core_mantle_boundary_depth = 0.0

        # environment parameters
        self.fO2 = fO2
        self.density_metal = 7800
        self.density_silicate = 3750

        # planet/planetesimal evolution parameters
        self.body_core_mass = body_core_mass  # the core mass of the modern body, set
        self.body_mass = body_mass  # the bulk mass of the modern body, set
        self.body_radius = body_radius
        self.volume_core = 0
        self.volume_mantle = (4 / 3) * pi * (self.body_radius) ** 3
        self.core_mass = 0.0  # the modeled core mass at time, do not set
        self.mantle_mass = self.volume_mantle * self.density_silicate
        self.metal_mass_added = 0.0
        self.core_fraction_at_time = 0.0
        self.core_fraction_per_timestep = core_fraction_per_timestep  # the fraction of metal to be added to the core
        self.core_fraction_at_time = 0

        # 182hf parameters
        self.conc_bulk_182hf = conc_bulk_182hf
        self.mass_bulk_182hf = self.conc_bulk_182hf * body_mass
        self.mass_silicate_182hf = self.mass_bulk_182hf
        self.moles_bulk_182hf = self.mass_bulk_182hf / self.__molar_mass_hf
        self.conc_silicate_182hf = 0
        self.__original_moles_bulk_182hf = copy(self.moles_bulk_182hf)
        self.__previous_timestep_moles_bulk_182hf = copy(self.moles_bulk_182hf)

        # 182w parameters
        self.conc_silicate_182w = 0.0
        self.moles_bulk_182w = 0.0
        self.mass_silicate_added_182w = 0.0
        self.mass_bulk_182w = 0.0
        self.__conc_metal_182w = 0.0
        self.conc_core_182w = 0.0
        self.conc_core_bound_metal_182w = 0.0
        self.mass_core_bound_metal_182w = 0.0
        self.mass_core_182w = 0.0
        self.mass_silicate_182w = 0.0

        # 184w parameters
        self.conc_bulk_184w = conc_bulk_184w
        self.conc_silicate_184w = (self.conc_bulk_184w * self.body_mass) / self.mantle_mass
        self.conc_core_184w = 0.0
        self.conc_core_bound_metal_184w = 0.0
        self.mass_core_bound_metal_184w = 0.0
        self.mass_bulk_184w = self.conc_bulk_184w * self.body_mass
        self.mass_silicate_184w = (self.conc_silicate_184w * self.mantle_mass)
        self.mass_core_184w = 0.0

    def adiabat(self):

        T = self.surface_temperature * (1 + (self.silicate_thermal_expansivity * self.gravitational_acceleration
                                             * self.core_mantle_boundary_depth) / self.silicate_heat_capacity)

        return T

    def hydrostat(self):

        P = (self.surface_pressure + (self.density_silicate * self.gravitational_acceleration *
                                     self.core_mantle_boundary_depth)) * 10**-9

        return P

    def partitioning(self, temperature, pressure):

        logD = self.alpha + (self.beta * self.fO2) + (self.chi * self.nbo_t) + (self.delta * (1 / temperature)) + \
               (self.epislon * (pressure / temperature))

        D = 10**logD

        return D

    def fractionate_metal_exponential(self):

        core_growth_step_function = 1.0
        if self.time <= self.max_time_core_formation:
            core_growth_step_function = 1.0
        else:
            core_growth_step_function = 1.0 / (1.0 - exp(self.__core_growth_constant * self.time))

        if self.time == 0.0:
            self.core_fraction_per_timestep = 0.0

        self.__core_growth_constant = log(0.01) / self.max_time_core_formation
        previous_core_fraction_per_timestep = copy(self.core_fraction_at_time)
        previous_core_mass = copy(self.core_mass)
        self.core_fraction_at_time = core_growth_step_function * (1.0 - exp(self.__core_growth_constant * self.time))
        self.core_fraction_per_timestep = self.core_fraction_at_time - previous_core_fraction_per_timestep
        self.core_mass = self.body_core_mass * self.core_fraction_at_time
        self.mantle_mass = self.body_mass - self.core_mass
        self.metal_mass_added = self.core_mass - previous_core_mass
        self.volume_core = (self.core_mass / self.density_metal)
        self.volume_mantle = self.volume_mantle - self.volume_core
        self.core_mantle_boundary_depth = self.body_radius - ((self.volume_core / ((4 / 3) * pi))**(1 / 3))


    def decay_182hf(self):

        self.__previous_timestep_moles_bulk_182hf = copy(self.moles_bulk_182hf)
        self.__previous_mass_bulk_182w = copy(self.mass_bulk_182w)
        self.__previous_mass_silicate_182w = copy(self.mass_silicate_182w)
        self.moles_bulk_182hf = self.__previous_timestep_moles_bulk_182hf * \
                                exp(self.__decay_const_182hf * self.timestep)
        self.moles_bulk_182w += self.__previous_timestep_moles_bulk_182hf - self.moles_bulk_182hf
        self.mass_bulk_182w += (self.__previous_timestep_moles_bulk_182hf - self.moles_bulk_182hf) * self.__molar_mass_w
        self.mass_silicate_added_182w = self.mass_bulk_182w - self.__previous_mass_bulk_182w
        self.mass_silicate_182w = self.__previous_mass_silicate_182w + self.mass_silicate_added_182w
        self.mass_bulk_182hf = self.moles_bulk_182hf * self.__molar_mass_hf
        self.mass_silicate_182hf = self.mass_bulk_182hf

        self.conc_silicate_182hf = self.mass_silicate_182hf / self.mantle_mass
        self.conc_bulk_182hf = self.mass_silicate_182hf / self.body_mass
        self.conc_silicate_182w = self.mass_silicate_182w / self.mantle_mass

        self.conc_silicate_184w = self.mass_silicate_184w / self.mantle_mass


    def exchange_mass(self, D, mass_w_silicate):

        """
        The exchange mass of the either 182W or 184W prior to equilibration.  This function is generic in order to
        support both 182W and 184W.
        :param D:
        :param conc_w_silicate:
        :return:
        """

        # conc_w_silicate should be the concentration of the isotope in the silicate at time t - \Delta t

        numerator = D * (mass_w_silicate / self.mantle_mass)
        denominator = (1.0 / self.metal_mass_added) + (D / self.mantle_mass)
        m = numerator / denominator

        return m

    def set_concentrations(self, exchange_mass_182w, exchange_mass_184w):

        if self.time <= self.max_time_core_formation:
            # 182w
            # self.conc_silicate_182w = self.conc_silicate_182w - (exchange_mass_182w / self.mantle_mass)
            self.mass_silicate_182w -= exchange_mass_182w
            self.conc_silicate_182w = self.mass_silicate_182w / self.mantle_mass
            self.conc_core_bound_metal_182w = exchange_mass_182w / self.metal_mass_added
            self.mass_core_182w += exchange_mass_182w
            self.conc_core_182w = self.mass_core_182w / self.core_mass

            # 184w
            # self.conc_silicate_184w = self.conc_silicate_184w - (exchange_mass_184w / self.mantle_mass)
            self.mass_silicate_184w -= exchange_mass_184w
            self.conc_silicate_184w = self.mass_silicate_184w / self.mantle_mass
            self.conc_core_bound_metal_184w = exchange_mass_184w / self.metal_mass_added
            self.mass_core_184w += exchange_mass_184w
            self.conc_core_184w = self.mass_core_184w / self.core_mass

    # def set_masses(self):
    #
    #     # 182w
    #     self.mass_core_182w += (self.mass_core_bound_metal_182w * self.conc_core_bound_metal_182w)
    #     self.mass_silicate_182w = self.mantle_mass * self.conc_silicate_182w
    #
    #     # 184w
    #     self.mass_core_184w += (self.mass_core_bound_metal_184w * self.conc_core_bound_metal_184w)
    #     self.mass_silicate_184w = self.mantle_mass * self.conc_silicate_184w

    def calculate_epsilon182w(self, sample_ratio, standard_ratio=0.86468):

        epsilon = ((sample_ratio / standard_ratio) - 1.0) * (10.0**4.0)

        return epsilon








class Model_SSI:

    def __init__(
            self, body_core_mass, body_mass, body_radius, timestep, max_time_core_formation, fO2, core_fraction_per_timestep,
            conc_bulk_182hf, conc_bulk_184w, alpha, beta, chi, delta, epsilon, surface_temperature, surface_pressure,
            silicate_thermal_expansivity, gravitational_acceleration, silicate_heat_capacity, w182_184_ssi
    ):

        # general modeling parameters
        self.timestep = timestep  # the modeling timestep, set
        self.time = 0.0  # the model time, do not set
        self.max_time_core_formation = max_time_core_formation  # the maximum timing of core formation, set
        self.__molar_mass_w = 183.84 / 1000  # kg/mol
        self.__molar_mass_hf = 178.49 / 1000  # kg/mol
        self.__half_life_182hf = 8.9 * 10 ** 6  # the half life of 182Hf, do not set
        self.__decay_const_182hf = log(0.5) / self.__half_life_182hf  # the decay constant of 182Hf, do not set
        self.alpha = alpha
        self.beta = beta
        self.chi = chi
        self.delta = delta
        self.epislon = epsilon

        self.nbo_t = 2.6

        self.surface_temperature = surface_temperature
        self.surface_pressure = surface_pressure

        self.gravitational_acceleration = gravitational_acceleration
        self.silicate_thermal_expansivity = silicate_thermal_expansivity
        self.silicate_heat_capacity = silicate_heat_capacity

        self.core_mantle_boundary_depth = 0.0

        # environment parameters
        self.fO2 = fO2
        self.density_metal = 7800
        self.density_silicate = 3750

        # planet/planetesimal evolution parameters
        self.body_core_mass = body_core_mass  # the core mass of the modern body, set
        self.body_mass = body_mass  # the bulk mass of the modern body, set
        self.body_radius = body_radius
        self.volume_core = 0
        self.volume_mantle = (4 / 3) * pi * (self.body_radius) ** 3
        self.core_mass = 0.0  # the modeled core mass at time, do not set
        self.mantle_mass = self.volume_mantle * self.density_silicate
        self.metal_mass_added = 0.0
        self.core_fraction_at_time = 0.0
        self.core_fraction_per_timestep = core_fraction_per_timestep  # the fraction of metal to be added to the core
        self.core_fraction_at_time = 0

        # 182hf parameters
        self.conc_bulk_182hf = conc_bulk_182hf
        self.mass_bulk_182hf = self.conc_bulk_182hf * body_mass
        self.mass_silicate_182hf = self.mass_bulk_182hf
        self.moles_bulk_182hf = self.mass_bulk_182hf / self.__molar_mass_hf
        self.conc_silicate_182hf = 0
        self.__original_moles_bulk_182hf = copy(self.moles_bulk_182hf)
        self.__previous_timestep_moles_bulk_182hf = copy(self.moles_bulk_182hf)

        # 184w parameters
        self.conc_bulk_184w = conc_bulk_184w
        self.conc_silicate_184w = (self.conc_bulk_184w * self.body_mass) / self.mantle_mass
        self.conc_core_184w = 0.0
        self.conc_core_bound_metal_184w = 0.0
        self.mass_core_bound_metal_184w = 0.0
        self.mass_bulk_184w = self.conc_bulk_184w * self.body_mass
        self.mass_silicate_184w = (self.conc_silicate_184w * self.mantle_mass)
        self.mass_core_184w = 0.0

        # 182w parameters
        self.moles_bulk_182w = self.conc_bulk_184w * w182_184_ssi
        self.mass_silicate_added_182w = 0.0
        self.mass_bulk_182w = self.moles_bulk_182w * self.__molar_mass_w
        self.__conc_metal_182w = 0.0
        self.conc_core_182w = 0.0
        self.conc_core_bound_metal_182w = 0.0
        self.mass_core_bound_metal_182w = 0.0
        self.mass_core_182w = 0.0
        self.mass_silicate_182w = self.mass_bulk_182w
        self.conc_silicate_182w = self.mass_silicate_182w / self.body_mass

    def adiabat(self):

        T = self.surface_temperature * (1 + (self.silicate_thermal_expansivity * self.gravitational_acceleration
                                             * self.core_mantle_boundary_depth) / self.silicate_heat_capacity)

        return T

    def hydrostat(self):

        P = (self.surface_pressure + (self.density_silicate * self.gravitational_acceleration *
                                     self.core_mantle_boundary_depth)) * 10**-9

        return P

    def partitioning(self, temperature, pressure):

        logD = self.alpha + (self.beta * self.fO2) + (self.chi * self.nbo_t) + (self.delta * (1 / temperature)) + \
               (self.epislon * (pressure / temperature))

        D = 10**logD

        return D

    def fractionate_metal_exponential(self):

        core_growth_step_function = 1.0
        if self.time <= self.max_time_core_formation:
            core_growth_step_function = 1.0
        else:
            core_growth_step_function = 1.0 / (1.0 - exp(self.__core_growth_constant * self.time))

        if self.time == 0.0:
            self.core_fraction_per_timestep = 0.0

        self.__core_growth_constant = log(0.01) / self.max_time_core_formation
        previous_core_fraction_per_timestep = copy(self.core_fraction_at_time)
        previous_core_mass = copy(self.core_mass)
        self.core_fraction_at_time = core_growth_step_function * (1.0 - exp(self.__core_growth_constant * self.time))
        self.core_fraction_per_timestep = self.core_fraction_at_time - previous_core_fraction_per_timestep
        self.core_mass = self.body_core_mass * self.core_fraction_at_time
        self.mantle_mass = self.body_mass - self.core_mass
        self.metal_mass_added = self.core_mass - previous_core_mass
        self.volume_core = (self.core_mass / self.density_metal)
        self.volume_mantle = self.volume_mantle - self.volume_core
        self.core_mantle_boundary_depth = self.body_radius - ((self.volume_core / ((4 / 3) * pi))**(1 / 3))


    def decay_182hf(self):

        self.__previous_timestep_moles_bulk_182hf = copy(self.moles_bulk_182hf)
        self.__previous_mass_bulk_182w = copy(self.mass_bulk_182w)
        self.__previous_mass_silicate_182w = copy(self.mass_silicate_182w)
        self.moles_bulk_182hf = self.__previous_timestep_moles_bulk_182hf * \
                                exp(self.__decay_const_182hf * self.timestep)
        self.moles_bulk_182w += self.__previous_timestep_moles_bulk_182hf - self.moles_bulk_182hf
        self.mass_bulk_182w += (self.__previous_timestep_moles_bulk_182hf - self.moles_bulk_182hf) * self.__molar_mass_w
        self.mass_silicate_added_182w = self.mass_bulk_182w - self.__previous_mass_bulk_182w
        self.mass_silicate_182w = self.__previous_mass_silicate_182w + self.mass_silicate_added_182w
        self.mass_bulk_182hf = self.moles_bulk_182hf * self.__molar_mass_hf
        self.mass_silicate_182hf = self.mass_bulk_182hf

        self.conc_silicate_182hf = self.mass_silicate_182hf / self.mantle_mass
        self.conc_bulk_182hf = self.mass_silicate_182hf / self.body_mass
        self.conc_silicate_182w = self.mass_silicate_182w / self.mantle_mass

        self.conc_silicate_184w = self.mass_silicate_184w / self.mantle_mass


    def exchange_mass(self, D, mass_w_silicate):

        """
        The exchange mass of the either 182W or 184W prior to equilibration.  This function is generic in order to
        support both 182W and 184W.
        :param D:
        :param conc_w_silicate:
        :return:
        """

        # conc_w_silicate should be the concentration of the isotope in the silicate at time t - \Delta t

        numerator = D * (mass_w_silicate / self.mantle_mass)
        denominator = (1.0 / self.metal_mass_added) + (D / self.mantle_mass)
        m = numerator / denominator

        return m

    def set_concentrations(self, exchange_mass_182w, exchange_mass_184w):

        if self.time <= self.max_time_core_formation:
            # 182w
            # self.conc_silicate_182w = self.conc_silicate_182w - (exchange_mass_182w / self.mantle_mass)
            self.mass_silicate_182w -= exchange_mass_182w
            self.conc_silicate_182w = self.mass_silicate_182w / self.mantle_mass
            self.conc_core_bound_metal_182w = exchange_mass_182w / self.metal_mass_added
            self.mass_core_182w += exchange_mass_182w
            self.conc_core_182w = self.mass_core_182w / self.core_mass

            # 184w
            # self.conc_silicate_184w = self.conc_silicate_184w - (exchange_mass_184w / self.mantle_mass)
            self.mass_silicate_184w -= exchange_mass_184w
            self.conc_silicate_184w = self.mass_silicate_184w / self.mantle_mass
            self.conc_core_bound_metal_184w = exchange_mass_184w / self.metal_mass_added
            self.mass_core_184w += exchange_mass_184w
            self.conc_core_184w = self.mass_core_184w / self.core_mass

    # def set_masses(self):
    #
    #     # 182w
    #     self.mass_core_182w += (self.mass_core_bound_metal_182w * self.conc_core_bound_metal_182w)
    #     self.mass_silicate_182w = self.mantle_mass * self.conc_silicate_182w
    #
    #     # 184w
    #     self.mass_core_184w += (self.mass_core_bound_metal_184w * self.conc_core_bound_metal_184w)
    #     self.mass_silicate_184w = self.mantle_mass * self.conc_silicate_184w

    def calculate_epsilon182w(self, sample_ratio, standard_ratio=0.86468):

        epsilon = ((sample_ratio / standard_ratio) - 1.0) * (10.0**4.0)

        return epsilon

