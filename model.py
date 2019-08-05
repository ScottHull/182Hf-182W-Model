from math import sqrt, log, exp, pi
from copy import copy

class Model:

    def __init__(
            self, body_core_mass, body_mass, body_radius, timestep, max_time_core_formation, fO2, core_fraction_per_timestep,
            conc_bulk_182hf, alpha, beta, chi, delta, epsilon, surface_temperature, surface_pressure,
            silicate_thermal_expansivity, gravitational_acceleration, silicate_heat_capacity
    ):

        self.__half_life_182hf = 8.9 * 10**6  # the half life of 182Hf, do not set
        self.__decay_const_182hf = log(0.5) / self.__half_life_182hf  # the decay constant of 182Hf, do not set
        self.body_core_mass = body_core_mass  # the core mass of the modern body, set
        self.body_mass = body_mass  # the bulk mass of the modern body, set
        self.body_radius = body_radius
        self.timestep = timestep  # the modeling timestep, set
        self.time = 0  # the model time, do not set
        self.max_time_core_formation = max_time_core_formation  # the maximum timing of core formation, set
        self.fO2 = fO2

        self.__molar_mass_182w = 183.84
        self.__molar_mass_182hf = 178.49

        self.core_mass = 0  # the modeled core mass at time, do not set
        self.silicate_mass = 0
        self.metal_mass_added = 0
        self.core_fraction_per_timestep = core_fraction_per_timestep  # the fraction of metal to be added to the core, set
        self.__silicate_mass = 0  # the mass of the silicate magma ocean at time, do not set

        self.conc_bulk_182hf = conc_bulk_182hf
        self.mass_bulk_182hf = self.conc_bulk_182hf * body_mass
        self.moles_bulk_182hf = self.mass_bulk_182hf / self.__molar_mass_182hf
        self.conc_silicate_182hf = 0
        self.__original_moles_bulk_182hf = copy(self.moles_bulk_182hf)
        self.__previous_timestep_moles_bulk_182hf = copy(self.moles_bulk_182hf)

        self.__conc_silicate_182w = 0
        self.moles_bulk_182w = 0
        self.mass_bulk_182w = 0
        self.__conc_metal_182w = 0
        self.conc_core_182w = 0

        self.__conc_silicate_184w = 0
        self.__conc_metal_184w = 0
        self.__conc_core_184w = 0

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

        self.core_mantle_boundary_depth = 0

        self.volume_core = 0
        self.density_metal = 7800
        self.density_silicate = 3750

    def adiabat(self):

        T = self.surface_temperature * (1 + (self.silicate_thermal_expansivity * self.gravitational_acceleration
                                             * self.core_mantle_boundary_depth) / self.silicate_heat_capacity)

        return T

    def hydrostat(self):

        P = self.surface_pressure + (self.density_silicate * self.gravitational_acceleration *
                                     self.core_mantle_boundary_depth) * 10**-9

        return P

    def partitioning(self, temperature, pressure):

        logD = self.alpha + (self.beta * self.fO2) + (self.chi * self.nbo_t) + (self.delta * (1 / temperature)) + \
               (self.epislon * (pressure / temperature))

        D = 10**logD

        return D


    def partition_w(self, D):

        pass


    def fractionate_metal_exponential(self):

        core_growth_step_function = 1
        if self.time <= self.max_time_core_formation:
            core_growth_step_function = 1
        else:
            core_growth_step_function = 1 / (1 - exp(self.__core_growth_constant * self.time))

        self.__core_growth_constant = log(0.01) / self.max_time_core_formation
        previous_core_fraction_per_timestep = copy(self.core_fraction_per_timestep)
        self.core_fraction_per_timestep = core_growth_step_function * (1 - exp(self.__core_growth_constant * self.time))
        self.core_mass = self.body_core_mass * self.core_fraction_per_timestep
        self.mantle_mass = self.body_mass - self.core_mass
        self.metal_mass_added = (self.core_fraction_per_timestep - previous_core_fraction_per_timestep) * \
                                  self.body_core_mass

        self.volume_core = (self.core_mass / self.density_metal)
        self.core_mantle_boundary_depth = (self.volume_core / ((4 / 3) * pi))**(1 / 3)


    def decay_182hf(self):

        self.__previous_timestep_moles_bulk_182hf = self.moles_bulk_182hf
        self.moles_bulk_182hf = self.moles_bulk_182hf * exp(self.__decay_const_182hf * self.timestep)
        self.moles_bulk_182w += self.__previous_timestep_moles_bulk_182hf - self.moles_bulk_182hf
        self.mass_bulk_182w = self.moles_bulk_182w * self.__molar_mass_182w
        self.mass_bulk_182hf = (self.__original_moles_bulk_182hf * self.__molar_mass_182hf) - self.mass_bulk_182w


    def set_concentrations(self):

        self.conc_silicate_182hf = self.mass_bulk_182hf / self.silicate_mass
        self.conc_bulk_182hf = self.conc_silicate_182hf
