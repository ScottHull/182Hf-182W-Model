from model import Model
import os
from math import pi
import numpy as np
import matplotlib.pyplot as plt

density_silicate_melt = 3750
density_metal = 7800
max_core_formation_time = 5.0 * 10**6
timestep = 5.0 * 10**5
core_fraction_per_timestep = max_core_formation_time / timestep
body_mass = 2.59 * (10**20)
body_radius = 262.7 * 1000

radius_vestian_core = 113.0 * 1000.0
volume_vestian_core = (4.0 / 3.0) * pi * (radius_vestian_core**3.0)
mass_vestian_core = volume_vestian_core * density_metal
max_modeling_time = 100 * 10**6

alpha = 1.11
beta = -1.18
chi = -0.85
delta = 1680
epsilon = 487

surface_temp = 2000.0
surface_pressure = 0
therm_expansivity = 6.0 * (10**(-5))
heat_capacity = 10**3
gravity = 0.25

fO2_and_starting_concs = {
    -1.5:
        {
            '182hf': 20.20881742605 * 10**-9,
            '184w': 1785.959333637 * 10**-9,
        },
    -1.0:
        {
            '182hf': 18.3219620165 * 10**-9,
            '184w': 92.4155533335 * 10**-9,
        },
    -0.5:
         {
             '182hf': 17.0951688211 * 10**-9,
             '184w': 30.605062812 * 10**-9,
         },
    0.0:
        {
            '182hf': 16.625296305 * 10**-9,
            '184w': 22.0148036578 * 10**-9,
        },
    0.5:
        {
            '182hf': 16.4886194373 * 10**-9,
            '184w': 20.1479869741 * 10**-9,
        },
    1.0:
        {
            '182hf': 16.4522891941 * 10**-9,
            '184w': 19.6887243759 * 10**-9,
        },
}

for fO2 in fO2_and_starting_concs.keys():
    time_range = []


    conc_bulk_182hf = fO2_and_starting_concs[fO2]['182hf']
    conc_bulk_184w = fO2_and_starting_concs[fO2]['184w']


    # build the model
    m = Model(
        body_core_mass=mass_vestian_core,
        body_mass=body_mass,
        body_radius=body_radius,
        timestep=timestep,
        max_time_core_formation=max_core_formation_time,
        fO2=fO2,
        core_fraction_per_timestep=core_fraction_per_timestep,
        conc_bulk_182hf=conc_bulk_182hf,
        conc_bulk_184w=conc_bulk_184w,
        alpha=alpha,
        beta=beta,
        chi=chi,
        delta=delta,
        epsilon=epsilon,
        surface_temperature=surface_temp,
        surface_pressure=surface_pressure,
        silicate_thermal_expansivity=therm_expansivity,
        gravitational_acceleration=gravity,
        silicate_heat_capacity=heat_capacity,
    )

    while m.time < max_modeling_time:

        # advance the modeled time
        m.time += m.timestep
        time_range.append(m.time)
        print("At time {} Myr".format(m.time / float(10**6)))

        # fractionate some metal
        m.fractionate_metal_exponential()

        # decay 182Hf to 182W
        m.decay_182hf()

        # calculate partitioning
        t = m.adiabat()
        p = m.hydrostat()
        d = m.partitioning(temperature=t, pressure=p)

        if m.time <= m.max_time_core_formation:
            # calculate equilibrium exchange
            mass_exchange_182w = m.exchange_mass(D=d, mass_w_silicate=m.mass_silicate_182w)
            mass_exchange_184w = m.exchange_mass(D=d, mass_w_silicate=m.mass_silicate_184w)
        else:
            mass_exchange_182w = 0
            mass_exchange_184w = 0

        # equilibrate
        m.set_concentrations(exchange_mass_182w=mass_exchange_182w, exchange_mass_184w=mass_exchange_184w)

        e_182w_bulk = m.calculate_epsilon182w(sample_ratio=(m.mass_bulk_182w / m.mass_bulk_184w))
        e_182w_mantle = m.calculate_epsilon182w(sample_ratio=(m.mass_silicate_182w / m.mass_silicate_184w))
        e_182w_core = m.calculate_epsilon182w(sample_ratio=(m.mass_core_182w / m.mass_core_184w))


