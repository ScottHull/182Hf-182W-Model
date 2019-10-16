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
conc_bulk_182hf = 20.692 * 10**-9
conc_bulk_184w = 23.878 * 10**-9
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

time_range = []


metal_added_at_time = []
core_mass = []
mantle_mass = []

bulk_mass_182w = []
core_mass_182w = []
mantle_mass_182w = []
core_mass_184w = []
mantle_mass_184w = []

bulk_moles_182hf = []
bulk_moles_182w = []

cmb_depth = []

bulk_conc_182hf = []
core_conc_184w = []
core_conc_182w = []
mantle_conc_182hf = []
mantle_conc_182w = []
mantle_conc_184w = []

moles_bulk_182hf = []
moles_bulk_182w = []

epsilon_182w_bulk = []
epsilon_182w_mantle = []
epsilon_182w_core = []

fO2_and_silicate_epsilon_infinite_time = []



# output_file = 'model_out.csv'
# if output_file in os.listdir(os.getcwd()):
#     os.remove(output_file)
# headers = ['time', 'fraction_of_core_accumulated', 'core_mass', 'mantle_mass', 'core_mass_added', 'cmb_depth',
#            'temperature', 'pressure', 'd', 'bulk_conc_182Hf', 'bulk_conc_182W', '182W_mass_exchange',
#            '184W_mass_exchange', 'silicate_conc_182hf', 'silicate_conc_182w', 'silicate_conc_184w']
# output_file = open(output_file, 'a')
# output_file.write(",".join(headers) + "\n")


for i in np.arange(-2.25, 0, 0.01):
    # build the model
    m = Model(
        body_core_mass=mass_vestian_core,
        body_mass=body_mass,
        body_radius=body_radius,
        timestep=timestep,
        max_time_core_formation=max_core_formation_time,
        fO2=i,
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
        metal_added_at_time.append(m.metal_mass_added)
        core_mass.append(m.core_mass)
        mantle_mass.append(m.mantle_mass)

        # decay 182Hf to 182W
        m.decay_182hf()
        bulk_mass_182w.append(m.mass_bulk_182w)
        bulk_moles_182hf.append(m.moles_bulk_182hf)
        bulk_moles_182w.append(m.moles_bulk_182w)

        # calculate partitioning
        t = m.adiabat()
        p = m.hydrostat()
        d = m.partitioning(temperature=t, pressure=p)
        cmb_depth.append(m.core_mantle_boundary_depth)

        if m.time <= m.max_time_core_formation:
            # calculate equilibrium exchange
            mass_exchange_182w = m.exchange_mass(D=d, mass_w_silicate=m.mass_silicate_182w)
            mass_exchange_184w = m.exchange_mass(D=d, mass_w_silicate=m.mass_silicate_184w)
        else:
            mass_exchange_182w = 0
            mass_exchange_184w = 0

        # equilibrate
        m.set_concentrations(exchange_mass_182w=mass_exchange_182w, exchange_mass_184w=mass_exchange_184w)
        # m.set_masses()

        e_182w_bulk = m.calculate_epsilon182w(sample_ratio=(m.mass_bulk_182w / m.mass_bulk_184w))
        e_182w_mantle = m.calculate_epsilon182w(sample_ratio=(m.mass_silicate_182w / m.mass_silicate_184w))
        e_182w_core = m.calculate_epsilon182w(sample_ratio=(m.mass_core_182w / m.mass_core_184w))

        bulk_conc_182hf.append(m.conc_bulk_182hf)

        core_conc_182w.append(m.conc_core_182w)
        core_conc_184w.append(m.conc_core_184w)

        mantle_conc_182hf.append(m.conc_silicate_182hf)
        mantle_conc_182w.append(m.conc_silicate_182w)
        mantle_conc_184w.append(m.conc_silicate_184w)

        core_mass_182w.append(m.mass_core_182w)
        mantle_mass_182w.append(m.mass_silicate_182w)
        core_mass_184w.append(m.mass_core_184w)
        mantle_mass_184w.append(m.mass_silicate_184w)

        moles_bulk_182hf.append(m.moles_bulk_182hf)
        moles_bulk_182w.append(m.moles_bulk_182w)

        epsilon_182w_bulk.append(e_182w_bulk)
        epsilon_182w_mantle.append(e_182w_mantle)
        epsilon_182w_core.append(e_182w_core)

    fO2_and_silicate_epsilon_infinite_time.append((i, m.calculate_epsilon182w(sample_ratio=(m.mass_silicate_182w / m.mass_silicate_184w))))



#     output_row = [m.time, m.core_fraction_at_time, m.core_mass, m.mantle_mass, m.metal_mass_added, m.core_mantle_boundary_depth, t, p, d, m.conc_bulk_182hf, (m.conc_core_182w + m.conc_silicate_182w),
#                   mass_exchange_182w, mass_exchange_184w, m.conc_silicate_182hf, m.conc_silicate_182w, m.conc_silicate_184w]
#     output_file.write(",".join([str(i) for i in output_row]) + "\n")
#
# output_file.close()
