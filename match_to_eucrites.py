from model import Model
import os
from math import pi
import numpy as np
import pandas as pd
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

fO2_and_bulk_epsilon_infinite_time = []
fO2_and_silicate_epsilon_infinite_time = []
fO2_and_core_epsilon_infinite_time = []

df_dict = {}



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

    fO2_and_bulk_epsilon_infinite_time.append(
        m.calculate_epsilon182w(sample_ratio=(m.mass_bulk_182w / m.mass_bulk_184w)))
    fO2_and_silicate_epsilon_infinite_time.append(
        m.calculate_epsilon182w(sample_ratio=(m.mass_silicate_182w / m.mass_silicate_184w)))
    fO2_and_core_epsilon_infinite_time.append(
        m.calculate_epsilon182w(sample_ratio=(m.mass_core_182w / m.mass_core_184w)))

    df_dict.update({
        fO2:
            {
                'bulk': m.calculate_epsilon182w(sample_ratio=(m.mass_bulk_182w / m.mass_bulk_184w)),
                'silicate': m.calculate_epsilon182w(sample_ratio=(m.mass_silicate_182w / m.mass_silicate_184w)),
                'core': m.calculate_epsilon182w(sample_ratio=(m.mass_core_182w / m.mass_core_184w)),
            }
    })


df_dict_to_csv = {}
for key in df_dict.keys():
    s_bulk = "{}-{}".format(key, 'bulk')
    s_silicate = "{}-{}".format(key, 'silicate')
    s_core = "{}-{}".format(key, 'core')
    df_dict_to_csv.update({
        s_bulk: [df_dict[key]['bulk']],
        s_silicate: [df_dict[key]['silicate']],
        s_core: [df_dict[key]['core']]
    })
pd.DataFrame(df_dict_to_csv).to_csv('fO2_matched_to_eucrites.csv')





fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(sorted(fO2_and_starting_concs.keys()),
         [fO2_and_starting_concs[key]['182hf'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())], linewidth=2.0, color='blue', label='182Hf')
ax1.plot(sorted(fO2_and_starting_concs.keys()),
         [fO2_and_starting_concs[key]['184w'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())], linewidth=2.0, color='red', label='184W')
ax1.grid()
ax1.set_xlabel("fO2")
ax1.set_ylabel("Concentration (ppb)")
ax1.set_title("Initial Concentrations in Vesta")
ax1.legend(loc='upper right')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(sorted(fO2_and_starting_concs.keys())[2:],
         [fO2_and_starting_concs[key]['182hf'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())][2:], linewidth=2.0, color='blue', label='182Hf')
ax2.plot(sorted(fO2_and_starting_concs.keys())[2:],
         [fO2_and_starting_concs[key]['184w'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())][2:], linewidth=2.0, color='red', label='184W')
ax2.grid()
ax2.set_xlabel("fO2")
ax2.set_ylabel("Concentration (ppb)")
ax2.set_title("Initial Bulk Concentrations in Vesta")
ax2.legend(loc='upper right')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
z = zip([fO2_and_starting_concs[key]['182hf'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())][2:],
        [fO2_and_starting_concs[key]['184w'] * 10**9 for key in sorted(fO2_and_starting_concs.keys())][2:])
ax3.plot(sorted(fO2_and_starting_concs.keys())[2:],
         [i[0] / i[1] for i in z], linewidth=2.0, color='black', label='182Hf/184W')
ax3.axhline(0.866555125, linewidth=2.0, linestyle='--', color='red', label='182W/184W (Vesta at infinite time)')
ax3.grid()
ax3.set_xlabel("fO2")
ax3.set_ylabel("182Hf/184W")
ax3.set_title("Initial Bulk 182Hf/184W in Vesta")
ax3.legend(loc='lower right')



# plt.show()




