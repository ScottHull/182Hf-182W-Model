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
fO2 = -2.25
body_mass = 2.59 * (10**20)
body_radius = 262.7 * 1000
conc_bulk_182hf = 20.692 * 10**-9
conc_bulk_184w = 23.878 * 10**-9
radius_vestian_core = 113.0 * 1000.0
volume_vestian_core = (4.0 / 3.0) * pi * (radius_vestian_core**3.0)
mass_vestian_core = volume_vestian_core * density_metal
max_modeling_time = 20 * 10**6

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



# output_file = 'model_out.csv'
# if output_file in os.listdir(os.getcwd()):
#     os.remove(output_file)
# headers = ['time', 'fraction_of_core_accumulated', 'core_mass', 'mantle_mass', 'core_mass_added', 'cmb_depth',
#            'temperature', 'pressure', 'd', 'bulk_conc_182Hf', 'bulk_conc_182W', '182W_mass_exchange',
#            '184W_mass_exchange', 'silicate_conc_182hf', 'silicate_conc_182w', 'silicate_conc_184w']
# output_file = open(output_file, 'a')
# output_file.write(",".join(headers) + "\n")


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
    m.set_masses()

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



#     output_row = [m.time, m.core_fraction_at_time, m.core_mass, m.mantle_mass, m.metal_mass_added, m.core_mantle_boundary_depth, t, p, d, m.conc_bulk_182hf, (m.conc_core_182w + m.conc_silicate_182w),
#                   mass_exchange_182w, mass_exchange_184w, m.conc_silicate_182hf, m.conc_silicate_182w, m.conc_silicate_184w]
#     output_file.write(",".join([str(i) for i in output_row]) + "\n")
#
# output_file.close()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.axhline(mass_vestian_core, linestyle="--", linewidth=1.0, color='red', label="Modern Vesta Core Mass")
ax1.plot([i / (10**6) for i in time_range], core_mass, linewidth=1.5, color='black', label="Vestian Core Mass")
ax1.grid()
ax1.set_xlabel("Time (Ma)")
ax1.set_ylabel("Core Mass Added at Time (kg)")
ax1.set_title("Exponential Core Growth Model on Vesta")

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot([i / (10**6) for i in time_range], bulk_conc_182hf, linewidth=1.5, color='red', label="Bulk [182Hf]")
ax2.plot([i / (10**6) for i in time_range], mantle_conc_182w, linewidth=1.5, color='blue', label="Mantle [182W]")
ax2.plot([i / (10**6) for i in time_range], mantle_conc_184w, linewidth=1.5, color='green', label="Mantle [184W]")
ax2.grid()
ax2.set_xlabel("Time (Ma)")
ax2.set_ylabel("Concentration (kg)")
ax2.set_title("Concentrations on Vesta")
ax2.legend(loc='upper right')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot([i / (10**6) for i in time_range], moles_bulk_182hf, linewidth=1.5, color='red', label="Bulk 182Hf")
ax3.plot([i / (10**6) for i in time_range], moles_bulk_182w, linewidth=1.5, color='blue', label="Bulk 182W")
ax3.plot([i / (10**6) for i in time_range], [sum(x) for x in zip(moles_bulk_182hf, moles_bulk_182w)], linewidth=1.5,
         color='green', label="Bulk Moles 182Hf + 182W")
ax3.grid()
ax3.set_xlabel("Time (Ma)")
ax3.set_ylabel("Moles")
ax3.set_title("Moles on Vesta")
ax3.legend(loc='lower right')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot([i / (10**6) for i in time_range], core_conc_182w, linewidth=1.5, color='blue', label="Core [182W]")
ax4.plot([i / (10**6) for i in time_range], core_conc_184w, linewidth=1.5, color='green', label="Core [184W]")
ax4.grid()
ax4.set_xlabel("Time (Ma)")
ax4.set_ylabel("Concentration (kg)")
ax4.set_title("Concentrations on Vesta")
ax4.legend(loc='upper right')

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot([i / (10**6) for i in time_range], epsilon_182w_bulk, linewidth=1.5, color='blue', label="Epsilon Bulk")
ax5.plot([i / (10**6) for i in time_range], epsilon_182w_mantle, linewidth=1.5, color='red', label="Epsilon Mantle")
ax5.plot([i / (10**6) for i in time_range], epsilon_182w_core, linewidth=1.5, color='green', label="Epsilon Core")
ax5.grid()
ax5.set_xlabel("Time (Ma)")
ax5.set_ylabel("Epsilon 182W")
ax5.set_title("Epsilon 182W on Vesta")
ax5.legend(loc='upper left')


plt.show()

















# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.axhline(mass_vestian_core, linestyle="--", linewidth=1.0, color='red', label="Modern Vesta Core Mass")
# ax1.plot([i / (10**6) for i in time_range], metal_added_at_time, linewidth=1.5, color='black', label="Vestian Core Mass")
# ax1.grid()
# ax1.set_xlabel("Time (Ma)")
# ax1.set_ylabel("Core Mass Added at Time (kg)")
# ax1.set_title("Exponential Core Growth Model on Vesta")
#
# # fig2 = plt.figure()
# # ax2 = fig2.add_subplot(111)
# # ax2.plot([i / (10**6) for i in time_range], bulk_mass_182w, linewidth=1.5, color='black', label="182W Bulk Mass")
# # ax2.plot([i / (10**6) for i in time_range], core_mass_182w, linewidth=1.5, color='red', label="182W Core Mass")
# # ax2.plot([i / (10**6) for i in time_range], mantle_mass_182w, linewidth=1.5, color='blue', label="182W Mantle Mass")
# # ax2.grid()
# # ax2.set_xlabel("Time (Ma)")
# # ax2.set_ylabel("Bulk Mass 182W (kg)")
# # ax2.set_title("Bulk Mass 182W in Vesta")
#
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111)
# ax3.plot([i / (10**6) for i in time_range], bulk_moles_182w, linewidth=1.5, color='black', label="182W Bulk Moles")
# ax3.plot([i / (10**6) for i in time_range], bulk_moles_182hf, linewidth=1.5, color='red', label="182Hf Bulk Moles")
# ax3.grid()
# ax3.set_xlabel("Time (Ma)")
# ax3.set_ylabel("Bulk Moles")
# ax3.set_title("Bulk Moles in Vesta")
# ax3.legend(loc='lower right')
#
# fig4 = plt.figure()
# ax4 = fig4.add_subplot(111)
# ax4.plot([i / (10**6) for i in time_range], [i / 1000 for i in cmb_depth], linewidth=1.5, color='black', label="CMB Depth")
# ax4.axhline(radius_vestian_core / 1000, linestyle="--", color='red', label='Modern Core Radius')
# ax4.grid()
# ax4.set_xlabel("Time (Ma)")
# ax4.set_ylabel("CMB Depth (km)")
# ax4.set_title("CMB Depth in Vesta")
# ax4.legend(loc='lower right')
#
# # fig5 = plt.figure()
# # ax5 = fig5.add_subplot(111)
# # ax5.plot([i / (10**6) for i in time_range], [i * 10**9 for i in core_conc_184w], linewidth=1.5, color='blue', label="184W")
# # ax5.plot([i / (10**6) for i in time_range], [i * 10**9 for i in core_conc_182w], linewidth=1.5, color='red', label="182W")
# # ax5.grid()
# # ax5.set_xlabel("Time (Ma)")
# # ax5.set_ylabel("W Isotope Concentration (ppb)")
# # ax5.set_title("W Isotope Concentration in Vesta")
# # ax5.legend(loc='lower right')
# #
# # fig6 = plt.figure()
# # ax6 = fig6.add_subplot(111)
# # ax6.axhline(mass_vestian_core, linestyle="--", linewidth=1.0, color='red', label="Modern Vesta Core Mass")
# # ax6.plot([i / (10**6) for i in time_range], core_mass, linewidth=1.5, color='black', label="Vestian Core Mass")
# # ax6.grid()
# # ax6.set_xlabel("Time (Ma)")
# # ax6.set_ylabel("Core Mass (kg)")
# # ax6.set_title("Exponential Core Growth Model on Vesta")
#
# plt.show()
