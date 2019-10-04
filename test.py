from model import Model
from math import pi
import numpy as np
import matplotlib.pyplot as plt

density_silicate_melt = 3750
density_metal = 7800
max_core_formation_time = 5 * 10**6
timestep = 5 * 10**6
core_fraction_per_timestep = max_core_formation_time / timestep
fO2 = -2.25
body_mass = 2.59 * 10**20
body_radius = 262.7 * 1000
conc_bulk_182hf = 20.692 * 10**-9
conc_bulk_184w = 23.878 * 10**-9
radius_vestian_core = 113 * 1000
volume_vestian_core = (4 / 3) * pi * (radius_vestian_core**3)
mass_vestian_core = volume_vestian_core * density_metal
max_modeling_time = 24 * 10**6

alpha = 1.11
beta = -1.18
chi = -0.85
delta = 1680
epsilon = 487

surface_temp = 2000
surface_pressure = 0
therm_expansivity = 6 * (10**(-5))
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

core_conc_184w = []
core_conc_182w = []

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
    print("At time {} Myr".format(m.time / 10**6))

    # decay 182Hf to 182W
    m.decay_182hf()
    bulk_mass_182w.append(m.mass_bulk_182w)
    bulk_moles_182hf.append(m.moles_bulk_182hf)
    bulk_moles_182w.append(m.moles_bulk_182w)

    # fractionate some metal
    metal_added = m.fractionate_metal_exponential()
    metal_added_at_time.append(m.metal_mass_added)
    core_mass.append(m.core_mass)
    mantle_mass.append(m.mantle_mass)

    # calculate partitioning
    t = m.adiabat()
    p = m.hydrostat()
    d = m.partitioning(temperature=t, pressure=p)
    cmb_depth.append(m.core_mantle_boundary_depth)

    if m.time <= m.max_time_core_formation:
        # calculate equilibrium exchange
        mass_exchange_182w = m.exchange_mass(D=d, conc_w_silicate=m.conc_silicate_182w)
        mass_exchange_184w = m.exchange_mass(D=d, conc_w_silicate=m.conc_silicate_184w)
    else:
        mass_exchange_182w = 0
        mass_exchange_184w = 0

        # equilibrate
    m.set_concentrations(exchange_mass_182w=mass_exchange_182w, exchange_mass_184w=mass_exchange_184w)
    mass_182w_added_to_core = m.conc_core_bound_metal_182w * m.metal_mass_added
    mass_184w_added_to_core = m.conc_core_bound_metal_184w * m.metal_mass_added
    mass_182w_in_core = m.mass_core_182w +




fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.axhline(mass_vestian_core, linestyle="--", linewidth=1.0, color='red', label="Modern Vesta Core Mass")
ax1.plot([i / (10**6) for i in time_range], metal_added_at_time, linewidth=1.5, color='black', label="Vestian Core Mass")
ax1.grid()
ax1.set_xlabel("Time (Ma)")
ax1.set_ylabel("Core Mass Added at Time (kg)")
ax1.set_title("Exponential Core Growth Model on Vesta")

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# ax2.plot([i / (10**6) for i in time_range], bulk_mass_182w, linewidth=1.5, color='black', label="182W Bulk Mass")
# ax2.plot([i / (10**6) for i in time_range], core_mass_182w, linewidth=1.5, color='red', label="182W Core Mass")
# ax2.plot([i / (10**6) for i in time_range], mantle_mass_182w, linewidth=1.5, color='blue', label="182W Mantle Mass")
# ax2.grid()
# ax2.set_xlabel("Time (Ma)")
# ax2.set_ylabel("Bulk Mass 182W (kg)")
# ax2.set_title("Bulk Mass 182W in Vesta")

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot([i / (10**6) for i in time_range], bulk_moles_182w, linewidth=1.5, color='black', label="182W Bulk Moles")
ax3.plot([i / (10**6) for i in time_range], bulk_moles_182hf, linewidth=1.5, color='red', label="182Hf Bulk Moles")
ax3.grid()
ax3.set_xlabel("Time (Ma)")
ax3.set_ylabel("Bulk Moles")
ax3.set_title("Bulk Moles in Vesta")
ax3.legend(loc='lower right')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot([i / (10**6) for i in time_range], [i / 1000 for i in cmb_depth], linewidth=1.5, color='black', label="CMB Depth")
ax4.axhline(radius_vestian_core / 1000, linestyle="--", color='red', label='Modern Core Radius')
ax4.grid()
ax4.set_xlabel("Time (Ma)")
ax4.set_ylabel("CMB Depth (km)")
ax4.set_title("CMB Depth in Vesta")
ax4.legend(loc='lower right')

# fig5 = plt.figure()
# ax5 = fig5.add_subplot(111)
# ax5.plot([i / (10**6) for i in time_range], [i * 10**9 for i in core_conc_184w], linewidth=1.5, color='blue', label="184W")
# ax5.plot([i / (10**6) for i in time_range], [i * 10**9 for i in core_conc_182w], linewidth=1.5, color='red', label="182W")
# ax5.grid()
# ax5.set_xlabel("Time (Ma)")
# ax5.set_ylabel("W Isotope Concentration (ppb)")
# ax5.set_title("W Isotope Concentration in Vesta")
# ax5.legend(loc='lower right')
#
# fig6 = plt.figure()
# ax6 = fig6.add_subplot(111)
# ax6.axhline(mass_vestian_core, linestyle="--", linewidth=1.0, color='red', label="Modern Vesta Core Mass")
# ax6.plot([i / (10**6) for i in time_range], core_mass, linewidth=1.5, color='black', label="Vestian Core Mass")
# ax6.grid()
# ax6.set_xlabel("Time (Ma)")
# ax6.set_ylabel("Core Mass (kg)")
# ax6.set_title("Exponential Core Growth Model on Vesta")

plt.show()
