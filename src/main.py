import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import CoolProp.CoolProp as cp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd
from math import pi, log10, sin, radians, floor
from decimal import Decimal

safety_factor = None
volheating = None
hpshield_id = None
hpshield_thickness = None
hpshield_sector = None
fluid_mass_flowrate = None
pipe_od = None
pipe_thickness = None

# Get values
def btn_enter_data_clicked():
    global safety_factor, volheating, hpshield_id, hpshield_thickness, hpshield_sector, fluid_mass_flowrate, pipe_od, pipe_thickness

    safety_factor = (box_safety_factor.get())
    volheating = (box_volheating.get())
    hpshield_id = (box_hpshield_id.get())
    hpshield_thickness = (box_hpshield_thickness.get())
    hpshield_sector = (box_hpshield_sector.get())
    fluid_mass_flowrate = (box_fluid_mass_flowrate.get())
    pipe_od = (box_pipe_od.get())
    pipe_thickness = (box_pipe_thickness.get())

    if not all([safety_factor, volheating, hpshield_id, hpshield_thickness, hpshield_sector, fluid_mass_flowrate, pipe_od, pipe_thickness]):
        messagebox.showerror("Error", "All fields must have data")
        return

# Create GUI for calculations
# Create window
window = tk.Tk()
window.title("HPS Coolant Calculation")
window.geometry('400x350')
window.tk.call('tk', 'scaling', 2.0)

# Close window button
def btn_exit_clicked():
    window.destroy()

# Define validation function
def validate_input(new_value):
    # Allow numbers with decimals or an empty input
    return new_value.replace('.', '', 1).isdigit() or new_value == ""

# Register validation function
validate_cmd = window.register(validate_input)

# Create labels & input boxes
lbl_safety_factor = tk.Label(window, text="Safety Factor")
lbl_safety_factor.grid(column=0, row=0)
box_safety_factor = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_safety_factor.grid(column=1, row=0)

lbl_volheating = tk.Label(window, text="volheating")
lbl_volheating.grid(column=0, row=1)
box_volheating = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_volheating.grid(column=1, row=1)

lbl_hpshield_id = tk.Label(window, text="HPS Inner Diameter")
lbl_hpshield_id.grid(column=0, row=2)
box_hpshield_id = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_hpshield_id.grid(column=1, row=2)

lbl_hpshield_thickness = tk.Label(window, text="HPS Thickness")
lbl_hpshield_thickness.grid(column=0, row=3)
box_hpshield_thickness = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_hpshield_thickness.grid(column=1, row=3)

lbl_hpshield_sector = tk.Label(window, text="Pipe Sector")
lbl_hpshield_sector.grid(column=0, row=4)
box_hpshield_sector = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_hpshield_sector.grid(column=1, row=4)

lbl_fluid_mass_flowrate = tk.Label(window, text="Mass Flowrate")
lbl_fluid_mass_flowrate.grid(column=0, row=5)
box_fluid_mass_flowrate = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_fluid_mass_flowrate.grid(column=1, row=5)

lbl_pipe_od = tk.Label(window, text="Pipe Outer Diameter")
lbl_pipe_od.grid(column=0, row=6)
box_pipe_od = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_pipe_od.grid(column=1, row=6)

lbl_pipe_thickness = tk.Label(window, text="Pipe Thickness")
lbl_pipe_thickness.grid(column=0, row=7)
box_pipe_thickness = ttk.Entry(window, validate="key", validatecommand=(validate_cmd, "%P"))
box_pipe_thickness.grid(column=1, row=7)

#Calculate button
btn_calculate = tk.Button(window, text="Calculate", command=btn_enter_data_clicked)
btn_calculate.grid(column=1, row=8)

#Exit button
btn_exit = tk.Button(window, text="Exit", command=btn_exit_clicked)
btn_exit.grid(column=0, row=8)

# Add a button to display the graph
#btn_graph = tk.Button(window, text="Show Graph", command=display_graph)
#btn_graph.grid(column=1, row=10)

window.mainloop()

# Unit converters
celsius_kelvin = 273.15
pascal_bar = 1e5
mega = 1e6

# HP shield, pipes, and initial fluid inputs at cryostat inlet
# Safety factor is applied to volumetric heating values

safety_factor = 1
#print('box safety factor = ',safety_factor)
#print(type(safety_factor))
volheating = 5.57e6 # W/m3
volheating_verticaldist_file = 'hpshield_volheating_verticaldist_pi6.xlsx'

hpshield_id = 3.074
hpshield_thickness = 0.163
hpshield_sector = 6 # degrees
hpshield_height = 0.5 # block

fluid_name = 'HeavyWater'
fluid_temperature = 200 # degrees C
fluid_pressure = 85 # bar
fluid_mass_flowrate = 11.3 # kg/s

# Coolant pipe through HP shield
pipe_od = 0.0334 # DN25
pipe_thickness = 0.001651 # Sch 5
pipe_roughness = 45e-6
pipe_stagger = 'y' # y/n, Staggered pipes = not aligned toroidally, Non-staggered = aligned toroidally

# Pipe length for each pipe between cryostat interface and bottom of central column
pipe_length_cryostat = 3

# Pipe bend at the top of the central column
pipe_bend_radius = 0.06
kl_bend_180 = 0.2 # Loss coefficients from Fluid Mechanics Fundamentals and Applications - Cengel & Cimbala

# PropsSI uses SI Units https://en.wikipedia.org/wiki/International_System_of_Units 
# PropsSI available properties http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function 

class fluid:
    def __init__(self, name, temperature, pressure):
        self.name = name
        self.temperature = temperature
        self.pressure = pressure

    @property
    def properties(self):
        temperature_si = self.temperature + celsius_kelvin
        pressure_si = self.pressure * pascal_bar
        
        self.density = cp.PropsSI('D','T',temperature_si,'P',pressure_si,self.name) # kg/m3
        self.dynamic_viscosity = cp.PropsSI('V','T',temperature_si,'P',pressure_si,self.name) # Pa.s
        self.specific_heat = cp.PropsSI('C','T',temperature_si,'P',pressure_si,self.name) # J/kg.K
        self.thermal_conductivity = cp.PropsSI('CONDUCTIVITY','T',temperature_si,'P',pressure_si,self.name) # W/m.K
        self.pr_number = cp.PropsSI('PRANDTL','T',temperature_si,'P',pressure_si,self.name)
        self.phase = cp.PhaseSI('T',temperature_si,'P',pressure_si,self.name)
        self.temperature_saturation = cp.PropsSI('T','P',pressure_si,'Q',0,self.name) - celsius_kelvin # C

        return self

# Calculate all pipe parameters
# Friction factor equation https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
# Gnielinkski correlation for Nusselt number https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation

class pipe:
    def __init__(self, od, thickness, roughness, length, fluid, mass_flowrate):
        self.od = od # m
        self.thickness = thickness # m
        self.length = length # m
        self.roughness = roughness # m
        self.fluid = fluid
        self.fluid.properties # Initialise fluid properties
        self.mass_flowrate = mass_flowrate # kg/s
        self.heatflux = 0.0 # W/m2
        self.heatflux_flomaster = 0.0 # W/m2
        self.power = 0.0 # W

    @property
    def static_properties(self):
        self.id = self.od - (2 * self.thickness) # m
        self.id_wall_area = pi * self.id * self.length # m2
        self.id_cs_area = 0.25 * pi * (self.id ** 2) # m2
        self.relative_roughness = self.roughness / self.id

        return self

    # Properties that change every iteration
    @property
    def dynamic_properties(self):
        self.volumetric_flowrate = self.mass_flowrate / self.fluid.density # m3/s
        self.velocity = self.volumetric_flowrate / self.id_cs_area # m/s
        self.re_number = self.fluid.density * self.velocity * self.id / self.fluid.dynamic_viscosity
        self.f = (1 / (-1.8 * log10(((self.relative_roughness / 3.7) ** 1.11) + (6.9 / self.re_number)))) ** 2
        self.dP = (self.f * self.length * (self.fluid.density / 2) * ((self.velocity ** 2) / self.id)) / pascal_bar # bar

        return self

    @property
    def gnielinski(self):
        if ((self.re_number >= 3000 and self.re_number <= 5e6) and (self.fluid.pr_number >= 0.5 and self.fluid.pr_number <= 2000)):
            self.nu_number = ((self.f / 8) * (self.re_number - 1000) * self.fluid.pr_number) / (1 + 12.7 * ((self.f / 8) ** 0.5) * ((self.fluid.pr_number ** (2/3)) - 1))
            self.htc = self.nu_number * self.fluid.thermal_conductivity / self.id # W/m2.K

            return self

        else:
            print("Pipe and fluid parameters not valid for Gnielinski correlation")

# Initialisations
# Coolant_next contains the fluid properties at the end of each pipe section

coolant = fluid(fluid_name, fluid_temperature, fluid_pressure)
coolant_next = fluid(fluid_name, fluid_temperature, fluid_pressure)

hpshield_inlet_pipe = pipe(pipe_od, pipe_thickness, pipe_roughness, hpshield_height, coolant_next, fluid_mass_flowrate)
hpshield_outlet_pipe = pipe(pipe_od, pipe_thickness, pipe_roughness, hpshield_height, coolant_next, fluid_mass_flowrate)
cryostat_pipe = pipe(pipe_od, pipe_thickness, pipe_roughness, pipe_length_cryostat, coolant_next, fluid_mass_flowrate)

pipe_bend_length = 2 * pi * pipe_bend_radius / 2
bend_pipe = pipe(pipe_od, pipe_thickness, pipe_roughness, pipe_bend_length, coolant_next, fluid_mass_flowrate)

cumulative_pipe_length = 0

# Results storage with initial values

results = {
    'Location': ['cryo_inlet'],
    'Total Length (m)': [cumulative_pipe_length],
    'Pressure (bar)': [coolant_next.pressure],
    'Temperature (C)': [coolant_next.temperature],
    'HTC (W/m2.K)': [0],
    'Input Pipe Power (W)': [0],
    'Input Pipe Heat Flux (W/m2)': [0]
}

# Results appender

def append_results(height, total_length, pressure, temperature, htc, pipe_power, pipe_heatflux_flomaster):
    results['Location'].append(height)
    results['Total Length (m)'].append(total_length)
    results['Pressure (bar)'].append(pressure)
    results['Temperature (C)'].append(temperature)
    results['HTC (W/m2.K)'].append(htc)
    results['Input Pipe Power (W)'].append(pipe_power)
    results['Input Pipe Heat Flux (W/m2)'].append(pipe_heatflux_flomaster)

# Calculate key dimensions, areas, volumens
# Calculate heating input due to volumetric loads
# Calculate equivalent heat flux into channel
# Outer ratio is worst case

hpshield_od = hpshield_id + (2 * hpshield_thickness)
hpshield_mid = (hpshield_od + hpshield_id) / 2

if pipe_stagger == 'y':
    sector_area_inlet = 0.25 * pi * ((hpshield_mid ** 2) - (hpshield_id ** 2)) * (hpshield_sector / 360)
    sector_vol_inlet = sector_area_inlet * hpshield_height # m3
    volheating_inlet_pipe = 2 * volheating / 2.5 # W/m3
    height_factor_inlet = 'Outer Ratio'

    sector_area_outlet = 0.25 * pi * ((hpshield_od ** 2) - (hpshield_mid ** 2)) * (hpshield_sector / 360)
    sector_vol_outlet = sector_area_outlet * hpshield_height # m3
    volheating_outlet_pipe = volheating_inlet_pipe * 1.5 # W/m3
    height_factor_outlet = 'Outer Ratio'

    coolant_percentage = (hpshield_inlet_pipe.static_properties.id_cs_area + hpshield_outlet_pipe.static_properties.id_cs_area) / (sector_area_inlet + sector_area_outlet)

elif pipe_stagger == 'n':
    sector_area_half = 0.5 * 0.25 * pi * ((hpshield_od ** 2) - (hpshield_id ** 2)) * (hpshield_sector / 360)

    sector_vol_inlet = sector_area_half * hpshield_height # m3
    volheating_inlet_pipe = volheating # W/m3
    height_factor_inlet = 'Outer Ratio'

    sector_vol_outlet = sector_vol_inlet # m3
    volheating_outlet_pipe = volheating_inlet_pipe # W/m3
    height_factor_outlet = 'Outer Ratio'

    coolant_percentage = (hpshield_inlet_pipe.static_properties.id_cs_area + hpshield_outlet_pipe.static_properties.id_cs_area) / (2 * sector_area_half)

hpshield_inlet_pipe.heatflux = (volheating_inlet_pipe * sector_vol_inlet) / hpshield_inlet_pipe.static_properties.id_wall_area # W/m2, for one pipe
hpshield_outlet_pipe.heatflux = (volheating_outlet_pipe * sector_vol_outlet) / hpshield_outlet_pipe.static_properties.id_wall_area # W/m2, for one pipe

# Calculate temperature, and pressure change through each pipe loop section
# Heights increments are used according to 3D neutronics analysis from PI-4
# Assumes cold inlet pipe on inboard side of HP shield for staggered case
# Volumetric heating vertical distributions are applied

volheating_verticaldist = pd.read_excel(volheating_verticaldist_file)
cc_bottom = volheating_verticaldist['Location'].min()
cc_top = volheating_verticaldist['Location'].max()
height_increment = int(hpshield_height * 100)

# 1. Pipe from cryostat inlet to central column inlet

coolant_next.temperature = ((safety_factor * 0 * cryostat_pipe.static_properties.id_wall_area) / (cryostat_pipe.mass_flowrate * coolant_next.properties.specific_heat)) + coolant_next.temperature
coolant_next.pressure -= cryostat_pipe.dynamic_properties.dP
cumulative_pipe_length += cryostat_pipe.length

append_results("cc_inlet", cumulative_pipe_length, coolant_next.pressure, coolant_next.temperature, cryostat_pipe.gnielinski.htc, cryostat_pipe.power, cryostat_pipe.heatflux_flomaster)

# 2. HP shield inlet pipe

for h in range(cc_bottom, (cc_top + 1), height_increment):
    row = volheating_verticaldist[volheating_verticaldist['Location'] == h].index
    height_factor = volheating_verticaldist.loc[row[0],height_factor_inlet]

    coolant_next.temperature = ((safety_factor * height_factor * hpshield_inlet_pipe.heatflux * hpshield_inlet_pipe.static_properties.id_wall_area) / (hpshield_inlet_pipe.mass_flowrate * coolant_next.properties.specific_heat)) + coolant_next.temperature
    coolant_next.pressure -= hpshield_inlet_pipe.dynamic_properties.dP

    hpshield_inlet_pipe.power = height_factor * volheating_inlet_pipe * sector_vol_inlet
    hpshield_inlet_pipe.heatflux_flomaster = height_factor * hpshield_inlet_pipe.heatflux

    cumulative_pipe_length += hpshield_inlet_pipe.length

    append_results(h, cumulative_pipe_length, coolant_next.pressure, coolant_next.temperature, hpshield_inlet_pipe.gnielinski.htc, hpshield_inlet_pipe.power, hpshield_inlet_pipe.heatflux_flomaster)

# 3. Pipe bend at top of HP shield between inlet and outlet pipes

coolant_next.temperature = ((safety_factor * 0 * bend_pipe.static_properties.id_wall_area) / (bend_pipe.mass_flowrate * coolant_next.properties.specific_heat)) + coolant_next.temperature
bend_dP = 0.5 * kl_bend_180 * coolant_next.density * (bend_pipe.dynamic_properties.velocity ** 2) / pascal_bar
coolant_next.pressure -= bend_dP
cumulative_pipe_length += bend_pipe.length

append_results("pipe_bend", cumulative_pipe_length, coolant_next.pressure, coolant_next.temperature, bend_pipe.gnielinski.htc, bend_pipe.power, bend_pipe.heatflux_flomaster)

# 4. HP shield outlet pipe

for h in range((cc_top - height_increment), (cc_bottom - 1), -height_increment):
    row = volheating_verticaldist[volheating_verticaldist['Location'] == h].index
    height_factor = volheating_verticaldist.loc[row[0],height_factor_outlet]

    coolant_next.temperature = ((safety_factor * height_factor * hpshield_outlet_pipe.heatflux * hpshield_outlet_pipe.static_properties.id_wall_area) / (hpshield_outlet_pipe.mass_flowrate * coolant_next.properties.specific_heat)) + coolant_next.temperature
    coolant_next.pressure -= hpshield_outlet_pipe.dynamic_properties.dP

    hpshield_outlet_pipe.power = height_factor * volheating_outlet_pipe * sector_vol_outlet
    hpshield_outlet_pipe.heatflux_flomaster = height_factor * hpshield_outlet_pipe.heatflux

    cumulative_pipe_length += hpshield_outlet_pipe.length

    append_results(h, cumulative_pipe_length, coolant_next.pressure, coolant_next.temperature, hpshield_outlet_pipe.gnielinski.htc, hpshield_outlet_pipe.power, hpshield_outlet_pipe.heatflux_flomaster)

# 5. Pipe from central column outlet to cryostat outlet

coolant_next.temperature = ((safety_factor * 0 * cryostat_pipe.static_properties.id_wall_area) / (cryostat_pipe.mass_flowrate * coolant_next.properties.specific_heat)) + coolant_next.temperature
coolant_next.pressure -= cryostat_pipe.dynamic_properties.dP
cumulative_pipe_length += cryostat_pipe.length

append_results("cc_inlet", cumulative_pipe_length, coolant_next.pressure, coolant_next.temperature, cryostat_pipe.gnielinski.htc, cryostat_pipe.power, cryostat_pipe.heatflux_flomaster)

# Table of results
# Export results to excel
results_df = pd.DataFrame(results)
results_df = results_df.set_index('Location')
print(results_df)

# Export results to excel
results_df.to_excel("hpshield_thermal_output.xlsx")

# Graph plot
def display_graph():
    fig = Figure(figsize=(3, 2), dpi=100)
    ax = fig.add_subplot(111)

    # Generate results diagram
    fig, ax1 = plt.subplots()
    plt.title("Temperature & Pressure vs Pipe Length")

    ax1.plot(results_df['Total Length (m)'], results_df['Temperature (C)'], color = 'r')
    ax2 = ax1.twinx()
    ax2.plot(results_df['Total Length (m)'], results_df['Pressure (bar)'], color = 'b', linestyle = 'dashed')

    ax1.set_xlabel("Total Length (m)")
    ax1.set_ylabel("Temperature (C)", color = 'r')
    ax2.set_ylabel("Pressure (bar)", color = 'b')

    ax1.grid(axis = 'x')
    plt.show()


# Length of chord calculation to check whether pipe fits in sector
sector_width_min = hpshield_id * sin(radians(hpshield_sector) / 2) 
if ((2 * pipe_od) < (sector_width_min + 0.05)):
    print("Pipe fits within sector")
else:
    print("Pipe may not fit")

print(f"Coolant is {coolant_percentage * 100:.1f}% of volume")
print(f"Inlet temperature = {coolant.temperature:.1f}°C")
print(f"Outlet temperature = {coolant_next.temperature:.1f}°C")
print(f"Inlet pressure = {coolant.pressure:.1f} bar")
print(f"Outlet pressure = {coolant_next.pressure:.1f} bar")
print(f"Pressure drop = {coolant.pressure - coolant_next.pressure:.1f} bar")
print(f"Final coolant phase = {coolant_next.phase}")
print(f"Saturation temperature at outlet pressure = {floor(coolant_next.properties.temperature_saturation)}°C")

final_velocity = cryostat_pipe.dynamic_properties.velocity
print(f"Fluid velocity = {final_velocity:.1f} m/s")

# Check power balance, Q (volumetric heating) = m * Cp * dT (power extracted by fluid)
total_power_volumetric = fluid_mass_flowrate * ((coolant.properties.specific_heat + coolant_next.properties.specific_heat) / 2) * (coolant_next.temperature - coolant.temperature)
total_power_fluid = results_df['Input Pipe Power (W)'].sum()
difference_percentage = abs(1 - (total_power_fluid / total_power_volumetric))
print(f"Total power into full HP shield = {total_power_volumetric * (360 / hpshield_sector) / mega:.1f} MW")
print(f"Total power extracted by coolant = {total_power_fluid * (360 / hpshield_sector) / mega:.1f} MW")
print(f"Total power discrepancy = {difference_percentage * 100:.1f}%")



