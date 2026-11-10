import openmc
import numpy as np


#volume calculation

PI = np.pi
r_fuel = 1.5
R = input("Enter F/M ratio: ")
mod_pitch = np.sqrt(((PI*r_fuel**2)/R)+PI*r_fuel**2)
rod_height = 30.0

v_fuel = PI*r_fuel**2*rod_height
v_mod = mod_pitch**2*rod_height - v_fuel
fm_ratio = v_fuel/v_mod

#-----
##MATERIALS
#-----

openmc.config['cross_sections'] = "endfb80/endfb-viii.1-hdf5/cross_sections.xml"

fuel = openmc.Material(name="Uranium Fuel")
fuel.add_nuclide('U235',0.1975, 'wo')
fuel.add_nuclide('U238',0.8025, 'wo')
fuel.set_density('g/cm3',10.97)

mod = openmc.Material(name="Be2C Moderator")
mod.add_element('Be',2.0)
mod.add_element('C',1.0)
mod.add_s_alpha_beta('c_Be_in_Be2C')
mod.add_s_alpha_beta('c_C_in_Be2C')
mod.set_density('g/cm3',1.9)

#export materials

materials = openmc.Materials([fuel,mod])
materials.cross_sections = openmc.config['cross_sections'] 
materials.export_to_xml()

#-----
##GEOMETRY
#-----

#surfaces

fuel_surface = openmc.ZCylinder(0.0,0.0,r=r_fuel)
left_mod = openmc.XPlane(x0=-mod_pitch/2,boundary_type='reflective')
right_mod = openmc.XPlane(x0=mod_pitch/2,boundary_type='reflective')
back_mod = openmc.YPlane(y0=mod_pitch/2,boundary_type='reflective')
front_mod = openmc.YPlane(y0=-mod_pitch/2,boundary_type='reflective')
top_mod = openmc.ZPlane(z0=rod_height/2,boundary_type='vacuum')
bot_mod = openmc.ZPlane(z0=-rod_height/2,boundary_type='vacuum')


#regions

fuel_region = -fuel_surface & -top_mod & +bot_mod
mod_region = +fuel_surface & +left_mod & -right_mod & -back_mod & +front_mod & -top_mod & +bot_mod

#cells

fuel_cell = openmc.Cell(fill=fuel,region=fuel_region)
mod_cell = openmc.Cell(fill=mod,region=mod_region)

#universe

universe = openmc.Universe(cells=[fuel_cell,mod_cell])

#geometry add

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

#-----
##SETTINGS
#-----

settings = openmc.Settings()
settings.particles=1000000
settings.inactive = 5
settings.batches=15

#source

r_dist = openmc.stats.PowerLaw(a=0.0, b=r_fuel, n=1.0)

phi_dist = openmc.stats.Uniform(a=0.0, b=2*np.pi)

z_dist =  openmc.stats.Uniform(a=-rod_height/2, b=rod_height/2)


settings.run_mode='eigenvalue'
source = openmc.IndependentSource()

source.space = openmc.stats.CylindricalIndependent(
    r=r_dist, 
    phi=phi_dist, 
    z=z_dist, 
    origin=(0.0, 0.0, 0.0)
)

source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([10.0e6],[1.0])
source.time = openmc.stats.Uniform(0,1e-6)
settings.source = source

#export settings
settings.output = {'summary': True, 'tallies': True, 'path': 'results'}
settings.export_to_xml()

#-----
##TALLIES
#------

flux_tally = openmc.Tally(name="flux")
flux_tally.scores = ['flux']

tallies = openmc.Tallies([flux_tally])
tallies.export_to_xml()

#-----
##RUN
#-----
print(fm_ratio)
openmc.run()



