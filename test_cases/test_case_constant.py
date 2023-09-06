import numpy as np
from datetime import datetime, timedelta
from opendrift.readers.reader_constant import Reader as ConstantReader
from opendrift.readers.reader_constant_2d import Reader as ConstantReader2D
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.bioplast import BioPlastDrift

def sea_water_density(T=10., S=35.):
	'''The function gives the density of seawater at one atmosphere
	pressure as given in :

	N.P. Fofonoff and R.C. Millard Jr.,1983,
	Unesco technical papers in marine science no. 44.

	S   = Salinity in promille of the seawater
	T   = Temperature of the seawater in degrees Celsius
	'''

	if np.atleast_1d(T).max() > 100:
		raise ValueError('Temperature should be in celcius, but is > 100')

	R4 = 4.8314E-04
	DR350 = 28.106331

	# Pure water density at atmospheric pressure
	# Bigg P.H. (1967) BR. J. Applied Physics pp.:521-537

	R1 = ((((6.536332E-09 * T - 1.120083E-06) * T + 1.001685E-04) *
		  T - 9.095290E-03) * T + 6.793952E-02) * T - 28.263737

	# Seawater density at atmospheric pressure
	# coefficients involving salinity :

	R2 = (((5.3875E-09 * T - 8.2467E-07) * T + 7.6438E-05) *
		  T - 4.0899E-03) * T + 8.24493E-01

	R3 = (-1.6546E-06*T+1.0227E-04)*T-5.72466E-03

	# International one-atmosphere equation of state of seawater :

	SIG = R1 + (R4*S + R3*np.sqrt(S) + R2)*S
	Dens0 = SIG + DR350 + 1000.
	return Dens0

def get_seawater_viscosity(T=10,S=35):
        return 0.001*(1.7915 - 0.0538*T+ 0.007*(T**(2.0)) - 0.0023*S)

# No horizontal movement, here only investigating vertical mixing and swimming
r = ConstantReader(
        {'x_sea_water_velocity': 0.01, 'y_sea_water_velocity': 0.05, 'x_wind': 0, 'y_wind': 0,
         'sea_water_temperature': 10, 'sea_water_salinity': 35,
         'land_binary_mask': 0, 'ocean_vertical_diffusivity': .02,
		 'mass_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water':0.00001,
         'mass_concentration_of_phytoplankton_expressed_as_chl_in_sea_water':0.00001,
		 'surface_net_downward_radiative_flux':0.1})
o = BioPlastDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.add_reader(r)
o.set_config('general:use_auto_landmask', False)

starting_diameter = 0.0014
starting_dense = sea_water_density(T=10,S=35)

o.seed_elements(lon=3, lat=60.5, z=-10, number=1000, radius=30000, time=datetime.now(), unfouled_diameter=starting_diameter, unfouled_density=starting_dense, total_diameter=starting_diameter, total_density=starting_dense)

# Adjusting some configuration
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
# Vertical mixing requires fast time step
o.set_config('vertical_mixing:timestep', 60.) # seconds

# Running model
o.run(duration=timedelta(hours=48), time_step=3600)

# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
o.animation(fast=True, color='z')
o.plot_property('z', mean=True)


"""
#Static 2D current field
#=======================
from opendrift.readers.reader_constant_2d import Reader

#%%
# Constructing a static, rotating ocean current field, 
lon, lat = np.meshgrid(np.linspace(2,6,30), np.linspace(59,62,30))
lon0 = 4
lat0 = 60.5
u = -(lat-lat0)/np.sqrt((lon-lon0)**2 + (lat-lat0)**2)
v = (lon-lon0)/np.sqrt((lon-lon0)**2 + (lat-lat0)**2)
lon = np.linspace(0,5,30)
lat = np.linspace(59,62,30)

r = ConstantReader(x=lon, y=lat, proj4='+proj=latlong',
           array_dict = {'x_sea_water_velocity': u, 'y_sea_water_velocity': v, 'sea_water_temperature': 10,
         'land_binary_mask': 0, 'ocean_vertical_diffusivity': .02, 'sea_water_salinity': 35})

o = OceanDrift(loglevel=20)
o.set_config('environment:fallback:land_binary_mask', 0)
o.add_reader(r)
o.seed_elements(lon=3, lat=60.5, number=1000, radius=30000, time=datetime.now())
o.run(duration=timedelta(hours=72))
o.animation(fast=True)
"""

"""
# Forcing with Topaz ocean model and MEPS atmospheric model
o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
    'https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc'])
"""

