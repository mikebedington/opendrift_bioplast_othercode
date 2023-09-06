import numpy as np
from datetime import datetime, timedelta
from opendrift.readers.reader_netCDF_CF_generic import Reader as reader_nc
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.bioplast import BioPlastDrift
import glob as gb


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


forcing_data = gb.glob('./get_data/*.nc')

o = BioPlastDrift(loglevel=20)  # Set loglevel to 0 for debug information

for this_file in forcing_data:
	o.add_reader(reader_nc(this_file))
o.set_config('general:use_auto_landmask', True)

starting_diameter = 0.0014
starting_dense = sea_water_density(T=10,S=35)

o.seed_elements(lon=-12, lat=40, z=-10, number=1000, radius=30000, time=datetime(2022,6,1,12), unfouled_diameter=starting_diameter, unfouled_density=starting_dense, total_diameter=starting_diameter, total_density=starting_dense)

# Adjusting some configuration
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
# Vertical mixing requires fast time step
o.set_config('vertical_mixing:timestep', 60.) # seconds

# Running model
o.run(duration=timedelta(hours=120), time_step=3600)

# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
o.animation(fast=True, color='z')
o.plot_property('z', mean=True)


