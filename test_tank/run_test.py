import numpy as np
from bioplast_sandbox import *

time_steps = 10
n_elements = 2

o = BioPlast(no_elements = n_elements)

store_history = {'unfouled_diameter':[], 'unfouled_density':[], 'biofilm_no_attached_algae':[], 'total_density':[],
        'total_diameter':[], 'terminal_velocity':[], 'd_brown':[np.zeros(n_elements)], 'd_settle':[np.zeros(n_elements)], 
        'd_shear':[np.zeros(n_elements)], 'growth':[np.zeros(n_elements)], 'grazing':[np.zeros(n_elements)],
        'respiration':[np.zeros(n_elements)], 'A':[np.zeros(n_elements)]}

o.append_history(store_history)
for i in np.arange(0, time_steps):	
	o.update_biofilm()
	o.append_history(store_history)
	o.update_density()
	o.update_terminal_velocity()
	o.append_history(store_history)

for this_var, this_data in store_history.items():
	store_history[this_var] = np.asarray(this_data)
