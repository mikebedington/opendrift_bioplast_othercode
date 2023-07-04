import numpy as np
import logging; logger = logging.getLogger(__name__)

from dotmap import DotMap


class DummyClass():
    def __init__(self, no_elements = 2, env_dict = {'sea_water_temperature':17, 'sea_water_salinity':35,
               'mass_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water':0.001}, ele_dict = None):

        self.no_elements = no_elements
        self.environment = DotMap()
        for this_key, this_val in env_dict.items():
            self.environment[this_key] = np.ones(no_elements)*this_val

        self.elements = DotMap()
        self.config = {}
        self.time_step = DotMap()
        self.time_step.total_seconds = 3600

        if ele_dict is not None:
            for this_var, this_data in ele_dict.items():
                self.elements[this_var] = np.ones(self.no_elements)*this_data

        self.debug = DotMap()

        self.CONFIG_LEVEL_ADVANCED = 4

    def _add_variables(self, var_data):
        for this_row in var_data:
            this_var = this_row[0]
            if this_var not in self.elements:
                self.elements[this_var] = np.ones(self.no_elements)*this_row[1]['default']

    def _add_config(self, conf_dict):
        for this_var, this_var_data in conf_dict.items():
            self.config[this_var] = this_var_data['default']    

    def set_config(self, conf_var, conf_val):
        self.config[conf_var] = conf_val

    def get_config(self, conf_var):
        return self.config[conf_var]

    def append_history(self, hist_dict):
        for this_key, this_item in self.elements.items():
            hist_dict[this_key].append(this_item)

        for this_key, this_item in self.debug.items():
            hist_dict[this_key].append(this_item)

        return hist_dict

# Defining the properties for plastic capable of biofouling
class BioPlast(DummyClass):
    """Extending Lagrangian3DArray with specific properties for biofoulable plastic
    """
    def __init__(self, *args, **kwargs):
        super(BioPlast, self).__init__(*args, **kwargs)

        self._add_variables([
        ('unfouled_diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # 
        ('unfouled_density', {'dtype': np.float32,
                                       'units':'kg/m^3',
                                       'default': self.sea_water_density(T=10, S=35)}),  # 
        ('biofilm_no_attached_algae', {'dtype': np.float32,
                                       'units': '',
                                       'default': 0}),
        ('total_density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': self.sea_water_density(T=10, S=35)}),
        ('total_diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),
        ('terminal_velocity', {'dtype': np.float32,
                      'units':'m/s',
                      'default': 0})])

        # Parameters for biofilm module
        self._add_config({
            'biofilm:biofilm_density':{'type':'float', 'default':1388.,
                                'min': None, 'max': None, 'units':'kg/m^3',
                                'description': 'Density (rho) of material collected on plastic',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:grazing_rate':{'type':'float', 'default':0.39,
                                'min': None, 'max': None, 'units': 'd-1',
                                'description': 'Grazing rate on biofilm',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:respiration_rate':{'type':'float', 'default':0.1,
                                'min': None, 'max': None, 'units': 'd-1',
                                'description': 'Respiration rate of algae on biofilm',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:temperature_coeff_respiration':{'type':'float', 'default':2,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Q10',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:algal_cell_volume':{'type':'float', 'default':2*10**-16,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Va',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:algal_cell_radius':{'type':'float', 'default':7.78*10**-6,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Ra',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:algal_carbon_per_cell':{'type':'float', 'default':2726*10**-9,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Used to convert phytoplankton C to no of cells',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'biofilm:shear':{'type':'float', 'default':1.7*10**-5,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Shear rate',
                                'level': self.CONFIG_LEVEL_ADVANCED}})

    def get_seawater_viscosity(self):
        return 0.001*(1.7915 - 0.0538*self.environment.sea_water_temperature+ 0.007*(self.environment.sea_water_temperature**(2.0)) - 0.0023*self.environment.sea_water_salinity)

    def update_biofilm(self):
        # This applies the 4 terms in eqn 11 of Kooi 2017

        ### Collision of algae with particle
        # Encounter rate is sum of brownian motion, differential settling, and shear
        r_a = self.get_config('biofilm:algal_cell_radius')
        algae_conc_cells = self.environment.mass_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water / self.get_config('biofilm:algal_carbon_per_cell')
        particle_surface_area = np.pi * self.elements.total_diameter**2

        d_brown = 4*np.pi*(self.get_diffusivity(self.elements.total_diameter/2) + self.get_diffusivity(r_a))*((self.elements.total_diameter/2) + r_a)

        d_settle = 0.5*np.pi*((self.elements.total_diameter/2)**2)*-self.elements.terminal_velocity

        d_shear = 1.3 * self.get_config('biofilm:shear') * ((self.elements.total_diameter/2) + r_a)**3

        d_tot = d_brown + d_settle + d_shear

        collision_term = (d_tot*algae_conc_cells)/particle_surface_area

        ### Growth on particle 
        # light and temp limited, needs light model and growth curve
        #growth = optimal_growth(light, temperature) - possible light model in sealice model
        growth = 0

        ### Grazing on particle
        grazing = self.get_config('biofilm:grazing_rate')

        ### Respiration
        # This is temperature dependent
        respiration = self.get_config('biofilm:temperature_coeff_respiration')**((self.environment.sea_water_temperature - 20)/10) * self.get_config('biofilm:respiration_rate')

        # Add up and convert to thickness
        day_seconds = 60*60*24 # To convert grazing and respiration from daily to per second
        A = self.elements.biofilm_no_attached_algae
        newAttached = A + (collision_term + growth*A - (grazing/day_seconds)*A - (respiration/day_seconds)*A)* self.time_step.total_seconds

        debug_dict = {'d_brown':d_brown, 'd_settle':d_settle, 'd_shear':d_shear, 'growth':np.ones(self.no_elements)*growth,
                        'grazing':np.ones(self.no_elements)*grazing, 'respiration':respiration, 'A':A}
        for this_key, this_val in debug_dict.items():
            self.debug[this_key] = this_val

        self.elements.biofilm_no_attached_algae = newAttached

    def get_diffusivity(self, radius):
        boltzmann = 1.0306*10**-13
        seawater_visc = self.get_seawater_viscosity()

        return boltzmann*(self.environment.sea_water_temperature) / (6*np.pi*seawater_visc*radius)

    def update_density(self):
        biofilm_volume = self.elements.biofilm_no_attached_algae*self.get_config('biofilm:algal_cell_volume')
        original_volume = (1/6)*np.pi*(self.elements.unfouled_diameter**3)

        new_volume = biofilm_volume + original_volume

        original_mass = original_volume * self.elements.unfouled_density

        self.elements.total_diameter = (6*new_volume/np.pi)**(1/3)
        self.elements.total_density = (biofilm_volume*self.get_config('biofilm:biofilm_density')+ original_mass) / new_volume

    def update_terminal_velocity(self):
        """Calculate terminal velocity for the plastic particle

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Get the properies which determine buoyancy
        partsize = self.elements.total_diameter

        DENSw = self.sea_water_density(T=self.environment.sea_water_temperature, S=self.environment.sea_water_salinity)
        DENSpart = self.elements.total_density
        dr = DENSw-DENSpart  # density difference

        # water viscosity
        my_w = self.get_seawater_viscosity()

        # terminal velocity for low Reynolds numbers
        W = (1.0/my_w)*(1.0/18.0)*g*partsize**2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W*1000*partsize/my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * self.environment.sea_water_temperature)  # in cm2/s
        d0 = (partsize * 100) - 0.4 * \
            (9.0 * my_w**2 / (100 * g) * DENSw / dr)**(1.0 / 3.0)  # cm
        W2 = 19.0*d0*(0.001*dr)**(2.0/3.0)*(my_w*0.001*DENSw)**(-1.0/3.0)
        # cm/s
        W2 = W2/100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    @staticmethod
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

