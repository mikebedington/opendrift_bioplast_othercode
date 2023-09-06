import configparser
import motuclient
import datetime as dt
from time import time

class MotuOptions:
    def __init__(self, attrs: dict):
        super(MotuOptions, self).__setattr__("attrs", attrs)

    def __setattr__(self, k, v):
        self.attrs[k] = v

    def __getattr__(self, k):
        try:
            return self.attrs[k]
        except KeyError:
            return None

cfg = configparser.ConfigParser()
cfg.read('motuclient-python.ini')

user = cfg.get('Main','user')
pwd = cfg.get('Main','pwd')

start_date = dt.datetime(2022,6,1)
end_date = dt.datetime(2022,7,1)
ll_box = [[-14,-10],[37,44]]

service_id = {'phys':'IBI_ANALYSISFORECAST_PHY_005_001-TDS', 'bio':'IBI_ANALYSISFORECAST_BGC_005_004-TDS'}

all_vars = {'phys':{'uo':'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m', 'vo':'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m',
    'thetao':'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m', 'so':'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m'},
	'bio':{'chl':'cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m', 'phyc':'cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m',
			'zooc':'cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m'}}

data_request_dict = {
    "service_id": service_id,
    "date_min": start_date.strftime('%Y-%m-%d %H:%M:%S'),
    "date_max": end_date.strftime('%Y-%m-%d %H:%M:%S'),
    "longitude_min": ll_box[0][0],
    "longitude_max": ll_box[0][1],
    "latitude_min": ll_box[1][0],
    "latitude_max": ll_box[1][1],
    "depth_min": 0,
    "depth_max": 5000,
    "motu": 'https://nrt.cmems-du.eu/motu-web/Motu',
    "out_dir": ".",
    "auth_mode": "cas",
    "user": user,
    "pwd": pwd
}

for this_type in all_vars.keys():
	data_request_dict['service_id'] = service_id[this_type]
	var_data = all_vars[this_type]
	for this_var, this_source in var_data.items():
		print(f'Getting {this_var}', flush=True)
		start_time = time()
		outfile = f'{this_var}_out.nc'
		data_request_dict['out_name'] = outfile
		data_request_dict['product_id'] = this_source
		data_request_dict['variable'] = [this_var]

		motuclient.motu_api.execute_request(MotuOptions(data_request_dict))
		print(f'Completed - {time() - start_time:0.2f}s', flush=True)

# The ibi output for phyto- and zooplankton are in mmol/m^3, need to convert to mg/m3
import netCDF4 as nc
import numpy as np
c_moleweight = 12.011

phyc = nc.Dataset('phyc_out.nc', 'r+')
c_grams_per_m3 = phyc['phyc'][:] * c_moleweight
new_var = phyc.createVariable('phyc_c', 'f4', ['time', 'depth', 'latitude', 'longitude']) 
new_var.setncattr('standard_name', 'mass_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water')
new_var.setncattr('unit', 'mg.m-3')
new_var[:] = c_grams_per_m3
phyc.close()

zooc = nc.Dataset('zooc_out.nc', 'r+')
c_grams_per_m3 = zooc['zooc'][:] * c_moleweight
new_var = zooc.createVariable('zooc_c', 'f4', ['time', 'depth', 'latitude', 'longitude'])
new_var.setncattr('standard_name', 'mass_concentration_of_zooplankton_expressed_as_carbon_in_sea_water')
new_var.setncattr('unit', 'mg.m-3')
new_var[:] = c_grams_per_m3
zooc.close()
	
