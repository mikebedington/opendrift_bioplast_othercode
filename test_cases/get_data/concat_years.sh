#!/bin/bash

for i in {2010..2020}
do
   ncks --mk_rec_dmn time era_data_${i}.nc part_${i}.nc
done

ncrcat part_*.nc era_out.nc

rm part_*.nc
