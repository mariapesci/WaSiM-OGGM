# Script for running OGGM with default dataset CRU

# # All glaciers in Gepatsch valid for the RGI 2003
# - run the simulations from 1969:
#        a) no initialization or spin-up
#        b) with dynamic spin-up
#        c) with initialization method (Eis et al., 2021)
# - use the climate data from CRU
# - use the glacier's outline from OGGM (RGI 6.0) for the year 2003

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import oggm
import os
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import inversion
from oggm.core import massbalance, climate
import xarray as xr
import seaborn as sns
import shapely.geometry as shpg


# Definition of parameters and paths for OGGM
cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border']= 160
cfg.PARAMS['baseline_climate'] = ''
cfg.CONFIG_FILE
cfg.PARAMS['use_multiprocessing'] = True

cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU' 

# Define the glaciers that are going to be used for the runs
fr = utils.get_rgi_region_file(11, version='62')
gdf = gpd.read_file(fr)
catchment = gpd.read_file(..)
in_bas = [catchment.geometry.contains(shpg.Point(x,y))[0] for
          (x,y) in zip(gdf.CenLon, gdf.CenLat)]
gdf_sel = gdf.loc[in_bas]
ax = catchment.plot();
gdf_sel.plot(ax=ax, edgecolor='k')

gdirs = workflow.init_glacier_directories(gdf_sel.RGIId, from_prepro_level=3, reset=True, force=True)
print(gdirs[0].rgi_date)

list_talks = [
         tasks.glacier_masks,
         tasks.compute_centerlines,
         tasks.initialize_flowlines,
         tasks.compute_downstream_line,
         tasks.catchment_area,
         tasks.catchment_width_geom,
         tasks.catchment_width_correction,
         tasks.compute_downstream_bedshape
         ]
for task in list_talks:
    workflow.execute_entity_task(task, gdirs)

from oggm.core import massbalance, climate
params_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'
workflow.download_ref_tstars(base_url=params_url)

# Now calibrate
cfg.PARAMS['mu_star_halfperiod'] = 3
workflow.execute_entity_task(tasks.local_t_star, gdirs);
workflow.execute_entity_task(tasks.mu_star_calibration, gdirs);

for i in range(len(gdirs)):
    print(gdirs[i].read_json('local_mustar')['t_star'])
    print(gdirs[i].read_json('local_mustar')['mu_star_glacierwide'])
    

# Compute the glacier mass balance:
years = np.arange(1969, 2020)
from oggm.core.massbalance import MultipleFlowlineMassBalance
for i in range(len(gdirs)):
    mbmod = MultipleFlowlineMassBalance(gdirs[i], use_inversion_flowlines=True)
    mb_ts = mbmod.get_specific_mb(year=years) # gets the annual mass balance
    plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)');
    mb = pd.DataFrame(index=years, data=mb_ts, columns={'mb':mb_ts})
    mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU/MBcru_'+gdirs[i].rgi_id+'.csv')


# Inversion
list_talks = [
         tasks.prepare_for_inversion,  # This is a preprocessing task
         tasks.mass_conservation_inversion,  # This does the actual job (glacier thickness along flowlines)
         tasks.filter_inversion_output  # This smoothes the thicknesses at the tongue a little
         ]
for task in list_talks:
    workflow.execute_entity_task(task, gdirs)

graphics.plot_inversion(gdirs, figsize=(8, 7))


# Convert the flowlines to a "glacier" for the ice dynamics module
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs);

workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, ys=1969,
                              #fixed_geometry_spinup_yr=gdirs[0].rgi_date,
                              climate_filename='climate_historical')

ds_run = utils.compile_run_output(gdirs)
ds_run

# Plot volume and length evolution of all the glaciers:
area = []
for i in range(len(gdirs)):
    area_i = gdirs[i].rgi_area_m2
    area.append(area_i)
  
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
ds_run.sum(dim='rgi_id').volume.plot.line(ax=ax1, hue='rgi_id');
ds_run.sum(dim='rgi_id').area.plot.line(ax=ax2, hue='rgi_id');
#ax2.set_ylim(2e7,7e7)
ax2.scatter(1969, 61359000, color='r');
ax2.scatter(1998, 52286000, color='r');
ax2.scatter(gdirs[0].rgi_date, sum(area), color='g')
ax2.scatter(2006, 49650000, color='r');
ax2.scatter(2015, 42566000, color='r');