#!/usr/bin/env python
# coding: utf-8
'''
# # All glaciers in Gepatsch valid for the RGI 2003
# - run the simulations from 1969:
        a) no initialization or spin-up
        b) with dynamic spin-up
        c) with initialization method (Eis et al., 2021)
# - use the climate data from CRU
# - use the glacier's outline from OGGM (RGI 6.0) for the year 2003

'''
# In[1]:
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

#%%
# Definition of parameters and paths for OGGM
cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border']= 160
cfg.PARAMS['baseline_climate'] = ''
cfg.CONFIG_FILE
cfg.PARAMS['use_multiprocessing'] = True

# Gepatsch
cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU' 
# Rofenache
#cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU'

init_glac = 'dyn_spinup'
cal_geod_mb = False

#%%
# Define the glaciers that are going to be used for the runs
# in this case, all the valid glaciers in Gepatsch for the year 2003
fr = utils.get_rgi_region_file(11, version='62')
gdf = gpd.read_file(fr)
# Gepatsch
catchment = gpd.read_file('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/shp/catchment_alm.shp')
# Rofenache
#catchment = gpd.read_file('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/rofenache_total.shp')
in_bas = [catchment.geometry.contains(shpg.Point(x,y))[0] for
          (x,y) in zip(gdf.CenLon, gdf.CenLat)]
gdf_sel = gdf.loc[in_bas]
# ax = catchment.plot();
# gdf_sel.plot(ax=ax, edgecolor='k')

#%%
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

#%%
from oggm.core import massbalance, climate
if cal_geod_mb:
    workflow.execute_entity_task(tasks.mu_star_calibration_from_geodetic_mb, gdirs,
                              ignore_hydro_months=True, min_mu_star=80, max_mu_star=300)
    for i in range(len(gdirs)):
        print(gdirs[i].read_json('local_mustar')['t_star'])
        print(gdirs[i].read_json('local_mustar')['mu_star_glacierwide'])

    workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs)
    workflow.calibrate_inversion_from_consensus(gdirs, apply_fs_on_mismatch=True)
    
else:    
    cfg.PARAMS['run_mb_calibration'] = True
    params_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'
    workflow.download_ref_tstars(base_url=params_url)

    # Now calibrate
    #cfg.PARAMS['mu_star_halfperiod'] = 3
    workflow.execute_entity_task(tasks.local_t_star, gdirs);
    workflow.execute_entity_task(tasks.mu_star_calibration, gdirs);

    t_mu_cru_ref = pd.DataFrame(index=range(len(gdirs)), columns=['ID', 'lon', 'lat','tstar', 'bias', 'mustar'])

    for i in range(len(gdirs)):
        print(gdirs[i].read_json('local_mustar')['t_star'])
        print(gdirs[i].read_json('local_mustar')['mu_star_per_flowline'])
        t_mu_cru_ref['ID'].iloc[i] = gdirs[i].rgi_id
        t_mu_cru_ref['lon'].iloc[i] = gdirs[i].cenlon
        t_mu_cru_ref['lat'].iloc[i] = gdirs[i].cenlat
        t_mu_cru_ref['tstar'].iloc[i] = gdirs[i].read_json('local_mustar')['t_star']
        t_mu_cru_ref['bias'].iloc[i] = gdirs[i].read_json('local_mustar')['bias']
        t_mu_cru_ref['mustar'].iloc[i] = gdirs[i].read_json('local_mustar')['mu_star_per_flowline']

    # Gepatsch
    t_mu_cru_ref.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU/t_mu_star_cru.csv')
    # Rofenache
    #t_mu_cru_ref.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU/t_mu_star_cru.csv')
    
#%%
# Compute the glacier mass balance:
years = np.arange(1969, 2020)
from oggm.core.massbalance import MultipleFlowlineMassBalance
for i in range(len(gdirs)):
    mbmod = MultipleFlowlineMassBalance(gdirs[i], use_inversion_flowlines=True)
    mb_ts = mbmod.get_specific_mb(year=years) # gets the annual mass balance
    plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)');
    mb = pd.DataFrame(index=years, data=mb_ts, columns={'mb':mb_ts})
    # Gepatsch
    mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU/MBcru_'+gdirs[i].rgi_id+'.csv')
    # Rofenache
    #mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU/MBcru_'+gdirs[i].rgi_id+'.csv')
    
#%%
# Convert the flowlines to a "glacier" for the ice dynamics module
init_year = 1969 #  gdirs[0].rgi_date # 
if init_glac == '':
    for gdir in gdirs:
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
        tasks.prepare_for_inversion(gdir), 
        tasks.mass_conservation_inversion(gdir), 
        tasks.filter_inversion_output(gdir), 
        tasks.init_present_time_glacier(gdir),
        tasks.run_from_climate_data(gdir,
                                    #ys=init_year+1,
                                    climate_filename='climate_historical')

elif init_glac == 'dyn_spinup':
    for gdir in gdirs:
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
        tasks.prepare_for_inversion(gdir),  
        tasks.mass_conservation_inversion(gdir),  
        tasks.filter_inversion_output(gdir), 
        tasks.init_present_time_glacier(gdir)
        tasks.run_dynamic_spinup(gdir,
                                 spinup_start_yr=init_year+1,
                                 minimise_for='area',
                                 precision_percent=17.0,
                                 output_filesuffix='',
                                 ye=gdir.get_climate_info()['baseline_hydro_yr_1']+1);
        

# workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs);

# workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, ys=1969,
#                               #fixed_geometry_spinup_yr=gdirs[0].rgi_date,
#                               climate_filename='climate_historical')

ds_run = utils.compile_run_output(gdirs)
ds_run

# Plot volume and length evolution of all the glaciers:
# time = years (in this example the total number of years was 200)
area = []
for i in range(len(gdirs)):
    area_i = gdirs[i].rgi_area_m2
    area.append(area_i)
  
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
ds_run.sum(dim='rgi_id').volume.plot.line(ax=ax1, hue='rgi_id');
ds_run.sum(dim='rgi_id').area.plot.line(ax=ax2, hue='rgi_id');
#ax2.set_ylim(2e7,7e7)
ax2.scatter(1969, 23092000, color='r'); # 61359000 for GEP, 23860000 all Gepatschalm
ax2.scatter(1998, 21070000, color='r'); # 52286000 for GEP, 21260000 all Gepatschalm
ax2.scatter(gdirs[0].rgi_date, sum(area), color='g')
ax2.scatter(2006, 20480000, color='r'); # 49650000 for GEP, 20600000 all Gepatschalm
ax2.scatter(2015, 18600000, color='r'); # 42566000 for GEP, 18620000 all Gepatschalm

#%%
# Get mass balance for every height on the tongue for Gepatschferner:
past_mb = massbalance.MultipleFlowlineMassBalance(gdirs[6])
heights = [2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525, 2575, 2625, 2675, 2725, 2775, 2825, 2875]
years_mb = list(np.arange(2012, 2020, 1))
mb_tongue = []
for y in years_mb:
    mb_tongue_aux = pd.DataFrame(past_mb.get_annual_mb(heights, y, fl_id=0))*cfg.SEC_IN_YEAR*cfg.PARAMS['ice_density']
    mb_tongue_aux.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU/MB_tongue/mb_'+str(y)+'.csv')
    mb_tongue.append(mb_tongue_aux)
