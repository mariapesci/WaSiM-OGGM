#!/usr/bin/env python
# coding: utf-8
'''
# # All glaciers in Gepatschalm
# - run the simulations from 1969:
        a) no initialization or spin-up
        b) with dynamic spin-up
        c) with initialization method (Eis et al., 2021)
# - use the climate data from WaSiM (INCA-kNN) and adjusted to CRU with anomalies if needed
# - use the glacier's outline from OGGM (RGI 6.2)

# The calibration of the mass balance (mu_star) will be done
# from 2003 onwards, since that is the date for which we have the available glacier's
# outline and DEM. A spin-up from 1969 to 2003 is used.
'''
# In[1]:
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import oggm
import os
import sys
from oggm import cfg, utils, workflow, tasks#, graphics
from oggm.core import massbalance, inversion
from oggm.core.massbalance import mb_calibration_from_scalar_mb, mb_calibration_from_geodetic_mb, mb_calibration_from_wgms_mb
import xarray as xr
import shapely.geometry as shpg
import gzip
import pickle
import glob
import shutil
import tarfile
import salem
import json

funct_oggm_path = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Scripts'
sys.path.append(funct_oggm_path)
import funct_oggm

#%%
# Main definitions from the user:
# if the simulations are to be initialized with dynamic spinup or according to Eis et al. (2021)
init_glac = 'dyn_spinup' #'' #'initialization_Eis' # 
# Do we calibrate on geodetic mass balances?
cal_geod_mb = False
# Use adjusted user climate dataset to W5E5?
adj_clim_data = False 
# Get a monthly mass balance from hydro output after RGI's date? 
monthly_mb_output = False
# Perform a sensitivity test or just calibrate mb?
mb_calib_option = 'calibrate' #'sensitivity' #'use_defaults' #

#%%
# Definition of parameters, paths and input data
cfg.initialize(logging_level='WARNING')
cfg.CONFIG_FILE
cfg.PARAMS['border']= 160
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['continue_on_error'] = True # Set to True for operational runs

cfg.PARAMS['prcp_fac'] = 1.0
cfg.PARAMS['baseline_climate'] = '' 
cfg.PARAMS['use_temp_bias_from_file'] = False
cfg.PARAMS['use_winter_prcp_fac'] = False
cfg.PARAMS['evolution_model'] = 'FluxBased'
cfg.PARAMS['downstream_line_shape']='parabola' # to make sure we avoid having bed shapes after the terminus

# Gepatschalm
cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/future_projections/dyn_spinup_calib/'#'Run6/' # Working directory: needs to be specified before each run
utils.mkdir(cfg.PATHS['working_dir'], reset=False)
cfg.PARAMS['climate_file'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/monthly_meteo_coupling.nc' #monthly_meteo_inca_1969-2019.nc' #''monthly_meteo_coupling_bs_1_45.nc' #'' # adjusted through regression
cfg.PATHS['climate_file']  = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/monthly_meteo_coupling.nc' #monthly_meteo_inca_1969-2019.nc' #monthly_meteo_coupling_bs_1_45.nc' #
# Provide shapefile of the catchment
catchment = gpd.read_file('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/catchment_alm.shp')

#%%
# Definition of the glaciers that are going to be used for the runs
# in this case, all the valid glaciers in Gepatsch for the year 2003
fr = utils.get_rgi_region_file(11, version='62') # Get glaciers outlines for specified region
gdf = gpd.read_file(fr)
in_bas = [catchment.geometry.contains(shpg.Point(x,y))[0] for
          (x,y) in zip(gdf.CenLon, gdf.CenLat)]
gdf_sel = gdf.loc[in_bas]
ax = catchment.plot();
gdf_sel.plot(ax=ax, edgecolor='k')
#rgi_ids = ['RGI60-11.00746'] #['RGI60-11.00897', 'RGI60-11.00787', 'RGI60-11.00719']

#%%
# # Create the glacier's directories
# gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=2, 
#                                             reset=False, force=True)
base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L1-L2_files/centerlines/'
gdirs = workflow.init_glacier_directories(gdf_sel.RGIId, prepro_base_url=base_url, from_prepro_level=2, 
                                          reset=False, force=True)
print(gdirs[0].rgi_date)

# Is possible to create a list of tasks to be appplied to the list of glaciers
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

from oggm import graphics
# # Some plots:
graphics.plot_domain(gdirs[6])
graphics.plot_centerlines(gdirs[6], figsize=(8, 7), use_flowlines=True, add_downstream=False, add_line_index=False)
graphics.plot_catchment_areas(gdirs[6], figsize=(8,7)) 
graphics.plot_catchment_width(gdirs[6], corrected=True, figsize=(8, 7))

#%%
# Climate tasks
# HERE IS climate_historical.nc CREATED!!!
workflow.execute_entity_task(tasks.process_custom_climate_data, gdirs);

if adj_clim_data:
    # Adjustment of climate data with anomalies
    # Gepatschalm
    path_in_def = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run1' , 'per_glacier')
    path_in_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6-adj' , 'per_glacier')
    temp_folder_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/INCA-adj')
    # Rofenache
    # path_in_cru = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU' , 'per_glacier')
    # path_in_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA' , 'per_glacier')
    # temp_folder_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA-adj')

    rgi_ids = gdf_sel.RGIId.tolist() #rgi_ids 

    funct_oggm.clim_adj_anom(rgi_ids, path_in_def, path_in_user, temp_folder_user)

    funct_oggm.replace_file(temp_folder_user, path_in_user, rgi_ids, 'climate_historical.nc')

    funct_oggm.copy_file(path_in_def, path_in_user, rgi_ids, 'mb_calib.json')
####### HERE THE NEW CLIMATE_HISTORICAL.nc IS REPLACED AND USED ########

#%%
fpath = gdirs[0].get_filepath('climate_historical')
ds = xr.open_dataset(fpath)
# Data is in hydrological years
#ds.prcp.resample(time='AS').sum()[1:-1].plot()
ds.temp.resample(time='AS').mean()[1:-1].plot()
#print(ds)

# In[36]:
if mb_calib_option == 'sensitivity':
    # Run sensitivity analysis of calibration parameters of the mass balance
    path_calib = cfg.PATHS['working_dir']+'calibration_sensitivity_init2003/'

    for gdir in gdirs:
        funct_oggm.sens_param_mb(gdir, path_calib, ref_period=cfg.PARAMS['geodetic_mb_period'],
                     default_param='melt_f', test_sens_param='prcp_fac', calib_param='temp_bias',
                     save_fig=True)         

elif mb_calib_option == 'calibrate':
    for gdir in gdirs:
        funct_oggm.calib_one_param_mb(gdir, default_param1='melt_f', default_param1_val=5, 
                               default_param2='temp_bias', default_param2_val=0.5,
                               calib_param='prcp_fac', ref_period=cfg.PARAMS['geodetic_mb_period'])

elif mb_calib_option == 'use_defaults':
    for gdir in gdirs:
        mbmod = massbalance.MonthlyTIModel(gdir, check_calib_params=False)

    
workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs)#, mb_model=mbmod)
# Calibrate inversion from consensus: finds the "best Glen A" to match all glaciers with a valid inverted vol.
workflow.calibrate_inversion_from_consensus(gdirs, apply_fs_on_mismatch=True)

#%%
# # Compute the glacier mass balance:
years = np.arange(1969, 2020)
from oggm.core.massbalance import MultipleFlowlineMassBalance
for gdir in gdirs:
    # mbmod already includes MonthlyTIModel which considers the calibrated parameters
    mbmod = MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
    mb_ts = mbmod.get_specific_mb(year=years) # gets the annual mass balance
    plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)');
    mb = pd.DataFrame(index=years, data=mb_ts, columns={'mb':mb_ts})
    # Gepatsch
    mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/MB_'+gdir.rgi_id+'.csv')
    # Rofenache
    #mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA/MB_'+gdirs[i].rgi_id+'.csv')

#%%
# Run the simulations
# For past simulations (strting before 2003), specify if we want to 
# initialize the model with Eis et al. (2021) or do a dynamic spin-up
obs_areas = pd.read_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/obs_areas_2006_2015.csv',
                        index_col='RGIId')
for gdir in gdirs:
    if mb_calib_option == 'sensitivity':
        tb, pf, _ = funct_oggm.get_best_calib_param(gdir, path_calib, 'sens_glac_tb_cal_pf_', obs_areas, melt_f_def_value=5, 
                                                prcp_fac_def_value= cfg.PARAMS['prcp_fac'],
                                                default_param='melt_f', test_sens_param='temp_bias', 
                                                calib_param='prcp_f')
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True,
                                                    temp_bias=tb, prcp_fac=pf)
    elif mb_calib_option == 'calibrate' or mb_calib_option == 'use_defaults':
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True,
                                                        check_calib_params=False)
        tasks.prepare_for_inversion(gdir), 
        tasks.mass_conservation_inversion(gdir), 
        tasks.filter_inversion_output(gdir), 
        tasks.init_present_time_glacier(gdir),
        
        if init_glac == '':
           init_year = gdirs[0].rgi_date
           if not monthly_mb_output:
               tasks.run_from_climate_data(gdir,
                                        mb_model=mbmod,
                                  #ys=init_year+1,
                                  climate_filename='climate_historical',
                                  #output_filesuffix='_historical',
                                  store_fl_diagnostics=True)
           elif monthly_mb_output:
               tasks.run_with_hydro(gdir,
                             run_task=tasks.run_from_climate_data,
                             store_monthly_hydro=True,
                             store_fl_diagnostics=True)

        elif init_glac == 'dyn_spinup':
            init_year = 1969
            if not monthly_mb_output:
                tasks.run_dynamic_spinup(gdir,
                                      spinup_start_yr=init_year+1,
                                      minimise_for='area',
                                      precision_percent=17.0,
                                      output_filesuffix='_historical',
                                      ye=2010, #gdir.get_climate_info()['baseline_yr_1']+1,
                                      store_fl_diagnostics=True);

            elif monthly_mb_output:
                tasks.run_with_hydro(gdir,
                              run_task=tasks.run_dynamic_spinup,
                              store_monthly_hydro=True,
                              spinup_start_yr=init_year+1,
                              minimise_for='area',
                              precision_percent=17.0,
                              output_filesuffix='',
                              ye=gdir.get_climate_info()['baseline_yr_1']+1,
                              store_fl_diagnostics=True);


    if init_glac == 'initialization_Eis':
        init_year = 1969
        path_init = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/Initialization' , 'per_glacier')

   #    Selecting best fit between simulated and observed areas at init_year
        for gdir in gdirs:
            init_dir = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA/Initialization/'
            glac_areas_1969 = pd.read_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/RGI_areas_1969.csv',
                                  index_col='RGIId')
            filename_mod_geom = str(init_year)+'_past_*'+'.nc'
            mod_geom = glob.glob(os.path.join(init_dir, gdir.rgi_id, str(init_year), 'model_geometry'+filename_mod_geom))
            area_candidates = pd.DataFrame(index=range(len(mod_geom)), columns=['candidate_name', 'area_1969', 'diff'])
            for i in glac_areas_1969.index:
                if i == gdir.rgi_id:
                    area_obs = glac_areas_1969['area_1969'].loc[i]
            for j in area_candidates.index:
                md_geom = oggm.flowline.FileModel(mod_geom[j])
                area_candidates.loc[j] = [mod_geom[j].split('1969/')[1], md_geom.area_km2, abs(md_geom.area_km2-area_obs)]
 
            best_candidate = area_candidates.sort_values('diff', ascending=True)['candidate_name'].iloc[0]
            file = glob.glob(os.path.join(init_dir, gdir.rgi_id, str(init_year), best_candidate))[0]
            shutil.copy(file, gdir.dir)
    
            mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
            tasks.prepare_for_inversion(gdir),  # This is a preprocessing task
            tasks.mass_conservation_inversion(gdir),  # This does the actual job (glacier thickness along flowlines)
            tasks.filter_inversion_output(gdir), 
            tasks.init_present_time_glacier(gdir)
        
            tasks.run_from_climate_data(gdir,
                                    ys=init_year+1,
                                    #ye=gdir.rgi_date,
                                    init_model_filesuffix=file.split('geometry')[1].split('.nc')[0],
                                    init_model_yr=init_year+1,
                                    output_filesuffix='')
        
        
ds_run = utils.compile_run_output(gdirs, input_filesuffix='_historical')

# ##### Plot volume and length evolution of all the glaciers:
# time = years (in this example the total number of years was 200)
area = []
for i in range(len(gdirs)):
    area_i = gdirs[i].rgi_area_m2
    area.append(area_i)
  
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
ds_run.sum(dim='rgi_id').volume.plot.line(ax=ax1, hue='rgi_id');
ds_run.sum(dim='rgi_id').area.plot.line(ax=ax2, hue='rgi_id');
#ax2.set_ylim(1.4e7,2.4e7)
# Gepatschalm
ax2.scatter(1969, 23092000, color='r'); # 61359000 for GEP, 23860000 all Gepatschalm
ax2.scatter(1998, 21070000, color='r'); # 52286000 for GEP, 21260000 all Gepatschalm
ax2.scatter(gdirs[0].rgi_date, sum(area), color='g')
ax2.scatter(2006, 20480000, color='r'); # 49650000 for GEP, 20600000 all Gepatschalm
ax2.scatter(2015, 18600000, color='r'); # 42566000 for GEP, 18620000 all Gepatschalm


# Get mass balance for every height on the tongue for Gepatschferner:
past_mb = massbalance.MultipleFlowlineMassBalance(gdirs[6])
heights = [2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525, 2575, 2625, 2675, 2725, 2775, 2825, 2875]
years_mb = list(np.arange(2012, 2020, 1))
mb_tongue = []
for y in years_mb:
    mb_tongue_aux = pd.DataFrame(past_mb.get_annual_mb(heights, y, fl_id=0))*cfg.SEC_IN_YEAR*1000
    mb_tongue_aux.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/dyn_spinup_calibrated/MB_tongue/mb_'+str(y)+'.csv')
    mb_tongue.append(mb_tongue_aux)


#%%
# From catchment widths to polygons to create the glacier outline
# For all flowlines together
years = np.arange(init_year+1, init_year+51) #+17
crs_out = 'EPSG:31254'
path_for_shp = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/dyn_spinup_calibrated/shp/outlines_shp_oggm/'
 
for gdir in gdirs:
    for year in years:
        funct_oggm.from_width_to_polygon(gdir, year, crs_out, path_for_shp, save_shp=True)

#%%
# Tasks for reading glacier's polygon shp, merging them into one and get the two input grids
# for WaSiM: glaciercells and glaciercodes.
path_oggm_outlines = path_for_shp
path_input_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/dyn_spinup_calibrated/shp/outlines_shp_for_wasim'
dem_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/dem_wasim_100m.tif'
grid_nodes_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/grid.shp'
subbasins = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/subbasins.shp'
resol_wasim = 100*100

for year in years:
    funct_oggm.prepare_grids_for_wasim(gdirs, year, crs_out, path_oggm_outlines, path_input_wasim, dem_wasim, grid_nodes_wasim, subbasins,
                                resol_wasim)

#%%
path_shp = path_for_shp

for gdir in gdirs: #[:10]: #[19:]: #[:18]: # (glacier 00783 doesn't exist after 2012 so it won't work)
    for year in years:
        # 1° create topo_smoothed and glacier_mask for each glacier
        topo_smoothed, glacier_mask, crs_mask = funct_oggm.topo_and_mask_per_glacier(gdir, year, path_shp)
        # 2° estimate the thickness distribution per glacier and year
        funct_oggm.glac_thickness(gdir, year, path_shp, glacier_mask)

#%%
# Get annual volumes and MB to be used in the coupling
output_vol = ds_run.sum(dim='rgi_id').volume
vol_yearly = pd.DataFrame(data=output_vol.data, index=output_vol.time, columns=['vol_m3'])
vol_yearly.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/dyn_spinup_calibrated/vol.csv')

output_area = ds_run.sum(dim='rgi_id').area
mb_yearly = pd.DataFrame(index=output_vol.time, columns=['vol_m3', 'area_m2', 'mb_mm'])
mb_yearly['vol_m3'] = vol_yearly['vol_m3']
mb_yearly['area_m2'] = output_area
mb_yearly['mb_mm'] = (vol_yearly['vol_m3'].diff()/mb_yearly['area_m2']) * 918

mb_yearly.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm_v16/Gepatschalm/Run6/dyn_spinup_calibrated/mb_oggm.csv')

