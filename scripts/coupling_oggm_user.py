#!/usr/bin/env python
# coding: utf-8
'''
# # All glaciers in Gepatsch
# - run the simulations from 1969:
        a) no initialization or spin-up
        b) with dynamic spin-up
        c) with initialization method (Eis et al., 2021)
# - use the climate data from WaSiM (INCA-kNN) adjusted to CRU with anomalies
# - use the glacier's outline from OGGM (RGI 6.0)

# The calibration of the mass balance (mu_star) will be done
# from 2003 onwards, since that is the date for which we have the available glacier's
# outline and DEM. A spin-up from 1969 to 2003 is used.
'''
# In[1]:
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import oggm
import os
import sys
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import inversion
from oggm.core import massbalance, climate, flowline
import xarray as xr
import shapely.geometry as shpg
import gzip
import pickle
import glob
import shutil
import tarfile
import salem

funct_oggm_path = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Scripts'
sys.path.append(funct_oggm_path)
import funct_oggm

#%%
# Main definitions from the user:
# if the simulations are to be initialized with dynamic spinup or according to Eis et al. (2021)
init_glac = 'dyn_spinup' #'' #'initialization_Eis' # 
# Do we calibrate on geodetic mass balances?
cal_geod_mb = True
# Use adjusted user climate dataset to CRU?
adj_clim_data = False
# Get a monthly mass balance from hydro output after RGI's date? 
monthly_mb_output = True 

#%%
# Definition of parameters, paths and input data
cfg.initialize(logging_level='WARNING')
cfg.CONFIG_FILE
cfg.PARAMS['border']= 160
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['continue_on_error'] = True # Set to True for operational runs
#cfg.PARAMS['prcp_scaling_factor']= 1.0  # 1.0: no precipitation correction
cfg.PARAMS['baseline_climate'] = '' 

# Gepatsch
cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/ice_thickness' # Working directory: needs to be specified before each run
utils.mkdir(cfg.PATHS['working_dir'], reset=False)
cfg.PARAMS['climate_file'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/monthly_meteo_inca_1969-2019.nc' # adjusted through regression
cfg.PATHS['climate_file']  = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/monthly_meteo_inca_1969-2019.nc'
# Provide shapefile of the catchment
catchment = gpd.read_file('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/shp/catchment_alm.shp')

# # Rofenache
# cfg.PATHS['working_dir'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA'
# utils.mkdir(cfg.PATHS['working_dir'], reset=False)
# cfg.PARAMS['climate_file'] = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/monthly_meteo_inca_rofenache.nc' # adjusted through regression
# cfg.PATHS['climate_file']  = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/monthly_meteo_inca_rofenache.nc'
# #Provide shapefile of the catchment
# catchment = gpd.read_file('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/rofenache_total.shp')


#%%
# Definition of the glaciers that are going to be used for the runs
# in this case, all the valid glaciers in Gepatsch for the year 2003
fr = utils.get_rgi_region_file(11, version='62') # Get glaciers outlines for specified region
gdf = gpd.read_file(fr)
in_bas = [catchment.geometry.contains(shpg.Point(x,y))[0] for
          (x,y) in zip(gdf.CenLon, gdf.CenLat)]
gdf_sel = gdf.loc[in_bas]
#ax = catchment.plot();
#gdf_sel.plot(ax=ax, edgecolor='k')
#rgi_ids = ['RGI60-11.00897', 'RGI60-11.00787', 'RGI60-11.00719']

#%%
# # Create the glacier's directories
# gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=2, 
#                                             reset=False, force=True)
gdirs = workflow.init_glacier_directories(gdf_sel.RGIId, from_prepro_level=2, 
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

# Some plots:
# graphics.plot_centerlines(gdirs, figsize=(8, 7), use_flowlines=True, add_downstream=False)
# graphics.plot_catchment_areas(gdirs, figsize=(8,7)) 
# graphics.plot_catchment_width(gdirs, corrected=True, figsize=(8, 7))

#%%
# Climate tasks
# HERE IS climate_historical.nc CREATED!!!
workflow.execute_entity_task(tasks.process_custom_climate_data, gdirs);

if adj_clim_data:
    # Adjustment of climate data with anomalies
    # Gepatsch
    path_in_cru = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU' , 'per_glacier')
    path_in_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA' , 'per_glacier')
    temp_folder_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA-adj')
    # Rofenache
    # path_in_cru = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU' , 'per_glacier')
    # path_in_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA' , 'per_glacier')
    # temp_folder_user = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA-adj')

    rgi_ids = gdf_sel.RGIId.tolist() #rgi_ids 

    funct_oggm.clim_adj_anom(rgi_ids, path_in_cru, path_in_user, temp_folder_user)

    funct_oggm.replace_file(temp_folder_user, path_in_user, rgi_ids, 'climate_historical.nc')

####### HERE THE NEW CLIMATE_HISTORICAL.nc IS REPLACED AND USED ########
# In[34]:
# fpath = gdirs[0].get_filepath('climate_historical')
# ds = xr.open_dataset(fpath)
# # Data is in hydrological years
# # ds.prcp.resample(time='AS').sum()[1:-1].plot()
# ds.temp.resample(time='AS').mean()[1:-1].plot()
# print(ds)

# In[36]:
# ### Calibration of the mass balance volume ==> was done with the CRU run,
# so here we take the t_star from that simulation
if cal_geod_mb:
    workflow.execute_entity_task(tasks.mu_star_calibration_from_geodetic_mb, gdirs,
                              ignore_hydro_months=True, min_mu_star=80, max_mu_star=300)
    for i in range(len(gdirs)):
        print(gdirs[i].read_json('local_mustar')['t_star'])
        print(gdirs[i].read_json('local_mustar')['mu_star_glacierwide'])

    workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs)
    workflow.calibrate_inversion_from_consensus(gdirs, apply_fs_on_mismatch=True)
    
    
else:
    cfg.PARAMS['run_mb_calibration'] = False
    cfg.PARAMS['check_calib_params'] = False
    # #massbalance.PastMassBalance(gdirs[0])
    # Gepatsch
    df_tstar = pd.read_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/CRU/ref_tstars.csv',
                         index_col=0) #1, dtype={'tstar': int})
    # Rofenache
    #df_tstar = pd.read_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/CRU/ref_tstars.csv',
    #                    index_col=0, dtype={'tstar': int})

    #df_tstar.drop(columns=['Unnamed: 0'], inplace=True)
    for gdir in gdirs:
        climate.mb_yearly_climate_on_glacier(gdir)

    cfg.PARAMS['mu_star_halfperiod'] = 3
    # # for i in range(len(gdirs)):
    #     #     workflow.execute_entity_task(tasks.local_t_star, gdirs[i],
    #     #                                   tstar=df_tstar['tstar'].iloc[i], bias=df_tstar['bias'].iloc[i])
        
 
    for i in range(len(gdirs)):
        workflow.execute_entity_task(tasks.local_t_star, gdirs, ref_df=df_tstar);
        workflow.execute_entity_task(tasks.mu_star_calibration, gdirs);
        print(gdirs[i].read_json('local_mustar')['t_star'])
        print(gdirs[i].read_json('local_mustar')['mu_star_glacierwide'])

    # mu_star_ref = df_tstar['mustar']
    
    # # Copy t_star and mu_star from CRU run to current run
    # for glac in range(len(rgi_ids)):
    #     region = rgi_ids[glac].split('.')[0]
    #     src = glob.glob(path_in_cru+'/'+region+'/'+region+'.00'+'/'+rgi_ids[glac]+'/local_mustar.json')[0]
    #     dst = glob.glob(path_in_user+'/'+region+'/'+region+'.00'+'/'+rgi_ids[glac]+'/local_mustar.json')[0]
    #     shutil.copyfile(src, dst)
    #     print(gdirs[i].read_json('local_mustar')['t_star'])
    #     print(gdirs[i].read_json('local_mustar')['mu_star_glacierwide'])
    #     workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs) #compute inversion
    #     # using the mu_star obtained from CRU run (only if adj_clim_data is True)
        
#%%
# Compute the glacier mass balance:
years = np.arange(1969+1, 2020)
from oggm.core.massbalance import MultipleFlowlineMassBalance
for i in range(len(gdirs)):
    mbmod = MultipleFlowlineMassBalance(gdirs[i], use_inversion_flowlines=True)
    mb_ts = mbmod.get_specific_mb(year=years) # gets the annual mass balance
    plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)');
    mb = pd.DataFrame(index=years, data=mb_ts, columns={'mb':mb_ts})
    # Gepatsch
    # mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA/MB_'+gdirs[i].rgi_id+'.csv')
    # Rofenache
    #mb.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/INCA/MB_'+gdirs[i].rgi_id+'.csv')
    

#%%
# Run the simulations
# For past simulations (strting before 2003), specify if we want to 
# initialize the model with Eis et al. (2021) or do a dynamic spin-up
init_year = 1969 #  gdirs[0].rgi_date #  
if init_glac == '':
    for gdir in gdirs:
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
        tasks.prepare_for_inversion(gdir), 
        tasks.mass_conservation_inversion(gdir), 
        tasks.filter_inversion_output(gdir), 
        tasks.init_present_time_glacier(gdir),
        if not monthly_mb_output:
            tasks.run_from_climate_data(gdir,
                                  #ys=init_year+1,
                                  climate_filename='climate_historical')
        elif monthly_mb_output:
            tasks.run_with_hydro(gdir,
                             run_task=tasks.run_from_climate_data,
                             store_monthly_hydro=True)

elif init_glac == 'dyn_spinup':
    #cfg.PARAMS['check_calib_params'] = False
    for gdir in gdirs:
        mbmod = massbalance.MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
        tasks.prepare_for_inversion(gdir),  
        tasks.mass_conservation_inversion(gdir),  
        tasks.filter_inversion_output(gdir), 
        tasks.init_present_time_glacier(gdir)
        
        if not monthly_mb_output:
            tasks.run_dynamic_spinup(gdir,
                                      spinup_start_yr=init_year+1,
                                      minimise_for='area',
                                      precision_percent=17.0,
                                      output_filesuffix='',
                                      ye=gdir.get_climate_info()['baseline_hydro_yr_1']+1);

        elif monthly_mb_output:
            tasks.run_with_hydro(gdir,
                              run_task=tasks.run_dynamic_spinup,
                              store_monthly_hydro=True,
                              spinup_start_yr=init_year+1,
                              minimise_for='area',
                              precision_percent=17.0,
                              output_filesuffix='',
                              ye=gdir.get_climate_info()['baseline_hydro_yr_1']+1);


            # tasks.run_dynamic_spinup(gdir,
            #                           spinup_start_yr=init_year+1,
            #                           minimise_for='area',
            #                           precision_percent=17.0,
            #                           output_filesuffix='_dyn_spinup',
            #                           ye=gdir.rgi_date+1);
            
            # tasks.run_with_hydro(gdir,
            #                  run_task=tasks.run_from_climate_data,
            #                  store_monthly_hydro=True,
            #                  init_model_filesuffix= '_dyn_spinup',
            #                  output_filesuffix= '_with_hydro')
            
            # total_run = tasks.merge_consecutive_run_outputs(gdir,
            #                                                 input_filesuffix_1='_dyn_spinup',
            #                                                 input_filesuffix_2='_with_hydro',
            #                                                 output_filesuffix='')
                           
elif init_glac == 'initialization_Eis':

    path_init = os.path.join('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/Initialization' , 'per_glacier')

    # for gdir in gdirs:
    #     region = gdir.rgi_id.split('.')[0]
    #     fp = path_init+'/'+region+'/'+region+'.00'+'/'+gdir.rgi_id+'.tar.gz'
    #     # Define the path where the "tar" files of the initialization run are stored
    #     init_dir = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA/Initialization/'
    #     # Extract the "tar" files
    #     with tarfile.open(fp, 'r') as tf:
    #         fp_from_tar = tf.extractall(os.path.dirname(init_dir))
    #     # Read the results folder (pickle) for each of the glaciers. It starts with the init_year
    #     result_init = init_dir+gdir.rgi_id+'/result'+str(init_year)+'.pkl'
    #     # Unpickle the results folder
    #     with gzip.open(result_init, 'rb') as f:
    #         result = pickle.load(f)
    #     # Select the glacier geometry with the smallest fitness value, which will be the 
    #     # initial geometry for the initialization in 1969.
    #     best_fitness = result.sort_values('fitness', ascending=True)['fitness'].iloc[0]
    #     temp_bias_init = result.sort_values('fitness', ascending=True)['temp_bias'].iloc[0]
    #     time_init = result.sort_values('fitness', ascending=True)['time'].iloc[0]
    #     # File with the best fit, that will be used afterwards
    #     filename_best_fit = str(init_year)+'_past_*'+str(temp_bias_init)+'_'+str(time_init)+'.nc'
    #     file = glob.glob(os.path.join(init_dir, gdir.rgi_id, str(init_year), 'model_geometry'+filename_best_fit))[0]
    #     shutil.copy(file, gdir.dir)
    #     # First run: from initialization until ???

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
            md_geom = flowline.FileModel(mod_geom[j])
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
        
        
ds_run = utils.compile_run_output(gdirs)

# In[49]:
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
# Gepatsch
# ax2.scatter(1969, 23092000, color='r'); # 61359000 for GEP, 23860000 all Gepatschalm
# ax2.scatter(1998, 21070000, color='r'); # 52286000 for GEP, 21260000 all Gepatschalm
# ax2.scatter(gdirs[0].rgi_date, sum(area), color='g')
# ax2.scatter(2006, 20480000, color='r'); # 49650000 for GEP, 20600000 all Gepatschalm
# ax2.scatter(2015, 18600000, color='r'); # 42566000 for GEP, 18620000 all Gepatschalm
# Rofenache
ax2.scatter(1969, 23387000, color='r'); #  52237000 HEF, KAS, VER
ax2.scatter(1998, 21499000, color='r'); # 37740000 HEF, KAS, VER
ax2.scatter(gdirs[0].rgi_date, sum(area), color='g')
ax2.scatter(2006, 19631000, color='r'); # 32960000 HEF, KAS, VER
ax2.scatter(2015, 16905000, color='r'); # 27850000 HEF, KAS, VER

#%%
# Get mass balance for every height on the tongue for Gepatschferner:
past_mb = massbalance.MultipleFlowlineMassBalance(gdirs[6])
heights = [2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525, 2575, 2625, 2675, 2725, 2775, 2825, 2875]
years_mb = list(np.arange(2012, 2020, 1))
mb_tongue = []
for y in years_mb:
    mb_tongue_aux = pd.DataFrame(past_mb.get_annual_mb(heights, y, fl_id=0))*cfg.SEC_IN_YEAR*1000
    mb_tongue_aux.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/INCA/MB_tongue/mb_'+str(y)+'.csv')
    mb_tongue.append(mb_tongue_aux)

#%%
# Plot only for Gepatschferner
# f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
# ds_run.sum(dim='rgi_id').volume.plot.line(ax=ax1, hue='rgi_id');
# ds_run.sum(dim='rgi_id').area.plot.line(ax=ax2, hue='rgi_id');
# ax2.set_ylim(1.4e7,1.8e7)
# ax2.scatter(1969, 17951000, color='r'); # 61359000 for GEP, 23860000 all Gepatschalm
# ax2.scatter(1998, 17160000, color='r'); # 52286000 for GEP, 21260000 all Gepatschalm
# ax2.scatter(gdirs[0].rgi_date, gdirs[0].rgi_area_m2, color='g')
# ax2.scatter(2006, 16620000, color='r'); # 49650000 for GEP, 20600000 all Gepatschalm
# ax2.scatter(2015, 15770000, color='r'); # 42566000 for GEP, 18620000 all Gepatschalm

# In[50]:
# Plot all the outputs given by OGGM for Gepatschferner:
# years = np.arange(2003,2004)
# for i in years:
#     plt.rcParams["figure.figsize"] = [12,5]
#     graphics.plot_modeloutput_map(gdirs[6], modelyr=i, vmax=350, filesuffix='_dyn_spinup')
#     plt.tight_layout()
#     plt.savefig('/home/maria/Documents/DIRT-X_Gepatsch/Figures/gep_thickness_dynsp/'+str(i)+'.png', dpi=500)

# f, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
# graphics.plot_modeloutput_map(gdirs[6], modelyr=1998, ax=ax1, vmax=350, filesuffix='_dyn_spinup')
# graphics.plot_modeloutput_map(gdirs[6], modelyr=2006, ax=ax2, vmax=350, filesuffix='_hist_dyn_spinup')
# plt.tight_layout();
#plt.savefig('/home/maria/Documents/DIRT-X_Gepatsch/oggm/glacier_thickness.png', dpi=500)

# In[52]:

# From catchment widths to polygons to create the glacier outline
# For all flowlines together
years = np.arange(init_year+1, init_year+51)
crs_out = 'EPSG:31254'
# path_for_shp = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/shp/outlines_shp_oggm/'
path_for_shp = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/ice_thickness/shp/outlines_shp_oggm/'
 
for gdir in gdirs:
    for year in years:
        funct_oggm.from_width_to_polygon(gdir, year, crs_out, path_for_shp)

# In[53]
# Tasks for reading glacier's polygon shp, merging them into one and get the two input grids
# for WaSiM: glaciercells and glaciercodes.
path_oggm_outlines = path_for_shp
path_input_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/ice_thickness/shp/outlines_shp_for_wasim'
# dem_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/dem_100m.tif'
dem_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/shp/dem_wasim_100m.tif'
# grid_nodes_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/grid_nodes.shp'
grid_nodes_wasim = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/shp/grid.shp'
# subbasins = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/Subbasins.shp'
subbasins = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/shp/subbasins.shp'
resol_wasim = 100*100

for year in years:
    funct_oggm.prepare_grids_for_wasim(gdirs, year, crs_out, path_oggm_outlines, path_input_wasim, dem_wasim, grid_nodes_wasim, subbasins,
                                resol_wasim)

#%%
path_shp = path_for_shp
# path_shp = '/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Rofenache/ice_thickness/shp_all_glaciers/outlines_shp_oggm/'
for gdir in gdirs: #[19:]: #[:18]: #
    for year in years:
        # 1° create topo_smoothed and glacier_mask for each glacier
        topo_smoothed, glacier_mask, crs_mask = funct_oggm.topo_and_mask_per_glacier(gdir, year, path_shp)
        # 2° estimate the thickness distribution per glacier and year
        funct_oggm.glac_thickness(gdir, year, path_shp, glacier_mask)

#%%
# Get annual volumes and MB to be used in the coupling
output_vol = ds_run.sum(dim='rgi_id').volume
vol_yearly = pd.DataFrame(data=output_vol.data, index=output_vol.time, columns=['vol_m3'])
vol_yearly.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/ice_thickness/vol.csv')

output_area = ds_run.sum(dim='rgi_id').area
mb_yearly = pd.DataFrame(index=output_vol.time, columns=['vol_m3', 'area_m2', 'mb_mm'])
mb_yearly['vol_m3'] = vol_yearly['vol_m3']
mb_yearly['area_m2'] = output_area
mb_yearly['mb_mm'] = (vol_yearly['vol_m3'].diff()/mb_yearly['area_m2']) * 918

mb_yearly.to_csv('/home/maria/Documents/DIRT-X_Gepatsch/coupling/oggm/Total_area/ice_thickness/mb_oggm.csv')
