## Coupling scheme workflow
This text summarizes the main steps needed to follow for running WaSiM-OGGM
_____

1. Run WaSiM

Daily (or sub-daily) simulations. Monthly mean and sum grids for temperature and precipitation have to be created, respectively. This is done in the interpolation section of the control file.

Variable      | Write code  | Output                
:-----------: |:-----------:| :--------------------:
temperature   | 73          | mean monthly values   
precipitation | 33          | total monthly values 

_____
2. Convert grids to nc files

Use the file: [Convert WaSiM grids to netCDF files](../main/scripts/convert_raster_to_nc.py) \
Results: one netCDF file (i.e.: 'monthly_meteo.nc') containing monthly values of temperature and precipitation.

_____
3. Run OGGM

In case the glaciers do not belong to the reference glaciers, an adjustement of the meteorological input data for OGGM is required. This is needed because the parameters are obtained through interpolation from neighboring glaciers, which were already calibrated using default climate datasets.
The following steps are required:

   * Run OGGM with default climate dataset (CRU) with the example: [coupling_CRU](../main/scripts/coupling_CRU.py)
   * Run OGGM with user's dataset: [coupling_user](...) 
      * The user dataset (x_user) is adjusted (x_adj) to the default dataset (y_def) based on anomalies
      * In both cases, the adjustement is based on the mean monthly values of the variables.
      * Define the initialization method: 
         * No initialization: if OGGM starts its simulation from the RGI inventory's date (e.g. 2003)
         * Dynamic spinup: if OGGM starts its simulation no more than 40 years before the RGI inventory's date (e.g. in 1970)
         * Initialization_Eis: if OGGM starts its simulations far away in time, the method developed by Eis et al. (2021) will be applied. In this case, the script [coupling_initialization](../main/scripts/coupling_initialization.py) needs to be run beforehand

The glacier's outlines are created for each glacier in the directory and for each year defined in the simulation period, based on the flowline model. They are then converted into polygons and saved as shapefiles. Two new files are created for WaSiM:
   * glaciercells_year: an ASCII format file (raster) containing the fraction of glacierization for each grid cell and for each of the years
   * glaciercodes_year: an ASCII format file (raster) containing the glacier codes needed to run the glacier module in WaSiM (which corresponds each year with the glaciercells file).


_____
4. Re-run WaSiM with initial states from OGGM outputs


_____
5. Optimize WaSiM results (calibration)
