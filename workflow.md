## Coupling scheme workflow
This text summarizes the main steps needed to follow for running WaSiM-OGGM
_____

1. Run WaSiM

Daily (or sub-daily) simulations. Monthly mean and sum grids for temperature and precipitation have to be created, respectively. This is done in the interpolation section of the control file.

Variable      | Write code  | Output                
:-----------: |:-----------:| :--------------------:
temperature   | 73          | mean monthly values   
precipitation | 33          | total monthly values 

Since we are going to update the glacier model with OGGM's output afterwards, the glacier model in WaSiM can be deactivated.
_____
2. Convert grids to nc files

Use the file: [Convert WaSiM grids to netCDF files](../main/scripts/convert_raster_to_nc.py) \
Results: one netCDF file (i.e.: 'monthly_meteo.nc') containing monthly values of temperature and precipitation.

_____
3. Run OGGM

OGGM calculates the monthly mass balance based on the Temperature index model ([Marzeion et al., 2012](https://tc.copernicus.org/articles/6/1295/2012/tc-6-1295-2012.html)). In this coupling scheme, we try to use the default parameters from OGGM (and avoid extra calibration steps). Since many of the glaciers are not reference glaciers (without observed mass balances), there are two options on how to run the model:
  - First option: with user's climate dataset (e.g. from WaSiM) calculate glacier mass balances and calibrate them with geodetic mass balance data [T-index model calibrated on geodetic MB data](https://docs.oggm.org/en/stable/mass-balance-2012-pergla.html).
  - Second option: adjust the climate input data for OGGM. This is needed because the parameters are obtained through interpolation from neighboring glaciers, which were already calibrated using default climate datasets. In case this option is selected, a "first" OGGM with the default climate dataset CRU is required: [coupling_CRU](../main/scripts/coupling_CRU.py)
  - Run OGGM with user's dataset: [coupling_user](...) 
    * If the second option is selected, the user dataset (x_user) is adjusted (x_adj) to the default dataset (y_def) based on anomalies (the adjustement is based on the mean monthly values of the variables)
    * Define the initialization method: 
      * No initialization: if OGGM starts its simulation from the RGI inventory's date (e.g. 2003)
      * Dynamic spinup: if OGGM starts its simulation no more than 40 years before the RGI inventory's date (e.g. in 1970)
      * Initialization_Eis: if OGGM starts its simulations far away in time, the method developed by Eis et al. (2021) will be applied. In this case, the script [coupling_initialization](../main/scripts/coupling_initialization.py) needs to be run beforehand

Based on the [ice dynamics flowline model](https://docs.oggm.org/en/stable/ice-dynamics.html), an outline is created for each of the glaciers and each year within the simulation period. The outlines are converted into polygons and saved as shapefiles for posterior use.\
Similarly, the glacier's thickness distribution within the previously defined outilne, is obtained for each of the glaciers and years, adapted from the function [distribute_thickness_per_altitude](https://github.com/OGGM/oggm/blob/9d173038862f36a21838034da07243bd189ef2d0/oggm/core/inversion.py#L747).\
\
Two new files are created for WaSiM:
   * glaciercells_year: an ASCII format file (raster) containing the volume of ice for each grid cell and for each of the years
   * glaciercodes_year: an ASCII format file (raster) containing the glacier codes needed to run the glacier module in WaSiM (which corresponds each year with the glaciercells file).

_____
4. Re-run WaSiM with initial states from OGGM outputs

In WaSiM, the dynamic glacier model is selected. This model is based on the volume-area scaling approach $V = b.A^f$ and depends on two empiric factors (b = mean glacier thickness of a 1 km2 glacier and f = scaling factor). \
Now, a new simulation with WaSiM can be performed. Each year, the integrated VA-scaling model in WaSiM is then "replaced" by OGGM's outputs (area x ice thickness = volume). \
WaSiM is run annually: at the beginning of each year, the glaciers' volume is known (OGGM outputs). At the end of each year, the outputs from WaSiM serve as initial states for the following year, thus ensuring a continuous simulation through the entire simulation period.\
WaSiM determines the mass balance of the glaciers on a daily basis and based on their own parameter set. The results differ from OGGM mass balances, since other parameters set is used (inherent to OGGM). Thus, a calibration based on the annual (or monthly) mass balances must be performed, in order to integrate OGGM's results. \
\

_____
5. Optimize WaSiM results (calibration)
