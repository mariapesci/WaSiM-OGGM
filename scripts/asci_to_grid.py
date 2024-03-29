# -*- coding: utf-8 -*-
"""
@author: María Pesci

This script takes the annual grid files containing the ice thicknesses (and areas) and glacier codes
created from OGGM's outputs and convert them into ASCII format.

Needed:
    - path: path with the glacier files for WaSIM:
            a) 'glaciercells*.asc' (where * refers to the year), contains the glacier thickness at each
                grid cell
            b) 'glaciercodes*.asc' contains the glacier codes at each grid cell
            c) the ascigrid.exe executable needs to be located in the same folder
"""

import glob
import os
import subprocess

path = r'E:\DIRT_X\Gepatsch\Coupling\wasim\init_states_glac_vol_oggm1.6\dyn_spinup'

glaccells = glob.glob(path+'\glaciercells*.asc')
glaccodes = glob.glob(path+'\glaciercodes*.asc')

for glac in glaccells:
    asc_name = glac.split(path)[1].split('\\')[1]
    grd_name = 'glcalm'+asc_name.split('glaciercells')[1].split('.asc')[0]+'.grd'
    cmd = ['%s' % path + '/ascigrid.exe', asc_name, grd_name]
    os.chdir(path)
    p = subprocess.Popen(cmd)

for glac in glaccodes:
    asc_name = glac.split(path)[1].split('\\')[1]
    grd_name = 'glidalm'+asc_name.split('glaciercodes')[1].split('.asc')[0]+'.grd'
    cmd = ['%s' % path + '/ascigrid.exe', asc_name, grd_name]
    os.chdir(path)
    p = subprocess.Popen(cmd)
