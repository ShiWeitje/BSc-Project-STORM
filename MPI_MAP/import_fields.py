# -*- coding: utf-8 -*-
"""
This module is part of the STORM model

For more information, please see 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

Functions described here are part of the data pre-processing and calculate the environmental
conditions + wind-pressure relationship.

Copyright (C) 2020 Nadia Bloemendaal. All versions released under the GNU General Public License v3.0
"""

import numpy as np
import os
import sys
import xarray as xr
dir_path=os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

####### GET FIELDS PER YEAR FOR 2008-2017 ###############################

def pressure_peryear(data):
    
    mslp=data.msl.values #mslp=data.msl.values, changed by Shirin, 29/09/2021
    # lon=data.longitude.values
    # lat=data.latitude.values
    
    for month in range(0,12):
        for t in range(28,38):    
            year = 1980 + t
            mslp_field = mslp[month+t*12,:,:]
            np.savetxt('MSLP_FIELDS/MSLP_month'+str(month+1)+'_year'+str(year)+'.txt',mslp_field)

def sst_peryear(data):
    
    sst=data.sst.values #mslp=data.msl.values, changed by Shirin, 29/09/2021
    # lon=data.longitude.values
    # lat=data.latitude.values
    
    for month in range(0,12):
        for t in range(28,38):    
            year = 1980 + t
            sst_field = sst[month+t*12,:,:]
            np.savetxt('SST_FIELDS/SST_month'+str(month+1)+'_year'+str(year)+'.txt',sst_field)

def t_peryear(data):
    
    t = data.t.values
    print("t loaded")
    # lon = data.longitude.values
    # print("lon loaded")
    # lat = data.latitude.values
    # print("lat loaded")
    levels = data.level.values
    print("levels loaded")
    
    for month in range(12):
        for y in range(0, 10):
            year = y + 2008
            t_field = t[month+y*12, ::-1, :, :]
            for level in range(len(levels)):
                t_field_level = t_field[level, :, :]
                np.savetxt("T_FIELDS/T_month"+str(month)+"_year"+str(year)+"_level"+str(level)+".txt", t_field_level)
                print("T Month %i, year %i, level %i saved" % (month, year, level))
    
def q_peryear(data):
    
    q = data.q.values
    print("q loaded")
    # lon = data.longitude.values
    # print("lon loaded")
    # lat = data.latitude.values
    # print("lat loaded")
    levels = data.level.values
    print("levels loaded")
    
    for month in range(12):
        for y in range(0, 10):
            year = y + 2008
            q_field = q[month+y*12, ::-1, :, :]
            for level in range(len(levels)):
                q_field_level = q_field[level, :, :]
                np.savetxt("Q_FIELDS/Q_month"+str(month)+"_year"+str(year)+"_level"+str(level)+".txt", q_field_level)
                print("Q Month %i, year %i, level %i saved" % (month, year, level))


#%%

data=xr.open_dataset('Monthly_mean_SST.nc')
sst_peryear(data)
data.close()

print("SST done")

data=xr.open_dataset('Monthly_mean_MSLP.nc')
pressure_peryear(data)
data.close()

print("MSLP done")

data=xr.open_dataset("Monthly_Mean_T.nc")
t_peryear(data)
data.close()

print("T done")


data=xr.open_dataset("Monthly_Mean_SH.nc")
q_peryear(data)
data.close()

print("Q done")
