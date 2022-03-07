#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 17:15:11 2022

@author: shiweiyuan
"""

from pcmin_fullterm import pcmin
import numpy as np
import preprocessing
from tqdm import tqdm
import matplotlib.pyplot as plt

#%%


def MPI_gridpoint(sst, mslp, p_grad, t_grad, q_grad):
    # formats:
    # sst: value
    # mslp: value
    # t_grad: array
    # p_grad: array
    # q_grad: array
    
    tmf = pcmin(sst, mslp, p_grad, t_grad, q_grad)
    
    PMIN = tmf[0]
    #VMAX = tmf[1]
    #TO = tmf[2]
    #IFL = tmf[3]
    #RAT = tmf[4]
    #CAPEMS = tmf[5]
    #CAPEM = tmf[6]
    #FAC = tmf[7]
    #CAPEA = tmf[8]
    return(PMIN)

def MPI_map(idx, months_override = []):
    
    p_grad = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800,
                       775, 750, 700, 650, 600, 550, 500, 450, 400,
                       350, 300, 250, 225, 200, 175, 150, 125, 100,
                       70, 50, 30, 20, 10, 7, 5, 3, 2, 1])
    
    months = np.array([[6,7,8,9,10,11],[6,7,8,9,10,11],[4,5,6,9,10,11],[1,2,3,4,11,12],[1,2,3,4,11,12],[5,6,7,8,9,10,11]])
    
    lat0,lat1,lon0,lon1 = preprocessing.BOUNDARIES_BASINS(idx)
    lat0 = (lat0 + 90) * 4
    lat1 = (lat1 + 90) * 4
    lon0 *= 4
    lon1 *= 4
    
    mpi_months = months[idx]
    
    if len(months_override) > 0:
        mpi_months = months_override
    
    for m in mpi_months:
        for y in range(2008, 2018):
            sst_globe = np.flip(np.loadtxt('SST_FIELDS/SST_month'+str(m)+'_year'+str(y)+'.txt'), 0)
            sst_basin = sst_globe[lat0:lat1, lon0:lon1]
            # plt.title("SST MAP")
            # plt.pcolormesh(sst_basin)
            # plt.show()
            
            np.savetxt("sst_basin_test.txt", sst_basin)
            
            mslp_globe = np.flip(np.loadtxt('MSLP_FIELDS/MSLP_month'+str(m)+'_year'+str(y)+'.txt'), 0)
            mslp_basin = mslp_globe[lat0:lat1, lon0:lon1]
            # plt.title("MSLP MAP")
            # plt.pcolormesh(mslp_basin)
            # plt.show()
            
            t_grads = np.zeros((len(p_grad), lat1-lat0, lon1-lon0))
            q_grads = np.zeros((len(p_grad), lat1-lat0, lon1-lon0))
            
            print("Loading in T values")
            for level in tqdm(range(len(p_grad))):
                tgrad_map_level = np.flip(np.loadtxt("T_FIELDS/T_month"+str(m)+"_year"+str(y)+"_level"+str(level)+".txt"),1)[lat0:lat1, lon0:lon1]
                t_grads[level,:,:] = tgrad_map_level
                #print("t_grad level %i loaded" % level)
            
            #plt.title("T Level 0")
            #plt.pcolormesh(t_grads[0])
            #plt.show()
            
            print("Loading in Q values")
            for level in tqdm(range(len(p_grad))):
                qgrad_map_level = np.flip(np.loadtxt("Q_FIELDS/Q_month"+str(m)+"_year"+str(y)+"_level"+str(level)+".txt"),1)[lat0:lat1, lon0:lon1]
                q_grads[level,:,:] = qgrad_map_level
                #print("q_grad level %i loaded" % level)
            
            #plt.title("SH Level 0")
            #plt.pcolormesh(q_grads[0])
            #plt.show()
            
            ########### SET TO CORRECT UNITS
            
            mslp_basin = mslp_basin / 1e2
            sst_basin = sst_basin - 273.15
            t_grads = t_grads - 273.15
            q_grads = q_grads * 1e3
            
            ###########
            
            mpi_map_month = np.zeros(sst_basin.shape)
            
            print("Calculating MPI Field, Month %i, Year %i" % (m,y))
            for lat in tqdm(range(sst_basin.shape[0])):
                for lon in range(sst_basin.shape[1]):
                    
                    t_grad = t_grads[:,lat,lon]
                    q_grad = q_grads[:,lat,lon]
                    sst = sst_basin[lat,lon]
                    mslp = mslp_basin[lat,lon]
                    """
                    if sst<-1000.: # LAND
                        mpi_map_month[lat,lon] = np.nan
                    else: # OCEAN
                        tmp = pcmin(sst, mslp, p_grad, t_grad, q_grad)
                        if ~np.isnan(tmp[0]) and tmp[0]!=0:
                            mpi_map_month[lat,lon] = tmp[0]
                        else:
                            mpi_map_month[lat,lon] = np.nan
                    """
                    if np.isnan(sst) == True:
                        #print("SST is nan")
                        mpi_map_month[lat,lon] = np.nan
                        continue
                    # COMMENT PRINTS OUT FOR SPEED
                    #print("tgrad: ", t_grad)
                    #print("qgrad: ", q_grad)
                    #print("sst: ", sst)
                    #print("mslp: ", mslp)
                    
                    mpi = MPI_gridpoint(sst, mslp, p_grad, t_grad, q_grad)
                    #print("mpi: ", mpi)
                    mpi_map_month[lat, lon] = mpi
                    
                    #print("MPI Calc. lat %i lon %i calculated" % (lat, lon))
            
            np.savetxt("MPI_MAPS/MPI_MAP_BASIN" +str(idx)+ "_MONTH" +str(m)+"_YEAR"+str(y)+".txt", mpi_map_month)
            print("Month %i, Year %i MPI MAP done." % (m,y))
        
#%% CREATE MPI MAPS FOR YEARS 2008-2017 FOR SPECIFIC BASIN

# Basin indices: 
# 0 = EP = Eastern Pacific
# 1 = NA = North Atlantic
# 2 = NI = North Indian
# 3 = SI = South Indian
# 4 = SP = South Pacific
# 5 = WP = Western Pacific

idx = 5

MPI_map(idx)

#%% PLOT WP MPI FIELDS FOR EACH MONTH,YEAR

idx = 5

for m in [5,6,7,8,9,10,11]:
    for y in range(2008,2018):
        mpi_map = np.loadtxt("MPI_MAP_BASIN"+str(idx)+"_MONTH" +str(m)+"_YEAR"+str(y)+".txt")
        plt.pcolormesh(mpi_map)
        plt.title("WP MPI Map, Month %i Year %i" % (m,y))
        plt.colorbar()
        plt.show()

    
