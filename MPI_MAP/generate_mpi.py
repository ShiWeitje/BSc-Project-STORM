#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 21:55:15 2022

@author: shiweiyuan
"""

import numpy as np
import matplotlib.pyplot as plt

#%% FUNCTIONS

def mpi_analysis(idx, m):
    
    mpi_maps = np.zeros(10)
    
    for i in range(10):
        y = 2008 + i
        mpi_map = np.loadtxt("MPI_MAP_BASIN"+str(idx)+"_MONTH"+str(m)+"_YEAR"+str(y)+".txt")
        mpi_maps[i] = mpi_map
    
    std_map = np.zeros(mpi_maps[0].shape)
    mean_map = np.zeros(mpi_maps[0].shape)
    
    for lat in range(mpi_maps[0].shape[0]):
        for lon in range(mpi_maps[0].shape[1]):
                mpis = mpi_maps[:,lat,lon]
                std = np.std(mpis)
                mean = np.mean(mpis)
                std_map[lat,lon] = std
                mean_map[lat,lon] = mean
    
    np.savetxt("MEAN_STD/MEAN_MAP_BASIN"+str(idx)+"_MONTH"+str(m)+".txt", mean_map)
    np.savetxt("MEAN_STD/STD_MAP_BASIN"+str(idx)+"_MONTH"+str(m)+".txt", std_map)
    
    mean_std = np.mean(std_map)
    mean_mpi = np.mean(mean_map)
    
    print("Average MPI across Grid: ", mean_mpi)
    print("Average STD across Grid: ", mean_std)
    return(mean_mpi, mean_std)

def mpi_generate(idx, m, no_years):
    
    MPI_means = np.loadtxt("MEAN_STD/MEAN_MAP_BASIN" + str(idx) +"_MONTH"+ str(m) + ".txt")
    MPI_stds = np.loadtxt("MEAN_STD/STD_MAP_BASIN" + str(idx) +"_MONTH"+ str(m) + ".txt")
    x, y = MPI_means.shape
    for year in range(no_years):
        MPI_gen = np.random.normal(MPI_means, MPI_stds)
        np.savetxt("GEN_MPI_FIELDS_BASIN" + str(idx) +"_MONTH"+ str(m) + "_YEAR"+ str(year) + ".txt", MPI_gen)
            


