README FOR MPI_MAP:

Code adapted from STORM model, see https://doi.org/10.1038/s41597-020-0381-2 
pcmin_fullterm.py and cape.py written by Shuai Wang, Imperial College London
Adapted by Shi Wei Yuan and Neil Patel, Imperial College London
Supervisor: Professor Ralf Toumi

Aims:
The code files in this map are used to generate MPI fields using SST, MSLP, T and Q (SH) data from ERA5. To download these data files:

Specific Humidity and Temperature
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=form

SST and MSLP
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form

Before running the modules, create these empty folders in the same directory as the modules:

“GEN_MPI”		will contain MPI fields generated from distribution
“MEAN_STD”		will contain mean and std per grid point across 10 years for reanalysed MPI fields
“MPI_MAPS”		will contain reanalysed MPI fields
“MSLP_FIELDS”		will contain MSLP ERA5 data in .txt
“Q_FIELDS”		will contain SH/Q ERA5 data in .txt
“SST_FIELDS”		will contain SST ERA5 data in .txt
“T_FIELDS”		will contain T ERA5 data in .txt

Please run the modules in the following order:

1. cape.py
2. pcmin_fullterm.py
3. preprocessing.py
4. import_fields.py  —> takes quite some time!
5. create_mpi_map.py —> takes hours to run, perhaps only run for 1 month, basin first to check running time
6. generate_mpi.py
