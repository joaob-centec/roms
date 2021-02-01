# roms
Scripts to do pre/post processing of ROMS model input/output:
1. Forcing/Download_ERA5_CDS_GRIB.py: this script downloads ROMS forcing fields from the ERA5 Reanalysis archive in GRIB format. The call sign of the script is python Download_ERA5_CDS_GRIB.py year0 month0 year1 month1, where year0, month0 are the first year and first month for which data is wanted and year1, month1 are the last year and month for which data is wanted. Limits of susbsetting area must be defined inside the script.

2. Forcing/ECMWF_Forcing_GRIB2ASCII.py: this script reads the GRIB files from the ERA5 archive and writes them in an intermediate ASCII format to be later processed by the ROMS Fortran forcing package. Data manipulations necessary before conversion to ROMS forcing field input are also performed by the script (e.g. unit conversion, transformation from accumulated to instantaneous fields). The script is a bit clumsy and difficult to follow, but AFAK it is relatively bug free.

3. Forcing/ROMS_forcing.py: this script writes the input file for the ROMS forcing package code. It uses dictionaries to configure the package's processing of the ASCII files produced by the ECMWF_Forcing_GRIB2ASCII.py script.

4. Forcing/Process_ECMWF_Forcing.sh: bash script to run the 3 scripts above for a single year. 

