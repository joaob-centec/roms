#!/bin/bash

#
# Script to process ECMWF forcing on an yearly basis.
# Since the transfer of accumulated fields to instantaneous needs the 
# field for the next time step, to process a month we need the fields of 
# the following month.
#
#
# Joao H Bettencourt, Jan 2021

# Parameters

DOWNLOAD=Download_ERA5_CDS_GRIB_Naval.py

TRANSFER=ECMWF_Forcing_GRIB2ASCII_Naval1181.py

WRITE=ROMS_forcing.py

REMOVE=/bin/rm

FORCING_DIR=/mnt/disk1/Observa.FISH/Forcing/ECNA124U

# Setup conda

CONDA_ENV=roms

CONDA_SOURCE=/home/joao/miniconda3/etc/profile.d/conda.sh

PYTHON=/home/joao/miniconda3/envs/roms/bin/python

source $CONDA_SOURCE

conda activate $CONDA_ENV

# Process forcing

cd $FORCING_DIR

YEAR=1999
# January
echo "Processing January $YEAR"
#$PYTHON $DOWNLOAD $YEAR 1 $YEAR 2
#$PYTHON $TRANSFER $YEAR 1 $YEAR 1
#$PYTHON $WRITE $YEAR 1
#$REMOVE *.dat
# February
echo "Processing February $YEAR"
#$PYTHON $DOWNLOAD $YEAR 3 $YEAR 3
#$PYTHON $TRANSFER $YEAR 2 $YEAR 2
#$PYTHON $WRITE $YEAR 2
#$REMOVE *.dat
# March 
echo "Processing March $YEAR"
#$PYTHON $DOWNLOAD $YEAR 4 $YEAR 4
#$PYTHON $TRANSFER $YEAR 3 $YEAR 3
#$PYTHON $WRITE $YEAR 3
#$REMOVE *.dat
# April
echo "Processing April $YEAR"
#$PYTHON $DOWNLOAD $YEAR 5 $YEAR 5
#$PYTHON $TRANSFER $YEAR 4 $YEAR 4
#$PYTHON $WRITE $YEAR 4
#$REMOVE *.dat
# May
echo "Processing May $YEAR"
#$PYTHON $DOWNLOAD $YEAR 6 $YEAR 6
#$PYTHON $TRANSFER $YEAR 5 $YEAR 5
#$PYTHON $WRITE $YEAR 5
#$REMOVE *.dat
# June 
echo "Processing June $YEAR"
#$PYTHON $DOWNLOAD $YEAR 7 $YEAR 7
#$PYTHON $TRANSFER $YEAR 6 $YEAR 6
#$PYTHON $WRITE $YEAR 6
#$REMOVE *.dat
# July
echo "Processing July $YEAR"
#$PYTHON $DOWNLOAD $YEAR 8 $YEAR 8
#$PYTHON $TRANSFER $YEAR 7 $YEAR 7
#$PYTHON $WRITE $YEAR 7
#$REMOVE *.dat
# August
echo "Processing August $YEAR"
#$PYTHON $DOWNLOAD $YEAR 9 $YEAR 9
#$PYTHON $TRANSFER $YEAR 8 $YEAR 8
#$PYTHON $WRITE $YEAR 8
#$REMOVE *.dat
# September
echo "Processing September $YEAR"
#$PYTHON $DOWNLOAD $YEAR 10 $YEAR 10
#$PYTHON $TRANSFER $YEAR 9 $YEAR 9
#$PYTHON $WRITE $YEAR 9
#$REMOVE *.dat
# October
echo "Processing October $YEAR"
#$PYTHON $DOWNLOAD $YEAR 11 $YEAR 11
#$PYTHON $TRANSFER $YEAR 10 $YEAR 10
#$PYTHON $WRITE $YEAR 10
#$REMOVE *.dat
# November
echo "Processing November $YEAR"
$PYTHON $DOWNLOAD $YEAR 12 $YEAR 12
#$PYTHON $TRANSFER $YEAR 11 $YEAR 11
#$PYTHON $WRITE $YEAR 11
#$REMOVE *.dat
# December
echo "Processing December $YEAR"
$PYTHON $DOWNLOAD $((YEAR+1)) 1 $((YEAR+1)) 1
$PYTHON $TRANSFER $YEAR 12 $YEAR 12
$PYTHON $WRITE $YEAR 12
$REMOVE *.dat

 
 
