# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# ROMS_forcing.py: Run ROMS forcing software to process bulk forcing data
# obtained from ECMWF's ERA-5 atmospheric reanalysis in text format
# and write it in NetCDF format to be used by the ROMS model.
#
# The input fields are stored in ASCII format.
# For details on this format, see:
#
#      http://www.myroms.org/index.php?page=forcing
# 
# The fields that are needed to force ROMS are the following:
#
# Table 1.
# 
# ID Index   Forcing Field 	        Variable Time 	   Units
# -----------------------------------------------------------
#  1 idUsms  surface U-momentum stress 	sustr 	 sms_time  Pa (N/m2)
#  2 idVsms  surface V-momentum stress 	svstr 	 sms_time  Pa (N/m2)
#  3 idUair  surface U-wind component 	Uair 	 wind_time m/s
#  4 idVair  surface V-wind component 	Vair 	 wind_time m/s
#  5 idWamp  wind-induced wave amp. 	Awave 	 wave_time m
#  6 idWdir  wind-induced wave dir. 	Dwave 	 wave_time degree
#  7 idWper  wind-induced wave Per. 	Pwave 	 wave_time s
#  8 idSden  bot sed grain density 	Ssize 	 - 	   m
#  9 idSsiz  bot sed grain diameter 	Sdens 	 - 	   kg/m3
# 10 idTdew  surf dew-point temp 	Tdew 	 tdew_time Kelvin
# 11 idTair  surf air temp       	Tair 	 tair_time Celsius
# 12 idPair  surf air press     	Pair 	 pair_time mb
# 13 idQair  surf air rel humidity 	Qair 	 qair_time percentage
# 14 idCfra  cloud fraction 	        cloud 	 cloud_time0 to 1
# 15 idrain  rain fall rate 	        rain 	 rain_time kg/m2/s
# 16 idSrad  net sw radiation flux 	swrad 	 srf_time  Watts/m2
# 17 idLrad  net lw radiation flux 	lwrad 	 lrf_time  Watts/m2
# 18 idShea  net sens heat flux 	sensible sen_time  Watts/m2
# 19 idfSSS  sea surf salinity 	        SSS 	 sss_time  PSU
# 20 idfSST  sea surf temperature 	SST 	 sst_time  Celsius
# 21 iddQdT  surf net heat flx sens SST dQdSST   sst_time  Watts/m2/Celsius
# 22 idTsur(itemp) net surf heat flux 	shflux 	 shf_time  Watts/m2
# 23 idTsur(isalt) net surf fhwter flux swflux 	 swf_time  cm/day 
# -----------------------------------------------------------
#
# Bulk forcing fields are 3, 4, 11, 12, 13, 14, 15, 16. 
#
# This script does the following:
#
#    1. Writes the forcing code input file 
#    2. Runs the forcing code
#
# The script processes one month of data each time it is called. The year and month
# are defined in the calling arguments.
#
#
#
#
# Revision list:
#
# Date          Changes                        Author
# ----          -------                        ------
# 2021/01       Created.                       Joao H. Bettencourt
#
#
# Module imports.

import numpy as np
import datetime
import glob
import argparse
import subprocess

# Parse input arguments

parser = argparse.ArgumentParser()

parser.add_argument("year0", help="starting year of data retrieval",
                    type=int)
parser.add_argument("month0", help="starting month of data retrieval",
                    type=int)
args = parser.parse_args()

#
# Define forcing executable
#

forcing="/home/joao/forcing/ECNA124U/forcing"

#
# Define forcing input file entries
#

job={
      "NAME":"ECNA124U ERA5 Bulk forcing",        # Job name
      "JOB":0,                       # JOB    : Processing job type: [0] create, [1] append               
      "IOTYPE":1,                   # IOTYPE : type of output floating-point data: [0] single, [1] double.  
      "INTOPT":0,                    # INTOPT : interpolation option: [0] linear, [1] Bessel.
      "GRIDFILE":"../../Malhas/ECNA124U/ecna124u_grd_202008191014_sponge.nc",  # ROMS grid file
      "FORCINGFILEHEADER":"ecna124u_bulk", # Forcing file header. Forcing file name is "FORCINGFILEHEADER+_YYYYMMM.nc"
      "LOPEN":"F",                    # Open input file with SEVERAL forcing fields, if any?
      "LOPENNAME":""         # If LOPEN=T, this is the name of the file to be read.
}

Lout={                        # Switch to process field [T/F].
      "idVsms":"F",
      "idVair":"T",
      "idWamp":"T",
      "idSden":"F",
      "idTdew":"T",
      "idTair":"T",
      "idPair":"T",
      "idQair":"T",
      "idCfra":"T",
      "idrain":"T",
      "idSrad":"T",
      "idLrad":"F",
      "idShea":"F",
      "idfSSS":"F",
      "idfSST":"T",
      "iddQdT":"F",
      "idTsur":"F",
      "idTswr":"F" 
}

Fdtype={                       # data type: [0] point, [1] grided.
      "idVsms":1,
      "idVair":1,
      "idWamp":1,
      "idSden":1,
      "idTdew":1,
      "idTair":1,
      "idPair":1,
      "idQair":1,
      "idCfra":1,
      "idrain":1,
      "idSrad":1,
      "idLrad":1,
      "idShea":1,
      "idfSSS":1,
      "idfSST":1,
      "iddQdT":1,
      "idTsur":1,
      "idTswr":1 
}

Fifrmt={                       # input data format: [1] ascii, [2] NetCDF.
      "idVsms":1,
      "idVair":1,
      "idWamp":1,
      "idSden":1,
      "idTdew":1,
      "idTair":1,
      "idPair":1,
      "idQair":1,
      "idCfra":1,
      "idrain":1,
      "idSrad":1,
      "idLrad":1,
      "idShea":1,
      "idfSSS":1,
      "idfSST":1,
      "iddQdT":1,
      "idTsur":1,
      "idTswr":1 
}

Fscale={                       # Scale factor.
      "idVsms":1.0,
      "idVair":1.0,
      "idWamp":1.0,
      "idSden":1.0,
      "idTdew":1.0,
      "idTair":1.0,
      "idPair":1.0,
      "idQair":1.0,
      "idCfra":1.0,
      "idrain":1.0,
      "idSrad":1.0,
      "idLrad":1.0,
      "idShea":1.0,
      "idfSSS":1.0,
      "idfSST":1.0,
      "iddQdT":1.0,
      "idTsur":1.0,
      "idTswr":1.0 
}

Fcycle={                       # Length (days) of time cycle.
      "idVsms":0.0,
      "idVair":0.0,
      "idWamp":0.0,
      "idSden":0.0,
      "idTdew":0.0,
      "idTair":0.0,
      "idPair":0.0,
      "idQair":0.0,
      "idCfra":0.0,
      "idrain":0.0,
      "idSrad":0.0,
      "idLrad":0.0,
      "idShea":0.0,
      "idfSSS":0.0,
      "idfSST":0.0,
      "iddQdT":0.0,
      "idTsur":0.0,
      "idTswr":0.0 
}

Frot={                       # Rotate input fields to specified grid [T/F].
      "idVsms":"T",
      "idVair":"T",
      "idWamp":"",
      "idSden":"",
      "idTdew":"",
      "idTair":"",
      "idPair":"",
      "idQair":"",
      "idCfra":"",
      "idrain":"",
      "idSrad":"",
      "idLrad":"",
      "idShea":"",
      "idfSSS":"",
      "idfSST":"",
      "iddQdT":"",
      "idTsur":"",
      "idTswr":"" 
}


Lopen={                        # Open input file?.
      "idVsms":"F",
      "idVair":"T",
      "idWamp":"T",
      "idSden":"F",
      "idTdew":"T",
      "idTair":"T",
      "idPair":"T",
      "idQair":"T",
      "idCfra":"T",
      "idrain":"T",
      "idSrad":"T",
      "idLrad":"F",
      "idShea":"F",
      "idfSSS":"F",
      "idfSST":"T",
      "iddQdT":"F",
      "idTsur":"F",
      "idTswr":"F" 
}

Lopenheader={                    # Open input file header. Input file name is "Lopenheader_YYYYMM.dat"
      "idVsms":"ERA5_cds_SMS",
      "idVair":"ERA5_cds_SWIND",
      "idWamp":"ERA5_cds_WAVE",
      "idSden":"",
      "idTdew":"ERA5_cds_TDEW",
      "idTair":"ERA5_cds_TAIR",
      "idPair":"ERA5_cds_PAIR",
      "idQair":"ERA5_cds_QAIR",
      "idCfra":"ERA5_cds_CLOUD",
      "idrain":"ERA5_cds_RAIN",
      "idSrad":"ERA5_cds_SWRAD",
      "idLrad":"ERA5_cds_LWRAD",
      "idShea":"ERA5_cds_SENSIBLE",
      "idfSSS":"",
      "idfSST":"ERA5_cds_SST",
      "iddQdT":"",
      "idTsur":"ERA5_cds_SHFLUX",
      "idTswr":"ERA5_cds_SWFLUX" 
}

#
# Write forcing file
#

yearMonth="_{0:4d}{1:02d}".format(args.year0,args.month0)

fFileName="forcing_bulk" + yearMonth

fFileId=open(fFileName,'w')

# Job parameters

fFileId.write(job["NAME"] + "\n")
fFileId.write("{:1d}".format(job["JOB"]) + "\n")
fFileId.write("{:1d}".format(job["IOTYPE"]) + "\n")
fFileId.write("{:1d}".format(job["INTOPT"]) + "\n")

# Field parameters

for key in Lout:

    fFileId.write(Lout[key]+"\n")
    fFileId.write("{:1d}".format(Fdtype[key]) + "\n")
    fFileId.write("{:1d}".format(Fifrmt[key]) + "\n")
    fFileId.write("{:<5.1f}".format(Fscale[key]) + "\n")
    if key != "idSden":
        fFileId.write("{:<5.1f}".format(Fcycle[key]) + "\n")
    if key == "idVsms" or key == "idVair":

      fFileId.write(Frot[key]+"\n")

# File parameters 

fFileId.write(job["GRIDFILE"] + "\n")
fFileId.write(job["FORCINGFILEHEADER"] + yearMonth + ".nc\n")
fFileId.write(job["LOPEN"] + "\n")
if job["LOPEN"] == "T":
  fFileId.write(job["LOPENNAME"] + "\n")
else:
  fFileId.write("/dev/null\n")

for key in Lopen:

    fFileId.write(Lopen[key] + "\n")
    if Lopen[key] == "T":
      fFileId.write(Lopenheader[key] + yearMonth + ".dat\n")
    else:
      fFileId.write("/dev/null\n")
       

# Close input file
fFileId.close()

#
# Run forcing 
#

inFileId = open(fFileName,'r')

outFile = "forcing" + yearMonth + ".out"

outFileId = open(outFile,'w')

subprocess.run(forcing, stdin=inFileId, stdout=outFileId, check=True)

outFileId.close()


