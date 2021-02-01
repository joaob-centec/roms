# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Process_ECMWF_Forcing.py: script to process ROMS forcing fields
# obtained from ECMWF's ERA-5 atmospheric reanalysis.
#
# ERA-5 fields are processed to the format that the ROMS forcing
# package specifies. The process fields are stored in ASCII format.
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
# Variables downloaded from ERA-5 CDS archive are:
#
# Table 2:
# -------
# 
# ROMS Var   CDS Var                           Code     Units     I/A
# ---------------------------------------------------------------------
# sustr      instantaneous_eastward_turbulent_ 128.229  N m^-2    I
#            surface_stress           
# svstr      instantaneous_northward_          128.230  N m^-2    I
#            turbulent_surface_stress
# Uair       10m_u_component_of_wind           128.165  m s^-1    I
# Vair       10m_v_component_of_wind           128.166  m s^-1    I
# Awave      significant_height_of_combined_   140.229  m         I
#            wind_waves_and_swell  
# Dwave      mean_wave_direction               140.230  degree    I
# Pwave      mean_wave_period                  140.232  s         I
# Ssize      N/A
# Sdens      N/A
# Tdew       2m_dewpoint_temperature           128.168  K         I
# Tair       2m_temperature                    128.167  K         I
# Pair       surface_pressure                  128.134  Pa        I
# Qair       N/A (1)
# cloud      total_cloud_cover                 128.164  (0-1)     I
# rain       total_precipitation               228      m         A
# swrad      surface_net_solar_radiation       176      J m^-2    A
# lwrad      surface_net_thermal_radiation     177      J m^-2    A
# sensible   surface_sensible_heat_flux        128.146  J m^-2    A
# SSS        N/A
# SST        sea_surface_temperature           128.34   K         I
# dQdSST     N/A
# shflux     N/A (2)
# swflux     N/A (3)
# ---------------------------------------------------------------------
# I - Instantaneous; A - Accumulated.
#
# (1) Qair is computed from the 2-meter dewpoint temperature, air 
# temperatur and surface pressure. See https://confluence.ecmwf.int/
# pages/viewpage.action?pageId=82870405#ERA5:datadocumentation-
# Computationofnear-surfacehumidityandsnowcover 
#
# (2) The surface net heat flux is:
# net heat = sensible + latent + net solar + net thermal. 
# The latent heat flux is thus downloaded as:
#        
#            surface_latent_heat flux          147     J m^-2     A 
#
# to be used in the calculation of the surface net heat flux.  
#
# (3) The surface net freshwater flux is:  
# net freshwater = precipitation + evaporation
# The evaporation is thus downloaded as:        
#
#            evaporation                      182      m         A 
#
# to be used in the calculation of the surface net heat flux. 
# 
#
# The present script performs two things:
#
# 1. Manipulates the ECMWF fields to put them in the format ROMS needs,
#    e.g., unit conversion, accumulated to instantaneous etc (see notes
#    to table 2). For details on the field manipulations, see comments
#    on the relevant section of the code.
#
# 2. Transfer the fields from the ECMWF netcdf format to ASCII format.
#
# Input ASCII format (https://www.myroms.org/index.php?page=forcing):
#
# The input ASCII format for each forcing field is as follows:
#
#      NTIME(id....)   Length of time series (number of records).
#      NLON (id....)   Number of longitude points.
#      NLAT (id....)   Number of latitude points.
#      TIME            Time of the data (days, Julian days).
#      LON             Longitude of the grided data, LON(NLON).
#      LAT             Latitude  of the grided data, LAT(NLAT).
#      VAL             Scaler forcing field, VAL(NLON*NLAT).
#      VALX            X-component forcing field, VALX(NLON*NLAT).
#      VALY            Y-component forcing field, VALY(NLON*NLAT).
#
# Grided Time-Series Data
#         read (iunit(id....),*) NTIME(id....),
#     &                          NLON (id....), NLAT (id....)
#         read (iunit(id....),*) (LON(i), i = 1, NLON(id....))
#         read (iunit(id....),*) (LAT(j), j = 1, NLAT(id....))
#         do n = 1, NTIME(id....)
#           read (iunit(id....),*) TIME
#           NP = NLON(id....) * NLAT(id....)
#           read (iunut(id....),*) (VAL(i), i = 1, NP)      scalar data
#      or
#           read (iunit(id....),*) TIME
#           NP = NLON(id....) * NLAT(id....)
#           read (iunit(id....),*) (VALX(i), i = NP)
#                                                           vector data
#           read (iunit(id....),*) TIME
#           NP = NLON(id....) * NLAT(id....)
#           read (iunit(id....),*) (VALY(i), i = NP)
#        enddo
#  
# Notice that all reading is free formatted and the input grid must 
# be ascending order:
#
#         LON(i) < LON(i+1)  and   LAT(j) < LAT(j+1)
#  
# and data must not have flagged values. Otherwise, the interpolation will 
# not work correctly!. Users can modify this input ASCII format easily. 
# The one provided above is just a guideline. Forcing Fields 
#
# The script does not interpolate from the ECMWF grid to the ROMS grid.
# This interpolation is performed in the ROMS forcing package.
#
#
# Revision list:
#
# Date          Changes                        Author
# ----          -------                        ------
# 2020/10       Created.                       Joao H. Bettencourt
#
#

# Module imports.

import numpy as np
import datetime
import glob
import pygrib

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("year0", help="starting year of data retrieval",
                    type=int)
parser.add_argument("month0", help="starting month of data retrieval",
                    type=int)
parser.add_argument("year1", help="ending year of data retrieval",
                    type=int)
parser.add_argument("month1", help="ending month of data retrieval",
                    type=int)
args = parser.parse_args()


# Functions

def Accu2Inst(fieldNext, fieldNow, dt):

    """ This function computes the instantaneous field Fi at time t 
    from the accumulated fields Fa at times t, t+dt, where dt
    is the time step between fields in seconds. """

#   The instantaneous field at t is computed as the average of the
#   accumulations at t and t+dt divided by dt:
#
#   Fi(t) = 1/(2*dt)*(FieldNext + FieldNow).
#
#   If FieldNext doesn't exist use FieldNow/dt.
#
#   Arguments:
#
#   fieldNext     : Fa(t+dt)
#   fieldNow      : Fa(t)
#   dt            : time step in seconds
#   
#   fieldI        : Fi 
#

    cff1 = 1.0/(2.0*dt)

    if (fieldNext == -1 ).all(): # Fa(t+dt) not available. 

      fieldI = 2.0 * cff1 * fieldNow

    else: # Use average. 

      fieldI = cff1 * (fieldNext + fieldNow)

    return fieldI

def surfRelHumidity( d2, t2):

    """ This function computes the surface relative humidity
    from the 2m dew point temperature, 2m air temperature
    and surface pressure. """
#
# As of June 2020, specific humidity is not available for retrieval from 
# the ERA Interim data server (see note by Coquart Laure, Jan 28 2020, 
# available at:
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214).
#
# To compute relative humidity we use the formula:
#
#  RH = 100 * es(Td)/es(T), where Td and T are the 2 m dew point
#  and air temperatures, and es is the saturation water vapour pressure 
# (documentation of the IFS for CY41R2, Ch. 7, sec. 7.2.1b, eq. 7.5  ):
#
#  es(T) = a1 * exp[a3*(T-T0)/(T-a4)],
#  with a1=611.21 Pa, a3=17.502, a4=32.19 K and T0=273.16 K. 

#   Set parameters

    T0 = 273.16

    a1 = 611.21
    a3 = 17.502
    a4 = 32.19

#   Saturation water vapour pressure

    esat_Td = a1 * np.exp(a3 * (d2 - T0)/(d2 - a4))

    esat_T = a1 * np.exp(a3 * (t2 - T0)/(t2 - a4))

#   Relative humidity

    relh = esat_Td/esat_T

    return relh

def writeASCIIHeader(fId, nTime, nLon, nLat, lon, lat):

    """ This function writes the header of the ASCII file 
    to be read by the ROMS forcing processing package. """

#   Data sizes

    fId.write('{0:<4d} {1:<4d} {2:<4d}\n'.format(nTime, nLon, nLat))

#   Latitudes and longitudes as numpy arrays using tofile method

    lon.tofile(fId, sep=" ", format="%9.4f")
    fId.write('\n')
    lat.tofile(fId, sep=" ", format="%9.4f")

    return 0

def writeASCIIField(fId, time, field):

    """ This function writes the data field at time time 
    to the ASCII file to be read by the ROMS forcing 
    processing package. """

#   Field time (days)

    fId.write('\n{0:<f}\n'.format(time))

#   Field data

    if isinstance(field, np.ma.MaskedArray):
    
      field.data.tofile(fId, sep=" ", format="%E")

    else:
         
      field.tofile(fId, sep=" ", format="%E")

    return 0

# EMCWF variable information.  Dictionaries keys are the variables
# names as define in Table 1.

Variables = { # Use these to get the grib message  
"sustr":"instantaneous_eastward_turbulent_surface_stress",           
#"svstr":"instantaneous_northward_turbulent_surface_stress",
"Uair":"10m_u_component_of_wind",         
#"Vair":"10m_v_component_of_wind",        
"Awave":"significant_height_of_combined_wind_waves_and_swell",  
#"Dwave":"mean_wave_direction",              
#"Pwave":"mean_wave_period",                
"Ssize":"N/A",
"Sdens":"N/A",
"Tdew":"2m_dewpoint_temperature",        
"Tair":"2m_temperature",                
"Pair":"surface_pressure",             
"Qair":"surface_air_relative_humidity", # Dummy name. Won't be used.
"cloud":"total_cloud_cover",           
"rain":"total_precipitation",       
"swrad":"surface_net_solar_radiation",    
"lwrad":"surface_net_thermal_radiation", 
"sensible":"surface_sensible_heat_flux",   
"SSS":"N/A",
"SST":"sea_surface_temperature",     
"dQdSST":"N/A",
"shflux":"surface_latent_heat_flux",
"swflux":"evaporation"
	}

Short_names = { # We use these to get input and output file names
"sustr":"iews",           
"svstr":"inss",
"Uair":"10u",         
"Vair":"10v",        
"Awave":"swh",  
"Dwave":"mwd",              
"Pwave":"mwp",                
"Ssize":"",
"Sdens":"",
"Tdew":"2d",        
"Tair":"2t",                
"Pair":"sp",             
"Qair":"2d", # Here 2d is a dummy variable. The Qair is computed with a formula.
"cloud":"tcc",           
"rain":"tp",       
"swrad":"ssr",    
"lwrad":"str", 
"sensible":"sshf",   
"SSS":"",
"SST":"sst",     
"dQdSST":"",
"shflux":"slhf",
"swflux":"e"
}

Accu_type = { # (i)nstantaneous or (a)ccumulated field
"sustr":"i",           
"svstr":"i",
"Uair":"i",         
"Vair":"i",        
"Awave":"i",  
"Dwave":"i",              
"Pwave":"i",                
"Ssize":"",
"Sdens":"",
"Tdew":"i",        
"Tair":"i",                
"Pair":"i",             
"Qair":"i",
"cloud":"i",           
"rain":"a",       
"swrad":"a",    
"lwrad":"a", 
"sensible":"a",   
"SSS":"",
"SST":"i",     
"dQdSST":"",
"shflux":"a",
"swflux":"a"
}

# Output files

outputHeader = 'ERA5_cds_'
outputTail   = '.dat'

# Years/months to process

yearStart = args.year0
yearEnd = args.year1
monthStart = args.month0
monthEnd = args.month1

# Averaging period

avePeriod = 3 # With respect to the time step of the download data, e.g., if hourly data then 
               # avePeriod = n provides averages every n hours

avePeriodDays = avePeriod/24 # Must make sure this is in days

sec2Day = 86400

Years = range(yearStart,yearEnd + 1)

# Set reference date

refDate = datetime.datetime(1950,1,1)

# Set 0ºC in Kelvin

Temp0 = 273.16;

# Loop Years/months and process data

for year in Years:

    print(" ")
    print("Process year: " + str(year))


#   Compute the first and last month of the year to process

    if Years.index(year) == 0 and Years.index(year) == len(Years) - 1 :

       month0 = monthStart
       month1 = monthEnd + 1

    elif Years.index(year) == 0  and Years.index(year) != len(Years) - 1 :

       month0 = monthStart # monthStart is the 1st month to be processed of yearStart
       month1 = 13

    elif Years.index(year) != 0  and Years.index(year) == len(Years) - 1 :

       month0 = 1
       month1 = monthEnd + 1 # monthEnd is the last month to be processed of yearEnd

    Months = range(month0,month1)

#   Process months

    for month in Months:

        print("  Process month: " + str(month))

#   Set the YYYYMM code string of the previous, present and next months

        presentCode = '{0:4d}{1:02d}'.format(year,month)

        print("   present code: " + presentCode)

        if month == 1 : # If 1st month, last month of previous year

           previousCode = '{0:4d}12'.format(year - 1)
           print("   previous code: " + previousCode)

        else:

           previousCode = '{0:4d}{1:02d}'.format(year,month-1)
           print("   previous code: " + previousCode)

        if month == 12 : # If last month, 1st month of next year

           nextCode = '{0:4d}01'.format(year + 1)
           print("   next code: " + nextCode)

        else:

           nextCode = '{0:4d}{1:02d}'.format(year,month+1)
           print("   next code: " + nextCode)

# Process monthly fields
                  
        for key in Variables:

            varTitle = " ".join(Variables[key].split('_')).swapcase()

            # Initialize temporary storage 

            #dewT = []
            #airT = []
            insTime = []
    #        raiPrevious = []
    #        raiNext = []
    #        swrPrevious = []
    #        swrNext = []
    #        lwrPrevious = []
    #        lwrNext = []
    #        senPrevious = []
    #        senNext = []
    #        latPrevious = []
    #        latNext = []
    #        evaPrevious = []
    #        evaNext = []

            if Variables[key]=="N/A":

               print(" EMCWF ERA 5 does not provide field for " + key)

            #elif Variables[key] == "C":

             #   print(" Variables " + varTitle + " to be computed from other variables" )

            else:

               print(" Processing " + varTitle)

               inName = outputHeader + Short_names[key].upper() + "_" + presentCode + ".grib"

               inGrib = pygrib.open(inName)  

               nMsgs = inGrib.messages

               print("    Messages in grib file: " + str(nMsgs))

               nFlds = nMsgs / avePeriod 

               # Get lat/lon of grib field

               inGrib.seek(0) # Go to first message

               grb = inGrib.read(1)[0] # Get first message

               lats, lons = grb.latlons()

               lats = np.flipud(lats) # ECMWF saves data north to south

               nLats, nLons = lats.shape

               # Prepare averaging aux variables

               aveStep = np.mod(np.arange(nMsgs) + 1, avePeriod) 

               aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

               

               if key == "sustr":

                  # Open y component of surface momentum stress

                  print(" Processing y component")

                  inName2 = outputHeader + "inss".upper() + "_" + presentCode + ".grib"
 
                  inGrib2 = pygrib.open(inName2)  

                  inGrib2.seek(0)
 
                  aveFieldAccum2 = np.zeros(lats.shape) # Accumulate fields here

                  outName = outputHeader + "SMS_" + presentCode + outputTail

                  outId = open(outName,'w')

                  # Print output file header

                  writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

               elif key == 'Uair':

                  # Open y component of surface wind velocity

                  print(" Processing y component")

                  inName2 = outputHeader + "10v".upper() + "_" + presentCode + ".grib"
 
                  inGrib2 = pygrib.open(inName2)  

                  inGrib2.seek(0)
 
                  aveFieldAccum2 = np.zeros(lats.shape) # Accumulate fields here

                  outName = outputHeader + "SWIND_" + presentCode + outputTail

                  outId = open(outName,'w')

                  # Print output file header

                  writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

               elif key == 'Awave':

                  # Open mean wave direction

                  print(" Processing mean wave direction ")

                  inName2 = outputHeader + "mwd".upper() + "_" + presentCode + ".grib"
 
                  inGrib2 = pygrib.open(inName2)  

                  inGrib2.seek(0)
 
                  aveFieldAccum2 = np.zeros(lats.shape) # Accumulate fields here

                  inName3 = outputHeader + "mwp".upper() + "_" + presentCode + ".grib"
 
                  inGrib3 = pygrib.open(inName3)  

                  inGrib3.seek(0)
 
                  aveFieldAccum3 = np.zeros(lats.shape) # Accumulate fields here

                  outName = outputHeader + "WAVE_" + presentCode + outputTail

                  outId = open(outName,'w')

                  # Print output file header

                  writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

               else:
               
                  if key == 'Tdew':
               
                    D2=[]
                  
                  elif key == 'Tair':
                  
                    T2=[]

                  outName = outputHeader + key.upper() + "_" + presentCode + outputTail

                  outId = open(outName,'w')

                  # Print output file header

                  writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

             
               # Cycle message times and process

               inGrib.seek(0)

               if Accu_type[key] == "i":

                  for j in range(nMsgs):

                      grb = inGrib.message(j+1)


                      msgTime = datetime.datetime.strptime('{0:d}{1:04d}'.format(
                                 grb.validityDate,grb.validityTime),
                                 '%Y%m%d%H%M') #- refDate

                      #print("       Message " + str(j+1) + 
                      #        " validityDate: " + '{0:d}H{1:04d}'.format(
                      #           grb.validityDate,grb.validityTime) +
                      #        " Time: " + str(msgTime.total_seconds()))  

                      if aveStep[j] == 1:

                         aveStartTime = msgTime

                      # Unit conversion or further processing
                      if key == "sustr" or key == "Uair":

                        outField = grb.values

                        outField2 = inGrib2.message(j+1)

                        aveFieldAccum2 = aveFieldAccum2 + outField2.values

                      if key == "Awave":

                        outField = grb.values.data

                        outField2 = inGrib2.message(j+1).values.data

                        aveFieldAccum2 = aveFieldAccum2 + outField2

                        outField3 = inGrib3.message(j+1).values.data

                        aveFieldAccum3 = aveFieldAccum3 + outField3

                      elif key == "Tdew" :

                        outField = grb.values

                        D2.append(grb.values)#dewT=grb.values # Saving for rel hum calculation
                        #print( "airDew=" + str(dewT[44,84]))

                      elif key == "Tair" :

                        outField = grb.values - Temp0 # Kelvin to centigrades

                        T2.append(grb.values)#airT=grb.values # Saving for rel hum calculation
                        #print( "airT=" + str(airT[44,84]))

                      elif key == "Qair":
                         #print( "airDew=" + str(dewT[44,84]))
                         #print( "airT=" + str(airT[44,84]))
                         outField =  surfRelHumidity( D2[j],T2[j] )#dewT, airT)
                         #print( "qair=" + str(outField[44,84]))

                      elif key == "Pair" :

                         outField = grb.values * 0.01 # Pascal to mb

                      elif key == "SST" :

                         outField = grb.values.data - Temp0  

                      else:

                         outField = grb.values


                      aveFieldAccum = aveFieldAccum + outField

                      insTime.append(msgTime)
                      

                      # Compute averages and write averaged field  
                      
                      if aveStep[j] == 0:


                         aveEndTime = msgTime

                         outField = aveFieldAccum / avePeriod
                         aveTime = aveStartTime + datetime.timedelta(days=0.5*avePeriodDays)


                         outTime = aveTime - refDate

                         outTimeDays = outTime.days + outTime.seconds/sec2Day

                         print("Averages from " + str(aveStartTime) +
                                " to " + str(aveEndTime) + " at " + str(aveTime) + 
                         " (" + str(outTimeDays) + " days since " + 
                         str(refDate) + ")")

                         writeASCIIField(outId, outTimeDays, np.flipud(outField))

                         aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here
                         
                         if key == "sustr" or key == "Uair":

                            outField2 = aveFieldAccum2 / avePeriod

                            writeASCIIField(outId, outTimeDays, np.flipud(outField2))

                            aveFieldAccum2 = np.zeros(lats.shape) # Accumulate fields here

                         elif key == "Awave":

                            outField2 = aveFieldAccum2 / avePeriod

                            writeASCIIField(outId, outTimeDays, np.flipud(outField2))

                            aveFieldAccum2 = np.zeros(lats.shape) # Accumulate fields here

                            outField3 = aveFieldAccum3 / avePeriod

                            writeASCIIField(outId, outTimeDays, np.flipud(outField3))

                            aveFieldAccum3 = np.zeros(lats.shape) # Accumulate fields here

               elif Accu_type[key] == "a": # Here we deal with accumulated fields

                   for j in range(nMsgs):


                      if j == 0:

                          grb = inGrib.message(1)

                         # Load next message

                          grbNext =  inGrib.message(2)

                      elif j == nMsgs - 1:

                          grb = inGrib.message(j+1)
 
                         # Load next month to get first message
                        
                          inTempName = outputHeader + Short_names[key].upper() + "_" + nextCode + ".grib"

                          inTempGrib = pygrib.open(inTempName)  

                          grbNext = inTempGrib.message(1)

                      else:  

                          grb = inGrib.message(j+1)

                          # Load next message

                          grbNext =  inGrib.message(j+2)


                      msgTime = datetime.datetime.strptime('{0:d}{1:04d}'.format(
                                 grb.validityDate,grb.validityTime),
                                 '%Y%m%d%H%M')
                      msgTimeNext = datetime.datetime.strptime('{0:d}{1:04d}'.format(
                                 grbNext.validityDate,grbNext.validityTime),
                                 '%Y%m%d%H%M')

                      dt = msgTimeNext - msgTime
                      
                      if aveStep[j] == 1:

                         aveStartTime = msgTime

                      # Unit conversions and temporary storage of data

                      if key == "rain":

                         outField = grb.values * 1000
                         outFieldNext = grbNext.values * 1000

                         raiNow = outField
                         raiNext=grbNext.values * 1000

                      elif key == "swrad":

                         outField = grb.values 
                         outFieldNext = grbNext.values

                         swrNow=outField
                         swrNext=grbNext.values

                      elif key == "lwrad":

                         outField = grb.values 
                         outFieldNext = grbNext.values

                         lwrNow=outField
                         lwrNext=grbNext.values

                      elif key == "sensible":

                         outField = grb.values
                         outFieldNext = grbNext.values

                         senNow=outField
                         senNext=grbNext.values

                      elif key == "shflux":

                         #outFieldPrevious = grbPrevious.values
                         #outFieldNext = grbNext.values

                         outField = grb.values 
                         
                         latNow = outField
                         latNext=grbNext.values


                         sens = Accu2Inst(senNext, senNow, dt.seconds )
                         late = Accu2Inst(latNext, latNow, dt.seconds )
                         sola = Accu2Inst(swrNext, swrNow, dt.seconds )
                         ther = Accu2Inst(lwrNext, lwrNow, dt.seconds )

                         shflux = late + sens + sola + ther

                      elif key == "swflux":

                         #outFieldPrevious = grbPrevious.values
                         #outFieldNext = grbNext.values

                         outField = grb.values
                         evaNow=grb.values
                         evaNext=grbNext.values

                         prec = Accu2Inst(raiNext, raiNow, dt.seconds )
                         evap = Accu2Inst(evaNext, evaNow, dt.seconds )

                         swflux = evap + prec # Evaporation has + sign because positive evaporation means condensation
                                              # and negative evaporation means evaporation

                      else:

                         outField = grb.values
                         #outFieldPrevious = grbPrevious.values
                         outFieldNext = grbNext.values
  
                     
                      # Compute instantanenous fields from accumulated data

                      if key == "shflux":

            #             dt = msgTime - msgTimePrevious
           
                         outFieldNow = shflux
                      
                         aveFieldAccum = aveFieldAccum + outFieldNow
                      
                         if aveStep[j] == 0:

                            aveEndTime = msgTime

                            outField = aveFieldAccum / avePeriod
                            aveTime = aveStartTime + datetime.timedelta(days=0.5*avePeriodDays)


                            outTime = aveTime - refDate

                            outTimeDays = outTime.days + outTime.seconds/sec2Day

                            print("Averages from " + str(aveStartTime) +
                                " to " + str(aveEndTime) + " at " + str(aveTime) +
                                " (" + str(outTimeDays) + " days since " + 
                                str(refDate) + ")")

                            writeASCIIField(outId, outTimeDays, np.flipud(outField))
                            aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here
                      
                      elif key == "swflux":

          #               dt = msgTime - msgTimePrevious
           
                         outFieldNow = swflux
                      
                         aveFieldAccum = aveFieldAccum + outFieldNow
                      
                         if aveStep[j] == 0:

                            aveEndTime = msgTime

                            outField = aveFieldAccum / avePeriod
                            aveTime = aveStartTime + datetime.timedelta(days=0.5*avePeriodDays)


                            outTime = aveTime - refDate

                            outTimeDays = outTime.days + outTime.seconds/sec2Day

                            print("Averages from " + str(aveStartTime) +
                                " to " + str(aveEndTime) + " at " + str(aveTime) +
                                " (" + str(outTimeDays) + " days since " + 
                                str(refDate) + ")")

                            writeASCIIField(outId, outTimeDays, np.flipud(outField))
                            aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here


                      else: #key != "shflux" or key != "swflux": 
                      
           #              dt = msgTime - msgTimePrevious
           
                         outFieldNow = Accu2Inst(outFieldNext, outField, dt.seconds )
                      
                         aveFieldAccum = aveFieldAccum + outFieldNow
                      
                         if aveStep[j] == 0:

                            aveEndTime = msgTime

                            outField = aveFieldAccum / avePeriod
                            aveTime = aveStartTime + datetime.timedelta(days=0.5*avePeriodDays)


                            outTime = aveTime - refDate

                            outTimeDays = outTime.days + outTime.seconds/sec2Day

                            print("Averages from " + str(aveStartTime) +
                                " to " + str(aveEndTime) + " at " + str(aveTime) +
                                " (" + str(outTimeDays) + " days since " + 
                                str(refDate) + ")")

                            writeASCIIField(outId, outTimeDays, np.flipud(outField))
                            
                            aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here
 
               outId.close()

        # Compute derived variables and write to file

    #    key = "Qair" # Air relative humidity at 2 m

    #   outName = outputHeader + key.upper() + "_" + presentCode + outputTail
    #    outId = open(outName,'w') 

    #    writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

    #    aveTimeAccum = 0 # We accumulate the averaging times here
    #    aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

    #    for k in range(len(insTime)):

    #       dew = dewT[k]
    #       air = airT[k]
    #       qair =  surfRelHumidity( d2, t2)/100

    #       aveFieldAccum = aveFieldAccum + qair 

    #       aveTimeAccum = aveTimeAccum + insTime[k]
                      
    #       if aveStep[k] == 0:

    #          outField = aveFieldAccum / avePeriod
    #          aveTime = aveTimeAccum / avePeriod   

    #          writeASCIIField(outId, aveTime, np.flipud(outField))
    #          aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

    #    outId.close()

#        key = "shflux" # Surface net heat flux

#        outName = outputHeader + key.upper() + "_" + presentCode + outputTail
#        outId = open(outName,'w')

#        writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

#        aveTimeAccum = 0 # We accumulate the averaging times here
#        aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

#        for k in range(len(insTime)):

#           sens = Accu2Inst(senPrevious[k], senNext[k], np.array(-1), dt.seconds )
#           late = Accu2Inst(latPrevious[k], latNext[k], np.array(-1), dt.seconds )
#           sola = Accu2Inst(swrPrevious[k], swrNext[k], np.array(-1), dt.seconds )
#           ther = Accu2Inst(lwrPrevious[k], lwrNext[k], np.array(-1), dt.seconds )

#           shflux = sens + late + sola + ther

#           aveFieldAccum = aveFieldAccum + shflux

#           aveTimeAccum = aveTimeAccum + insTime[k]
                      
#           if aveStep[k] == 0:

#              outField = aveFieldAccum / avePeriod
#              aveTime = aveTimeAccum / avePeriod   

#              writeASCIIField(outId, aveTime, np.flipud(outField))
#              aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here
           #writeASCIIField(outId, insTime[k], np.flipud(shflux))
 
#        outId.close()

#        key = "swflux" # Surface net fresh flux

#        outName = outputHeader + key.upper() + "_" + presentCode + outputTail
#        outId = open(outName,'w')

#        writeASCIIHeader(outId, int(nFlds), nLons, nLats, lons[0,:], lats[:,0])

#        aveTimeAccum = 0 # We accumulate the averaging times here
#        aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

#        for k in range(len(insTime)):

#           prec = Accu2Inst(raiPrevious[k], raiNext[k], np.array(-1), dt.seconds )
#           evap = Accu2Inst(evaPrevious[k], evaNext[k], np.array(-1), dt.seconds )

#           swflux = prec + evap # Evaporation has + sign because positive evaporation means condensation
                                # and negative evaporation means evaporation

#           aveFieldAccum = aveFieldAccum + swflux 

#           aveTimeAccum = aveTimeAccum + insTime[k]
                      
#           if aveStep[k] == 0:

#              outField = aveFieldAccum / avePeriod
#              aveTime = aveTimeAccum / avePeriod   

#              writeASCIIField(outId, aveTime, np.flipud(outField))
#              aveFieldAccum = np.zeros(lats.shape) # Accumulate fields here

#        outId.close()
#           writeASCIIField(outId, insTime[k], np.flipud(swflux))



