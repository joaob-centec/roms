# ---------------------------------------------------------------------
# Download_ERA5_CDS_GRIB.py: script to download ERA5 data from the CDS
# archive at cds.climate.copernicus.eu. Currently, the script downloads
# data that is needed to force the ROMS model. Forcing data is stored
# in monthly data files, one for each variable. The variables that the
# script downloads are the following:
#
# Table 1:
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
# Dwave      mean_wave_direction               140.230  º         I
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
#
#
# Revision list:
#
# Date          Changes                        Author
# ----          -------                        ------
# 2020/11       Created.                       Joao H. Bettencourt
#
#

import calendar
import sys
import cdsapi
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



# c = cdsapi.Client()

def era5_request(c, param, year, month, day, time, area, name_grib):
    """      
        An ERA 5 request for data.
        Change the keywords below to adapt it to your needs.
    """
    c.retrieve(
    'reanalysis-era5-single-levels',
    {
        "product_type": "reanalysis",
        "variable": param,
        "year": year,
        "month": month,
        "day": day,
        "time": time,
        "area": area,
        "format": "grib",
    },
    name_grib)

# ERA5 CDS Variables

Variables = { # Input to data retrieval request
"sustr":"instantaneous_eastward_turbulent_surface_stress",           
"svstr":"instantaneous_northward_turbulent_surface_stress",
"Uair":"10m_u_component_of_wind",         
"Vair":"10m_v_component_of_wind",        
"Awave":"significant_height_of_combined_wind_waves_and_swell",  
"Dwave":"mean_wave_direction",              
"Pwave":"mean_wave_period",                
"Ssize":"N/A",
"Sdens":"N/A",
"Tdew":"2m_dewpoint_temperature",        
"Tair":"2m_temperature",                
"Pair":"surface_pressure",             
"Qair":"N/A",
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

Short_names = { # We use these to build the output file names
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
"Qair":"",
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

# Download parameters

yearStart = args.year0 
yearEnd = args.year1
monthStart = args.month0
monthEnd = args.month1


days31=[] # List with day number as 2 digit string

for i in range(31):

    days31.append("{:02d}".format(i+1)) 

hours24=[]

for i in range(24):

    hours24.append("{:02d}".format(i)+":00") 

# Retrieval area

North = 50.
West = -39.
South = 23.
East = 1.5

Area = [North, West, South, East]
 
# Loop years and months and download the data

C = cdsapi.Client()

for year in list(range(yearStart, yearEnd + 1)):
    for month in list(range(monthStart, monthEnd + 1)):

       for key in Variables:

        if Variables[key] == "N/A":

           print(" *** Skipping Field *** ")

        else:

          numberOfDays = calendar.monthrange(year, month)[1]
        
          Day = days31[:numberOfDays]
        
          Year = "{:4d}".format(year)
          Month = "{:02d}".format(month)

          Name_grib = "ERA5_cds_" + Short_names[key].upper() + "_" + Year + Month + ".grib" 
        
          Param = Variables[key]
          Time  = hours24
        
          print("Processing request: ")
          print("param:"+Param)
          print("year:"+Year)        
          print("month:"+Month)
          print(" from day: "+Day[0]+" to day: "+Day[-1])        

          # Submit request
          era5_request(C, Param, Year, Month, Day, Time, Area, Name_grib)


 
