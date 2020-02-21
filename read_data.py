##########################################################################################
#
##########################################################################################
from matplotlib import pyplot as plt
import os,fnmatch,netCDF4,numpy as np
import datetime
from ftplib import FTP

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)
##########################################################################################
# Configuration
##########################################################################################
# Data location
dirData = "/Projects/ATOMIC/data/clavrx/2km_01min/"

# What data to read in? [yyyy,mm,dd,hh,mm]
t_start = [2020,01,31,23,20]
t_stop  = [2020,02,01,01,40]

##########################################################################################
#
##########################################################################################
# Compute day-of-year (doy), used in filenaming convention.
deltaDays = datetime.date(t_start[0],t_start[1],t_start[2]) - datetime.date(t_start[0],1,1)
doy = deltaDays.days + 1

# How many timesteps were requested?
to = datetime.datetime(t_start[0], t_start[1], t_start[2], t_start[3], t_start[4])
tf = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2],  t_stop[3],  t_stop[4])
dt = tf - to
max_tsteps = dt.days*24*60 + dt.seconds/60
print("Maximum number of time-steps (dt=1min): "+str(max_tsteps))

# Data is stored in daily directories. How many need to be parsed?
to1 = datetime.datetime(t_start[0], t_start[1], t_start[2])
tf1 = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2])
nDirsToReadFrom = abs(tf1-to1).days+1
print("Requested time period spans "+str(nDirsToReadFrom)+" days")

##########################################################################################
# Loop over daily directories and read in data
##########################################################################################
print("Looping over all possible days to read in...")
countDAY = 0
for single_date in daterange(to1,tf1):
    # What directory is the data stored in?
    dirDaily = single_date.strftime("%Y_%m_%d")+'_'+str(doy+countDAY).zfill(3)
    
    # Does this directory exist?
    try:
        # If so, then read in the data....
        fileList = sorted(os.listdir(dirData+dirDaily))
        nFiles   = len(fileList)
        dir = dirData+dirDaily
        print("   "+dir+" contains "+str(nFiles)+" files.")

        # Does the directory contain the requested starting time?
        if (countDAY == 0):
            fname = 'clavrx_goes16_'+str(t_start[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_start[3]).zfill(2)+str(t_start[4]).zfill(2)+'_BARBADOS-2KM-FD.level2.nc'
            try:        
                t_begin = fileList.index(fname)
                t_end   = nFiles
            except:
                print("Requested starting time not part of dataset")
        # Starting indices for subsequent days are always t=0, except in final timestep.
        if (countDAY > 0):
            t_begin = 0
            t_end   = nFiles
        # Final timestep.
        if (countDAY == nDirsToReadFrom-1):            
            fname = 'clavrx_goes16_'+str(t_stop[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_stop[3]).zfill(2)+str(t_stop[4]).zfill(2)+'_BARBADOS-2KM-FD.level2.nc'
            try:        
                t_begin = 0
                t_end   = fileList.index(fname)+1
            except:
                print("Requested ending time not part of dataset")
                              
        
        # Loop over files in directory, read in
        for ij in range(t_begin,t_end):
            print("      Reading in "+dir+"/"+fileList[ij])
            dataIN          = netCDF4.Dataset(dir+"/"+fileList[ij],'r')
            time            = dataIN.variables['scan_line_time'][:]
            lat             = dataIN.variables['latitude'][:,:]
            lon             = dataIN.variables['longitude'][:,:]
            cloud_mask      = dataIN.variables['cloud_mask'][:,:]
            cld_height_acha = dataIN.variables['cld_height_acha']
            ao              = cld_height_acha.getncattr('add_offset')
            sf              = cld_height_acha.getncattr('scale_factor')
            cld_height_acha = cld_height_acha*sf+ao
    except:
        # If not, squack
        print("Missing day: "+dirDaily+" does not exist")

    # Increment daily counter
    countDAY = countDAY + 1
    
##########################################################################################
# END PROGRAM
##########################################################################################

