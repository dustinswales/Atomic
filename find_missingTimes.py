import os,netCDF4
import datetime
import read_data
import numpy as np

##########################################################################################
##########################################################################################
def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)

##########################################################################################
# Configuration
##########################################################################################
# Data location (OpenDap)
dirData   = '/Projects/ATOMIC/data/clavrx/2km_01min/'
dirData10 = '/Projects/ATOMIC/data/clavrx/2km_10min/'

# What data to read in? [year,month,day,hour,minute]
t_start = [2020,1,14,15,52]
t_stop  = [2020,1,14,18,54]

# File name format (file_prefix)YYYY_DOY_HHMM(file_suffix) (01-min data)
file_prefix = 'clavrx_goes16_'
file_suffix = '_BARBADOS-2KM-FD.level2.nc'

# File name format (file10_prefix)YYYYDOYHHMMXXX(file10_suffix) (10-min data)
file10_prefix = 'clavrx_OR_ABI-L1b-RadF-M6C01_G16_s'
file10_suffix = '_BARBADOS-2KM-FD.level2.nc'

##########################################################################################
# Determine problem size...
##########################################################################################
# Compute day-of-year (doy), used in filenaming convention.
deltaDays = datetime.date(t_start[0],t_start[1],t_start[2]) - datetime.date(t_start[0],1,1)
doy       = deltaDays.days + 1

# How many timesteps (1-minute data) were requested?
to = datetime.datetime(t_start[0], t_start[1], t_start[2], t_start[3], t_start[4])
tf = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2],  t_stop[3],  t_stop[4])
dt = tf - to
max_tsteps = dt.days*24*60 + dt.seconds/60
print('######################################################################')
print('Maximum number of time-steps (dt=1min): '+str(max_tsteps))
    
# Data is stored in daily directories. How many need to be parsed?
to1 = datetime.datetime(t_start[0], t_start[1], t_start[2])
tf1 = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2])
nDirsToReadFrom = abs(tf1-to1).days+1
print('Requested time period spans '+str(nDirsToReadFrom)+' days')
    
##########################################################################################
# Create file list
##########################################################################################
fileList = ['' for x in range(nDirsToReadFrom*24*60)]
yy = np.zeros(nDirsToReadFrom*24*60,dtype=int)
mm = np.zeros(nDirsToReadFrom*24*60,dtype=int)
dd = np.zeros(nDirsToReadFrom*24*60,dtype=int)
hh = np.zeros(nDirsToReadFrom*24*60,dtype=int)
mmm = np.zeros(nDirsToReadFrom*24*60,dtype=int)

# Loop over all days spanning data request...
countDAY    = 0
countRecord = 0
for single_date in daterange(to1,tf1):
    # What directory is the data stored in?
    dirDaily = single_date.strftime('%Y_%m_%d')+'_'+str(doy+countDAY).zfill(3)
    dir      = dirData+dirDaily
    
    # Create file names for that day
    for iHour in range(0,24):
        for iMinute in range(0,60):
            fileList[countRecord] = dir + '/' + file_prefix+str(t_start[0])+'_'+ \
                str(doy+countDAY).zfill(3)+'_'+ str(iHour).zfill(2)+str(iMinute).zfill(2)+file_suffix
            yy[countRecord] = single_date.year
            mm[countRecord] = single_date.month
            dd[countRecord] = single_date.day
            hh[countRecord] = iHour
            mmm[countRecord] = iMinute
            countRecord = countRecord + 1
            
    # Increment daily counter
    countDAY = countDAY + 1

# What is the first file to read in?
dirI   = dirData + str(t_start[0])+'_'+str(t_start[1]).zfill(2)+'_'+\
    str(t_start[2]).zfill(2)+'_'+str(doy).zfill(3)+'/'
fnameI = dirI + file_prefix+str(t_start[0])+'_'+str(doy).zfill(3)+'_'+ \
    str(t_start[3]).zfill(2)+str(t_start[4]).zfill(2)+file_suffix

# What is the last file to read in?
dirF   = dirData + str(t_stop[0])+'_'+str(t_stop[1]).zfill(2)+'_'+\
    str(t_stop[2]).zfill(2)+'_'+str(doy+countDAY-1).zfill(3)+'/'
fnameF = dirF + file_prefix+str(t_stop[0])+'_'+str(doy+countDAY-1).zfill(3)+'_'+ \
    str(t_stop[3]).zfill(2)+str(t_stop[4]).zfill(2)+file_suffix

# Pull out indices for first/last files.
fileiI = fileList.index(fnameI)
fileiF = fileList.index(fnameF)

##########################################################################################
# Try to read in field from file, if no file, look for 10min data.
##########################################################################################
init  = 1
count = 0
for iTime in range(fileiI,fileiF+1):
    try:
        # See if file exists                                                                                                                                                                    
        dataIN = netCDF4.Dataset(fileList[iTime],'r')
    except:
        # What directory is the 10-min data stored in?                                                                                                                                          
        dirDaily = single_date.strftime('%Y_%m_%d')+'_'+str(doy).zfill(3)
        dir      = dirData10+dirDaily

        # If file doesn't exist, see if 10min data is available, if thats available copy over
        # Compute day-of-year (doy), used in filenaming convention.
        deltaDays = datetime.date(yy[iTime],mm[iTime],dd[iTime]) - datetime.date(yy[iTime],1,1)
        doy       = deltaDays.days + 1
        fileRoot  = dir + '/' + file10_prefix + str(yy[iTime]) + str(doy).zfill(3) + \
            str(hh[iTime]).zfill(2) + str(10*int(np.floor(mmm[iTime]/10))).zfill(2) + file10_suffix
        print('##########################################################################')
        print(fileList[iTime]+' is missing')
        print(fileRoot+' in its place')
        try:
            dataIN2 = netCDF4.Dataset(fileRoot,'r')

        except:
            print(fileRoot+' is missing as well')
        count = count + 1
