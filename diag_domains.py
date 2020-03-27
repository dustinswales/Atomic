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
dirData   = '/Projects/ATOMIC/preliminary/data/clavrx/2km_01min/'
dirData10 = '/Projects/ATOMIC/preliminary/data/clavrx/2km_10min/'

# What data to read in? [year,month,day,hour,minute]
t_start = [2020,1,14,12,00]
t_stop  = [2020,2,16,12,00]

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
  
# Data is stored in daily directories. How many need to be parsed?
to1 = datetime.datetime(t_start[0], t_start[1], t_start[2])
tf1 = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2])
nDirsToReadFrom = abs(tf1-to1).days+1
    
##########################################################################################
# Create file list of 01-minute data. 
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
nTime  = fileiF - fileiI + 1

##########################################################################################
# Try to read in field from file, if no file, look for 10min data.
##########################################################################################
# Create output file
fileIDout = netCDF4.Dataset('clavrx_01min_10min_domain_status.nc', 'w', format='NETCDF4')
# Dimensions
fileIDout.createDimension('time', None)
# Varaibles
yearout = fileIDout.createVariable('year', 'int', ('time',))
yearout.long_name = 'Year'
monthout = fileIDout.createVariable('month', 'int', ('time',))
monthout.long_name = 'Month'
dayout = fileIDout.createVariable('day', 'int', ('time',))
dayout.long_name = 'Day'
hourout = fileIDout.createVariable('hour', 'int', ('time',))
hourout.long_name = 'Hour'
minuteout = fileIDout.createVariable('minute', 'int', ('time',))
minuteout.long_name = 'Minute'
flag01out = fileIDout.createVariable('flag01', 'int', ('time',))
flag01out.long_name = '01-minute data present?'
flag02out = fileIDout.createVariable('flag02', 'int', ('time',))
flag02out.long_name = '10-minute data present?'
flag03out = fileIDout.createVariable('flag03', 'int', ('time',))
flag03out.long_name = '01-minute data subsets 10-minute data?'

# Create flags for domain status.
# Loop over all requested times...
count = 0
for iTime in range(fileiI,fileiF+1):
    # Compute day-of-year (doy), used in filenaming convention.
    deltaDays = datetime.date(yy[iTime],mm[iTime],dd[iTime]) - datetime.date(yy[iTime],1,1)
    doy       = deltaDays.days + 1
    
    # What directory is the 10-min data stored in?
    dirDaily = str(yy[iTime])+"_"+ str(mm[iTime]).zfill(2)+"_"+ str(dd[iTime]).zfill(2)+'_'+str(doy).zfill(3)
    dir      = dirData10+dirDaily
    fileRoot  = dir + '/' + file10_prefix + str(yy[iTime]) + str(doy).zfill(3) + \
        str(hh[iTime]).zfill(2) + str(10*int(np.floor(mmm[iTime]/10))).zfill(2) + file10_suffix

    # Output time information
    yearout[count]   = yy[iTime]
    monthout[count]  = mm[iTime]
    dayout[count]    = dd[iTime]
    hourout[count]   = hh[iTime]
    minuteout[count] = mmm[iTime] 
    
    # Read in 01-minute geo data
    if (os.path.exists(fileList[iTime])):
        dataIN = netCDF4.Dataset(fileList[iTime],'r')
        lat01  = dataIN.variables['latitude'][:,:]
        lon01  = dataIN.variables['longitude'][:,:]
        nlon01 = len(lon01[:,0])
        nlat01 = len(lat01[0,:])
        flag01out[count] = 1
    else:
        flag01out[count] = 0
        
    # Read in 10-minute geo data
    if (os.path.exists(fileRoot)):
        dataIN2  = netCDF4.Dataset(fileRoot,'r')
        lat10    = dataIN2.variables['latitude'][:,:]
        lon10    = dataIN2.variables['longitude'][:,:]
        nlon10   = len(lon10[:,0])
        nlat10   = len(lat10[0,:])
        flag02out[count] = 1
    else:
        flag02out[count] = 0

    # Does 01-minute degree data subset the 10-minute?
    if (flag01out[count] and flag02out[count]):
        if (np.min(lat01) > np.min(lat10) and \
            np.max(lat01) < np.max(lat10) and \
            np.min(lon01) > np.min(lon10) and \
            np.max(lon01) < np.max(lon10)):
            flag03out[count] = 1
        else:
            flag03out[count] = 0
    else:
        flag03out[count] = -1
            
    # Increment counter
    count = count + 1

# Close file
fileIDout.close()
