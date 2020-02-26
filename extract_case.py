import os,netCDF4
import datetime
import read_data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

##########################################################################################
##########################################################################################
def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)

##########################################################################################
# Animation function
##########################################################################################
def animate(i): 
    z = cld_hgt2d[:,:,i]
    cont = plt.contourf(lon, lat, z, levels, cmap='YlGnBu')
    t1 = datetime.datetime(year[i],month[i],day[i],\
                           hour[i],minute[i])
    ax.set(title=t1, ylabel='latitude',xlabel='longitude')
    return cont

levels = np.arange(0, 15000, 1000)

##########################################################################################
# Configuration
##########################################################################################

# Data location (OpenDap)
dirData = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ATOMIC/data/clavrx/2km_01min/'

# What data to read in? [year,month,day,hour,minute]
t_start = [2020,2,9,21,0]
t_stop  = [2020,2,9,23,59]

# Subset the domain? [lon1,lat1,lon2,lat2]
sub_extent = [-60,15,-58,17]
#sub_extent = [] #Using the full domain.

# Name for output movie
fileOUT = 'hot_towers'       + '_'+ \
    str(t_start[0]).zfill(4) + \
    str(t_start[1]).zfill(2) + \
    str(t_start[2]).zfill(2) + \
    str(t_start[3]).zfill(2) + \
    str(t_start[4]).zfill(2) + '_'+ \
    str(t_stop[0]).zfill(4)  + \
    str(t_stop[1]).zfill(2)  + \
    str(t_stop[2]).zfill(2)  + \
    str(t_stop[3]).zfill(2)  + \
    str(t_stop[4]).zfill(2)

##########################################################################################
# Determine problem size...
##########################################################################################
# File name format (file_prefix)YYYY_DOY_HHMM(file_suffix)
file_prefix = 'clavrx_goes16_'
file_suffix = '_BARBADOS-2KM-FD.level2.nc'

# Compute day-of-year (doy), used in filenaming convention.
deltaDays = datetime.date(t_start[0],t_start[1],t_start[2]) - datetime.date(t_start[0],1,1)
doy = deltaDays.days + 1

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
# Read in cloud fields
##########################################################################################
init  = 1
count = 0
for iTime in range(fileiI,fileiF+1):
    [cld_hgt, cld_mask, lon, lat, error] = read_data.read_data(fileList[iTime], sub_extent)
    if not (error):        
        print('      Reading in '+fileList[iTime])

        # Pull out time information from filename
        tempStr = fileList[iTime]
        si      = len(tempStr)-len(file_suffix)-4
        sf      = len(tempStr)-len(file_suffix)
        hour0   = int(tempStr[si:si+2])
        minute0 = int(tempStr[si+2:sf])
        si      = len(dirData)
        year0   = int(tempStr[si:si+4])
        month0  = int(tempStr[si+5:si+7])
        day0    = int(tempStr[si+8:si+10])

        if (init):
            nlon      = len(lon[:,0])
            nlat      = len(lat[0,:])
            nPts      = nlon*nlat
            year      = year0
            month     = month0
            day       = day0
            hour      = hour0
            minute    = minute0
            cld_hgt2d = np.dstack((cld_hgt,cld_hgt))
            init      = 0            
        else:
            year   = np.vstack((year,   year0))
            month  = np.vstack((month,  month0))
            day    = np.vstack((day,    day0))
            hour   = np.vstack((hour,   hour0))
            minute = np.vstack((minute, minute0))
            if (count == 1): cld_hgt2d[:,:,1] = cld_hgt
            if (count > 1):  cld_hgt2d = np.append(cld_hgt2d,  np.atleast_3d(cld_hgt), axis=2)

        # Increment counter
        count = count + 1
    else:
        # If not, squack
        print('      Missing day: '+fileList[iTime]+' does not exist')

# Make plot
fig = plt.figure()
ax  = plt.axes()
anim = animation.FuncAnimation(fig, animate, frames=count-1)
writer = animation.writers['ffmpeg'](fps=1)
anim.save(fileOUT+'.mp4',writer=writer,dpi=100)