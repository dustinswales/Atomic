import os,netCDF4
import datetime
import read_data
import numpy as np
import matplotlib.pyplot as plt

##########################################################################################
##########################################################################################
def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)

##########################################################################################
# Configuration
##########################################################################################

# Set true to make plots cloud-fraction when "hot-towers" are detected.
mkplot = 0

# Data location (OpenDap)
dirData = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ATOMIC/data/clavrx/2km_01min/'

# What data to read in? [year,month,day,hour,minute]
t_start = [2020,1,14,15,40]
#t_stop  = [2020,1,14,15,45]
t_stop  = [2020,2,13,23,59]

# Subset the domain? [lon1,lat1,lon2,lat2]
sub_extent = [-60,15,-58,17]
#sub_extent = [] #Using the full domain.

# Here "hot-towers" are defined by a sudden change in the number of pixels with cloud-top
# height exceeding some threshold.
# Cloud height threshold
cld_hgt_thresh   = 4000.
# Relative change in pixels containing high-clouds.
changeInPixCount = 0.5

fileOUT = 'hot_towers'       + '_' + \
    str(t_start[0]).zfill(4) + \
    str(t_start[1]).zfill(2) + \
    str(t_start[2]).zfill(2) + \
    str(t_start[3]).zfill(2) + \
    str(t_start[4]).zfill(2) + '_' + \
    str(t_stop[0]).zfill(4)  + \
    str(t_stop[1]).zfill(2)  + \
    str(t_stop[2]).zfill(2)  + \
    str(t_stop[3]).zfill(2)  + \
    str(t_stop[4]).zfill(2)  + '_' + \
    str(cld_hgt_thresh/1000.).zfill(2)+'km_dp'+str(changeInPixCount*100.).zfill(2)+'.nc'

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
    [cld_hgt, cld_temp, cld_pres, lw_clrCh, lon, lat, error] = \
        read_data.read_data(fileList[iTime], sub_extent)
    if not (error):        
        print('      Reading in '+fileList[iTime])
        
        # What fraction of the domain is covered in high clouds?
        mask1 = np.where(cld_hgt >= cld_hgt_thresh,1,0)
        cfp0  = np.sum(mask1)
        
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
            nlon   = len(lon[:,0])
            nlat   = len(lat[0,:])
            nPts   = nlon*nlat
            year   = year0
            month  = month0
            day    = day0
            hour   = hour0
            minute = minute0
            cfp    = cfp0
            init   = 0
        else:
            year   = np.vstack((year,   year0))
            month  = np.vstack((month,  month0))
            day    = np.vstack((day,    day0))
            hour   = np.vstack((hour,   hour0))
            minute = np.vstack((minute, minute0))
            cfp    = np.vstack((cfp,    cfp0))
            cld_hgt2d  = np.dstack((cld_hgtm1t,cld_hgt))
            cld_temp2d = np.dstack((cld_tempm1t,cld_temp))
            cld_pres2d = np.dstack((cld_presm1t,cld_pres))
        cld_hgtm1t  = cld_hgt
        cld_tempm1t = cld_temp
        cld_presm1t = cld_pres

        # If number of pixels jumps drastically, say 20%?, from subsequent timesteps, make plot
        if (count > 0):
            t1 = datetime.datetime(year[count-1],month[count-1],day[count-1],\
                                   hour[count-1],minute[count-1])
            t2 = datetime.datetime(year[count],  month[count],  day[count],\
                                   hour[count],  minute[count])
            if (cfp[count] > 0):
                dcfp = float(cfp0-cfp[count-1])/float(cfp[count])
            else:
                dcfp = 0
            if (dcfp > changeInPixCount):
                print("   Possible hot-tower identified @ ",t2 )
                print("      Pixels exceeding threshold increased by "+str(dcfp*100.)+\
                      "% , from "+str(cfp[count-1])+" to "+str(cfp[count])+" pixels")
                if (mkplot):
                    # Cloud-top height
                    levels = np.arange(0, 15, 1)
                    fig = plt.figure(figsize=(11, 7))                
                    ax1 = fig.add_subplot(2,3,1)
                    s1 = ax1.contourf(lon, lat, cld_hgt2d[:,:,0]/1000., levels, cmap='YlGnBu')
                    ax1.set(title=t1, ylabel='latitude')
                    clb1 = fig.colorbar(s1, ax=ax1, shrink=0.9)
                    clb1.set_label('(km)') 
                    ax2 = fig.add_subplot(2,3,4)        
                    s2 = ax2.contourf(lon, lat, cld_hgt2d[:,:,1]/1000., levels, cmap='YlGnBu')
                    ax2.set(title=t2, ylabel='latitude',xlabel='longitude')
                    clb2 = fig.colorbar(s2, ax=ax2, shrink=0.9)
                    clb2.set_label('(km)')
                    
                    # Cloud-top temperature
                    levelsT = np.arange(200, 300, 2)
                    ax3 = fig.add_subplot(2,3,2)        
                    s3 = ax3.contourf(lon, lat, cld_temp2d[:,:,0], levelsT, cmap='YlGnBu')
                    ax3.set(title=t1, yticklabels=[])
                    clb3 = fig.colorbar(s3, ax=ax3, shrink=0.9)
                    clb3.set_label('(K)')
                    ax4 = fig.add_subplot(2,3,5)        
                    s4 = ax4.contourf(lon, lat, cld_temp2d[:,:,1], levelsT, cmap='YlGnBu')
                    ax4.set(title=t2, xlabel='longitude', yticklabels=[])
                    clb4 = fig.colorbar(s4, ax=ax4, shrink=0.9)
                    clb4.set_label('(K)')

                    # Cloud-top pressure
                    levelsP = np.arange(100, 1100, 50)
                    ax5 = fig.add_subplot(2,3,3)        
                    s5 = ax5.contourf(lon, lat, cld_pres2d[:,:,0], levelsP, cmap='YlGnBu')
                    ax5.set(title=t1, yticklabels=[])
                    clb5 = fig.colorbar(s5, ax=ax5, shrink=0.9)
                    clb5.set_label('(hPa)')
                    ax6 = fig.add_subplot(2,3,6)        
                    s6 = ax6.contourf(lon, lat, cld_pres2d[:,:,1], levelsP, cmap='YlGnBu')
                    ax6.set(title=t2, xlabel='longitude', yticklabels=[])
                    clb6 = fig.colorbar(s6, ax=ax6, shrink=0.9)
                    clb6.set_label('(hPa)')
                    
                    #
                    plt.show()                    

        # Increment counter
        count = count + 1
    else:
        # If not, squack
        print('      Missing day: '+fileList[iTime]+' does not exist')

##########################################################################################
# Create output file
##########################################################################################
ncfile = netCDF4.Dataset(fileOUT,'w')
ncfile.createDimension('nTime',count)
yearOUT   = ncfile.createVariable('year',       np.dtype('int32').char, ('nTime'))
monthOUT  = ncfile.createVariable('month',      np.dtype('int32').char, ('nTime'))
dayOUT    = ncfile.createVariable('day',        np.dtype('int32').char, ('nTime'))
hourOUT   = ncfile.createVariable('hour',       np.dtype('int32').char, ('nTime'))
minuteOUT = ncfile.createVariable('minute',     np.dtype('int32').char, ('nTime'))
npixelOUT = ncfile.createVariable('nPixel',     np.dtype('int32').char, ('nTime'))
ntotpOUT  = ncfile.createVariable('nTotalPixel',np.dtype('int32').char, ('nTime'))
yearOUT[:]   = year
monthOUT[:]  = month
dayOUT[:]    = day
hourOUT[:]   = hour
minuteOUT[:] = minute
npixelOUT[:] = cfp
ntotpOUT[:]  = nPts
ncfile.close()
