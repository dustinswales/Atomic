##########################################################################################
#
##########################################################################################
import os,netCDF4,numpy as np
import datetime

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)
        
##########################################################################################
# Configuration
##########################################################################################
debug     = 1
storeData = 0

# Data location
#dirData = '/Projects/ATOMIC/data/clavrx/2km_01min/'
dirData = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ATOMIC/data/clavrx/2km_01min/'

# What data to read in? [year,month,day,hour,minute]
t_start = [2020,1,18,13,55]
t_stop  = [2020,1,22,0,3]

# Subset the domain? [lon1,lat1,lon2,lat2]
sub_extent = [-59,14,-55,17]

# File name format (file_prefix)YYYY_DOY_HHMM(file_suffix)
file_prefix = 'clavrx_goes16_'
file_suffix = '_BARBADOS-2KM-FD.level2.nc'

##########################################################################################
# Determine problem size...
##########################################################################################
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

# Are we only reading in a portion of the domain?
if (len(sub_extent) == 4):
    print('Using subdomain')
    subDomain = 1
else:
    print('Using full domain')
    subDomain = 0

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
dirI   = dirData + str(t_start[0])+'_'+str(t_start[1]).zfill(2)+'_'+str(t_start[2]).zfill(2)+'_'+str(doy).zfill(3)+'/'
fnameI = dirI + file_prefix+str(t_start[0])+'_'+str(doy).zfill(3)+'_'+ \
    str(t_start[3]).zfill(2)+str(t_start[4]).zfill(2)+file_suffix

# What is the last file to read in?
dirF   = dirData + str(t_stop[0])+'_'+str(t_stop[1]).zfill(2)+'_'+str(t_stop[2]).zfill(2)+'_'+str(doy+countDAY-1).zfill(3)+'/'
fnameF = dirF + file_prefix+str(t_stop[0])+'_'+str(doy+countDAY-1).zfill(3)+'_'+ \
    str(t_stop[3]).zfill(2)+str(t_stop[4]).zfill(2)+file_suffix

# Pull out indices for first/last files.
fileiI = fileList.index(fnameI)
fileiF = fileList.index(fnameF)

##########################################################################################
# Loop over all files and read in data
##########################################################################################
init        = 1
countRecord = 0
for iFile in range(fileiI,fileiF+1):

    # Only read in file if it exists.
    try:
        # Open file
        dataIN = netCDF4.Dataset(fileList[iFile],'r')
        time   = dataIN.variables['scan_line_time'][:]

        # Initialize
        if (init):
            lat    = dataIN.variables['latitude'][:,:]
            lon    = dataIN.variables['longitude'][:,:]
            nlon   = len(lon[:,0])
            nlat   = len(lat[0,:])
            # Are we only reading in a portion of the full domain?
            if (subDomain):
                mask1 = np.where(lon >= sub_extent[0],1,0)
                mask2 = np.where(lon <= sub_extent[2],1,0)
                mask3 = np.where(lat >= sub_extent[1],1,0)
                mask4 = np.where(lat <= sub_extent[3],1,0)
                mask  = mask1*mask2*mask3*mask4
                # Since the data is not on a regular grid, we will use a grid that covers the
                # entire domain, but with some extra bits on the end.
                Xextent = np.where(np.sum(mask,axis=0) > 0)
                Yextent = np.where(np.sum(mask,axis=1) > 0)
                xi0     = np.min(Xextent)
                xi1     = np.max(Xextent)
                yi0     = np.min(Yextent)
                yi1     = np.max(Yextent)
                # Reset lon/lat arrays
                lon  = lon[xi0:xi1,yi0:yi1]
                lat  = lat[xi0:xi1,yi0:yi1]
                nlon = len(lon[:,0])
                nlat = len(lat[0,:])
            
        # Pull out hour/minute information from filename
        tempStr = fileList[iFile]
        si      = len(tempStr)-len(file_suffix)-4
        sf      = len(tempStr)-len(file_suffix)
        hour0   = int(tempStr[si:si+2])
        minute0 = int(tempStr[si+2:sf])
                
        # Read in fields
        if (debug): print('      Reading in '+fileList[iFile])
        var       = dataIN.variables['cloud_mask']
        ao        = var.getncattr('add_offset')
        sf        = var.getncattr('scale_factor')
        cld_mask0 = var[xi0:xi1,yi0:yi1]*sf + ao        
        var       = dataIN.variables['cld_height_acha']
        ao        = var.getncattr('add_offset')
        sf        = var.getncattr('scale_factor')
        cld_hgt0  = var[xi0:xi1,yi0:yi1]*sf + ao
        
        # Store data (Append along record dimension)
        if (init):
            year     = single_date.year
            month    = single_date.month
            day      = single_date.day
            hour     = hour0
            minute   = minute0
            init     = 0
        if (countRecord > 0):
            year   = np.vstack((year,   single_date.year))
            month  = np.vstack((month,  single_date.month))
            day    = np.vstack((day,    single_date.day))
            hour   = np.vstack((hour,   hour0))
            minute = np.vstack((minute, minute0))
        if (storeData):
            if (countRecord == 0):
                cld_mask = np.dstack((cld_mask0,cld_mask0))
                cld_hgt  = np.dstack((cld_hgt0, cld_hgt0))
            if (countRecord == 1):
                cld_mask[:,:,countRecord] = cld_mask0
                cld_hgt[:,:,countRecord]  = cld_hgt0
            if (countRecord > 1):
                cld_mask = np.append(cld_mask, np.atleast_3d(cld_mask0), axis=2)
                cld_hgt  = np.append(cld_hgt,  np.atleast_3d(cld_hgt0), axis=2)

        # If not storing the data, do processing here...
        if not(storeData):
            if (countRecord==0): print('Not saving data, doing some analysis and tossing')
                                        
        # Increment counter
        countRecord = countRecord + 1
         
    except:
        # If not, squack
        print('Missing day: '+fileList[iFile]+' does not exist')
        
#
if (countRecord < max_tsteps and countRecord > 0):
    print('######################################################################')
    print('Warning: Some of the requested data is missing...')
    print('   Maximum time-steps avaibale for requested period:  ' + str(max_tsteps))
    print('   Number of time-steps ingested:                     ' + str(countRecord))
    print('######################################################################')
    
##########################################################################################
# END PROGRAM
##########################################################################################
