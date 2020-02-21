##########################################################################################
#
##########################################################################################
#from matplotlib import pyplot as plt
import os,netCDF4,numpy as np
import datetime
from ftplib import FTP

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + datetime.timedelta(n)
        
##########################################################################################
# Configuration
##########################################################################################
debug     = 0
storeData = 0

# Data location
dirData = "/Projects/ATOMIC/data/clavrx/2km_01min/"

# What data to read in? [yyyy,mm,dd,hh,mm]
t_start = [2020,02,03,11,50]
t_stop  = [2020,02,03,14,52]

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
print("######################################################################")
print("Maximum number of time-steps (dt=1min): "+str(max_tsteps))

# Data is stored in daily directories. How many need to be parsed?
to1 = datetime.datetime(t_start[0], t_start[1], t_start[2])
tf1 = datetime.datetime(t_stop[0],  t_stop[1],  t_stop[2])
nDirsToReadFrom = abs(tf1-to1).days+1
print("Requested time period spans "+str(nDirsToReadFrom)+" days")

# Are we only reading in a portion of the domain?
if (len(sub_extent) == 4):
    print("Using subdomain")
    subDomain = 1
else:
    print("Using full domain")
    subDomain = 0

##########################################################################################
# Loop over daily directories and read in data
##########################################################################################
countDAY    = 0
countRecord = 0
init        = 1
print("Looping over all possible days to read in...")
for single_date in daterange(to1,tf1):
    # What directory is the data stored in?
    dirDaily = single_date.strftime("%Y_%m_%d")+'_'+str(doy+countDAY).zfill(3)
    
    # Does this directory exist?
    if (os.path.exists(dirData+dirDaily)): 
        # If so, then read in the data....
        fileList = sorted(os.listdir(dirData+dirDaily))
        nFiles   = len(fileList)
        dir      = dirData+dirDaily
        print("Reading in data...")
        
        ##################################################################################
        # Determine starting/ending poitns in daily directory.
        # In the cases where the requested times do not cover an entire day, we only want
        # read in a subset of the data.
        ##################################################################################
        # A) Does the directory contain the requested starting time?
        if (countDAY == 0):
            fname = file_prefix+str(t_start[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_start[3]).zfill(2)+str(t_start[4]).zfill(2)+file_suffix
            try:        
                t_begin = fileList.index(fname)
                t_end   = nFiles
            except:
                print("Requested starting time not present. Reading entire day...")
                t_begin = 0
                t_end   = nFiles
                
        # B) Starting indices for subsequent days are always t=0, except in final timestep.
        if (countDAY > 0):
            t_begin = 0
            t_end   = nFiles
            
        # C) Final timestep.
        if (countDAY == nDirsToReadFrom-1):            
            fname = file_prefix+str(t_stop[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_stop[3]).zfill(2)+str(t_stop[4]).zfill(2)+file_suffix
            try:        
                t_begin = 0
                t_end   = fileList.index(fname)+1
            except:
                print("Requested ending time not present. Reading entire day...")
                t_begin = 0
                t_end   = nFiles
                
        # D) Requested time period is subset of single-day
        if (countDAY == 0 and nDirsToReadFrom == 1):
            fname1 = file_prefix+str(t_start[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_start[3]).zfill(2)+str(t_start[4]).zfill(2)+file_suffix
            fname2 = file_prefix+str(t_stop[0])+'_'+str(doy+countDAY).zfill(3)+'_'+ \
                str(t_stop[3]).zfill(2)+str(t_stop[4]).zfill(2)+file_suffix
            if     (os.path.exists(dir+'/'+fname1)): t_begin = fileList.index(fname1)
            if     (os.path.exists(dir+'/'+fname2)): t_end   = fileList.index(fname2)
            if not (os.path.exists(dir+'/'+fname1)): t_begin = 0
            if not (os.path.exists(dir+'/'+fname2)): t_end   = nFiles
                
        ##################################################################################
        # Loop over files in directory, read in requested fields, store data.
        ##################################################################################
        for ij in range(t_begin,t_end):
            # Open file
            dataIN = netCDF4.Dataset(dir+"/"+fileList[ij],'r')
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
            tempStr = fileList[ij]
            si      = len(tempStr)-len(file_suffix)-4
            sf      = len(tempStr)-len(file_suffix)
            hour0   = int(tempStr[si:si+2])
            minute0 = int(tempStr[si+2:sf])
 
                
            # Read in fields
            if (debug): print("      Reading in "+dir+"/"+fileList[ij])
            var       = dataIN.variables["cloud_mask"]
            ao        = var.getncattr('add_offset')
            sf        = var.getncattr('scale_factor')
            cld_mask0 = var[xi0:xi1,yi0:yi1]*sf + ao        
            var       = dataIN.variables["cld_height_acha"]
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
                if (countRecord==0): print("Not saving data, doing some analysis and tossing")
                                        
            # Increment counter
            countRecord = countRecord + 1
                
    else:
        # If not, squack
        print("Missing day: "+dirDaily+" does not exist")

    # Increment daily counter
    countDAY = countDAY + 1
      
#
if (countRecord < max_tsteps and countRecord > 0):
    print("######################################################################")
    print("Warning: Some of the requested data is missing...")
    print("   Maximum time-steps avaibale for requested period:  " + str(max_tsteps))
    print("   Number of time-steps ingested:                     " + str(countRecord))
    for ij in range(0,countRecord):
        print("      "+str(datetime.datetime(year[ij],month[ij],day[ij],hour[ij],minute[ij])))

    print("######################################################################")
    
##########################################################################################
# END PROGRAM
##########################################################################################

