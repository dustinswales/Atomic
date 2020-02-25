##########################################################################################
#
##########################################################################################
import netCDF4,numpy as np

def read_data(fileIN, sub_extent):

    # Are we only reading in a portion of the domain?
    if (len(sub_extent) == 4):
        subDomain = 1
    else:
        subDomain = 0

    # Only read in file if it exists. 
    try:
        # Open file
        dataIN = netCDF4.Dataset(fileIN,'r')
        time   = dataIN.variables['scan_line_time'][:]
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
            nlon = len(lon[:,0])
            nlat = len(lat[0,:])
        else:
            xi0 = 0
            xi1 = nlon-4
            yi0 = 0
            yi1 = nlat-4
        lon  = lon[xi0:xi1,yi0:yi1]
        lat  = lat[xi0:xi1,yi0:yi1]
            
        # Read in fields
        var      = dataIN.variables['cloud_mask']
        cld_mask = dataIN.variables['cloud_mask'][xi0:xi1,yi0:yi1]        
        var      = dataIN.variables['cld_height_acha']
        cld_hgt  = dataIN.variables['cld_height_acha'][xi0:xi1,yi0:yi1]
        error    = 0
    except:
        error    = 1
    
    return cld_hgt,cld_mask,lon,lat,error
        
        
# Check to see if this file is being executed as the "Main" python
# script instead of being used as a module by some other python script
# This allows us to use the module which ever way we want.
if __name__ == '__main__':
    read_data()

##########################################################################################
# END PROGRAM
##########################################################################################
