##########################################################################################
#
# About: This script imports CLAVRx (Cloud from AVHRR extended processing system) from
#        SSEC to NOAA PSD.
#
##########################################################################################
from ftplib import FTP
import os
from datetime import timedelta, date
import os.path
from os import path

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days+1)):
        yield start_date + timedelta(n)
        
##########################################################################################
# Import configuration
##########################################################################################
# What is the name of the FTP server?
ftpServer = 'ftp.ssec.wisc.edu'

# Where is the data stored on the FTP server?
dirFTP = '/pub/clavrx/barbados/'

# Which configuration to download?
configID = '2km_01min'

# Where should the data be stored locally?
dirLocal = '/Projects/ATOMIC/data/clavrx/'

# What date range to import?
start_date = date(2020, 1, 24)
end_date   = date(2020, 1, 26)

# Compute day-of-year, used in filenaming convention.
deltaDays = start_date - date(start_date.year,1,1)
countDOY  = deltaDays.days + 1

##########################################################################################
# Import data
##########################################################################################
# Make local root directory (only if it doesn't exist)
try:
    os.mkdir(dirLocal+configID)
except:
    print("Root directory already exists")

# Loop over all times and pull over data, when available
for single_date in daterange(start_date, end_date):
    # On the FTP site the data is archived in directories with the following format.
    NewTempDir = single_date.strftime("%Y_%m_%d")+'_'+str(countDOY).zfill(3)
    countDOY   = countDOY+1

    # Login to FTP server
    ftp = FTP(ftpServer)
    ftp.login()
    # Is the data we requested on the server?
    try:
        ftp.cwd(dirFTP+configID+'/'+NewTempDir+'/')
        if not (os.path.exists(dirLocal+configID+'/'+NewTempDir+'/')):
            os.mkdir(dirLocal+configID+'/'+NewTempDir+'/')
        os.chdir(dirLocal+configID+'/'+NewTempDir+'/')
        ls = ftp.nlst()
        count = len(ls)
        curr = 0
        print "found {} files".format(count)
        for fn in ls:
            curr += 1
            print 'Processing file {} ... {} of {} ...'.format(fn, curr, count)
            ftp.retrbinary('RETR ' + fn, open(fn, 'wb').write)
    except:
        print("Data missing on FTP site: "+NewTempDir)
        
ftp.quit()
print "download complete."

    
##########################################################################################
# END PROGRAM
##########################################################################################
