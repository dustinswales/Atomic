#########################################################################################
#
##########################################################################################
from ftplib import FTP
import os
from datetime import timedelta, date

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

start_date = date(2020, 1, 14)
end_date = date(2020, 2, 14)
count = 14
for single_date in daterange(start_date, end_date):
    localDir=single_date.strftime("%Y_%m_%d")+'_'+str(count).zfill(3)
    count = count+1

    ftp = FTP('ftp.ssec.wisc.edu')
    print "Welcome: ", ftp.getwelcome()
    ftp.login()
    ftp.cwd('/pub/clavrx/barbados/2km_01min/'+localDir+'/')
    ftp.retrlines('LIST')
    os.chdir('/data/dswales/Atomic/'+localDir+'/')
    ls = ftp.nlst()
    count = len(ls)
    curr = 0
    print "found {} files".format(count)
    for fn in ls:
        curr += 1
        print 'Processing file {} ... {} of {} ...'.format(fn, curr, count)
        ftp.retrbinary('RETR ' + fn, open(fn, 'wb').write)
    
ftp.quit()
print "download complete."

    
##########################################################################################
# END PROGRAM
##########################################################################################
