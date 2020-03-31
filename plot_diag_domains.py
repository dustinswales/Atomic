import os,netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime

# Read in data
fileIN = 'clavrx_01min_10min_domain_status.nc'
dataIN = netCDF4.Dataset(fileIN,'r')
year   = dataIN.variables['year'][:]
month  = dataIN.variables['month'][:]
day    = dataIN.variables['day'][:]
hour   = dataIN.variables['hour'][:]
minute = dataIN.variables['minute'][:]
flag01 = dataIN.variables['flag01'][:]
flag02 = dataIN.variables['flag02'][:]
flag03 = dataIN.variables['flag03'][:]
ntime  = len(year)

# Create arrays for plotting
time = np.empty(ntime,dtype=float)
for iTime in range(0,ntime):
    # Time
    deltaDays   = datetime.date(year[iTime],month[iTime],day[iTime]) - datetime.date(year[iTime],1,1)
    doy         = deltaDays.days + 1
    time[iTime] = year[iTime] + float(doy)/365. + hour[iTime]/(24.*365.) + minute[iTime]/(24.*60.*365.)

#########################################################################
# 01-minute data
#########################################################################
# Where is the first time w/ 01-minute data?
a = np.min(np.where(flag01 == 1))

# Is the data continuous? If not, how many discontinuities?
count_01discon = 0
init = 1
status = 0 # Initially assume always in a continuous data block...
for iTime in range(a+1,ntime):
    # Found a data gap...
    if (flag01[iTime] != flag01[iTime-1] and not flag01[iTime]):
        #print("Found a data gap, starting at ",str(iTime))
        count_01discon = count_01discon + 1
        if (init):
            xrange01 = [a,iTime]
            xrange01 = np.stack((xrange01,[iTime+1,ntime]))
            init = 0
        else:
            xrange01[count_01discon-1,1] = iTime - 1
            xrange01 = np.vstack((xrange01,[iTime+1,ntime]))
        status = 1
        
    # Exiting data gap...
    if (flag01[iTime] != flag01[iTime-1] and flag01[iTime] and status == 1):
        #print("                  ending at   ",str(iTime))
        xrange01[count_01discon,0] = iTime
        t0 = iTime
        status = 0
# The plotting routine expects the start:stride, vs start:stop
for ij in range(0,count_01discon):
    xrange01[ij,1] = xrange01[ij,1]-xrange01[ij,0] + 1
    #print(flag01[xrange01[ij,0]-1:xrange01[ij,1]+1])

#########################################################################
# 10-minute data
#########################################################################
# Where is the first time w/ 10-minute data?
a = np.min(np.where(flag02 == 1))

# Is the data continuous? If not, how many discontinuities?
count_10discon = 0
init = 1
status = 0 # Initially assume always in a continuous data block...
for iTime in range(a+1,ntime):
    # Found a data gap...
    if (flag02[iTime] != flag02[iTime-1] and not flag02[iTime]):
        #print("Found a data gap, starting at ",str(iTime))
        count_10discon = count_10discon + 1
        if (init):
            xrange10 = [a,iTime]
            xrange10 = np.stack((xrange10,[iTime+1,ntime]))
            init = 0
        else:
            xrange10[count_10discon-1,1] = iTime
            xrange10 = np.vstack((xrange10,[iTime+1,ntime]))
        status = 1
        
    # Exiting data gap...
    if (flag02[iTime] != flag02[iTime-1] and flag02[iTime] and status == 1):
        #print("                  ending at   ",str(iTime))
        xrange10[count_10discon,0] = iTime
        t0 = iTime
        status = 0

for ij in range(0,count_10discon):
    xrange10[ij,1] = xrange10[ij,1]-xrange10[ij,0] + 1
    #print(flag02[xrange10[ij,0]-1:xrange10[ij,1]+1])

#########################################################################
# 01-minute data subsets 10-minute?
#########################################################################
# Where is the first time w/ 10-minute data?
a = np.min(np.where(flag03 == 1))

# Is the data continuous? If not, how many discontinuities?
count_10discon = 0
init = 1
status = 0 # Initially assume always in a continuous data block...
for iTime in range(a+1,ntime):
    # Found a data gap...
    if (flag03[iTime] != flag03[iTime-1] and not flag03[iTime]):
        #print("Found a data gap, starting at ",str(iTime))
        count_10discon = count_10discon + 1
        if (init):
            xrangeSS = [a,iTime]
            xrangeSS = np.stack((xrangeSS,[iTime+1,ntime]))
            init = 0
        else:
            xrangeSS[count_10discon-1,1] = iTime
            xrangeSS = np.vstack((xrangeSS,[iTime+1,ntime]))
        status = 1
        
    # Exiting data gap...
    if (flag03[iTime] != flag03[iTime-1] and flag03[iTime] and status == 1):
        #print("                  ending at   ",str(iTime))
        xrangeSS[count_10discon,0] = iTime
        t0 = iTime
        status = 0

for ij in range(0,count_10discon):
    xrangeSS[ij,1] = xrangeSS[ij,1]-xrangeSS[ij,0] + 1


    
fig, ax = plt.subplots()
ax.broken_barh(xrange01, (10,9), facecolors='tab:blue')
ax.broken_barh(xrange10, (20,9), facecolors='tab:green')
ax.broken_barh(xrangeSS, (30,9), facecolors='tab:red')
ax.set_ylim(5, 45)
ax.set_xlim(0,ntime)
ax.set_yticks([15,25,35])
ax.set_yticklabels(['01-minute', '10-minute','subsets'])

# Create x tick labels
stride = 5000
xti    = np.linspace(0, ntime, np.ceil(ntime/stride))
nxt    = len(xti)
xtl    = ["" for x in range(nxt)]
count  = 0
for ij in range(0,ntime,stride):
    xtl[count] = str(month[ij])+"/"+str(day[ij])
    count = count + 1
#
ax.set_xticks(xti)
ax.set_xticklabels(xtl)
plt.show()
