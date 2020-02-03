"""
Script calculates time series mean for oisstv2 data
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 20 January 2020
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import calc_Utilities as UT
import scipy.stats as sts

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Figures/'
years = np.arange(1982,2019+1,1)

datasst = Dataset(directorydata + 'sst.mnmean.nc')
sstNO = datasst.variables['sst'][1:] # degrees C
lat1 = datasst.variables['lat'][:]
lon1= datasst.variables['lon'][:]
datasst.close()

dataland = Dataset(directorydata + 'lsmask.nc')
mask = dataland.variables['mask'][:]
dataland.close()
sstn = sstNO * mask
sstn[np.where(sstn == 0.)] = np.nan

### Reshape for [year,month,lat,lon]
sst = np.reshape(sstn,(sstn.shape[0]//12,12,lat1.shape[0],lon1.shape[0]))

### Calculate average for polar cap and mask data
def aveSST(data,lat,lon,mask,period):
    if mask == 'polar_cap':
        latpolar = 67.
        latq = np.where((lat >= latpolar))[0]
        sstnew = sst[:,:,latq,:]
        latt = lat[latq]
        lonn = lon
    elif mask == 'north_barents':
        latq = np.where((lat >= 76.4) & (lat <= 79.4))[0]
        lonq = np.where((lon >= 38) & (lon <= 60))[0]
        latt = lat[latq]
        lonn = lon[lonq]
        sstnew1 = sst[:,:,latq,:]
        sstnew = sstnew1[:,:,:,lonq]
    elif mask == 'chukchi':
        latq = np.where((lat >= 68) & (lat <= 74))[0]
        lonq = np.where((lon >= 180) & (lon <= 200))[0]
        latt = lat[latq]
        lonn = lon[lonq]
        sstnew1 = sst[:,:,latq,:]
        sstnew = sstnew1[:,:,:,lonq]
    elif mask == 'bering':
        latq = np.where((lat >= 54) & (lat <= 64))[0]
        lonq = np.where((lon >= 180) & (lon <= 200))[0]
        latt = lat[latq]
        lonn = lon[lonq]
        sstnew1 = sst[:,:,latq,:]
        sstnew = sstnew1[:,:,:,lonq]
    elif mask == 'global':
        sstnew = sst
        latt = lat
        lonn = lon
        
    ### Mask lat/lon for 2d
    lon2,lat2 = np.meshgrid(lonn,latt)
    
    ### Calculated weighted average
    sstave = UT.calc_weightedAve(sstnew,lat2)
    plt.contourf(sstnew[0,0,:,:])
    plt.show()
    
    ### Calculate temporal period
    if period == 'annual':
        sstperiod = np.nanmean(sstave,axis=1)
    elif period == 'aug':
        sstperiod = sstave[:,7]
    elif period == 'sep':
        sstperiod = sstave[:,8]
        
    ### Calculate anomalies using 1982-2010 baseline
    years = np.arange(1982,2019+1,1)
    yearq = np.where((years >= 1982) & (years <= 2010))[0]
    
    climo = np.nanmean(sstperiod[yearq])
    anomperiod = sstperiod - climo
    
    return sstave,anomperiod,latt,lonn

### Use functions
periodmonth = 'aug'
meansst_global,period_global,latt,lonn = aveSST(sst,lat1,lon1,'global',periodmonth)
meansst_polar,period_polar,latt,lonn = aveSST(sst,lat1,lon1,'polar_cap',periodmonth)
meansst_barents,period_barents,latt,lonn = aveSST(sst,lat1,lon1,'north_barents',periodmonth)
meansst_chukchi,period_chukchi,latt,lonn = aveSST(sst,lat1,lon1,'chukchi',periodmonth)
meansst_bering,period_bering,latt,lonn = aveSST(sst,lat1,lon1,'bering',periodmonth)

### Calculate trends
regions = ['Global','Arctic','Barents','Chukchi','Bering']
dataall = [period_global,period_polar,period_barents,period_chukchi,period_bering]
print('\n----------------------------')

slope = np.empty((len(dataall)))
intercept = np.empty((len(dataall)))
rvalue = np.empty((len(dataall)))
pvalue = np.empty((len(dataall)))
std_err = np.empty((len(dataall)))
for i in range(len(dataall)):
    time = np.arange(dataall[i].shape[0])
    slope[i],intercept[i],rvalue[i],pvalue[i],std_err[i] = sts.linregress(time,dataall[i])
    print('Trend is *%s* per year for *%s*' % (np.round(slope[i],3),regions[i]))
print('----------------------------\n')  

### Confidence interval for 95% level
std_err = std_err*1.96

#### Saving files for time series data
#for i in range(len(dataall)):
#    np.savetxt(directorydata + 'SST_%s_%s_1982-2019.txt' % (regions[i],periodmonth),
#               dataall[i])
#    print('Completed: Saved data for %s!' % regions[i])
