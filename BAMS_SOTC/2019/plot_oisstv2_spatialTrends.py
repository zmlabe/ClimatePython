"""
Script calculates time series mean for oisstv2 data
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 20 January 2020
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
from netCDF4 import Dataset
import calc_Utilities as UT
import scipy.stats as sts

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOC/2019/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOC/2019/Figures/'
years = np.arange(1982,2019+1,1)
time = np.arange(years.shape[0])

datasst = Dataset(directorydata + 'sst.mnmean.nc')
sstn = datasst.variables['sst'][1:] # degrees C
lat1 = datasst.variables['lat'][:]
lon1= datasst.variables['lon'][:]
datasst.close()

### Reshape for [year,month,lat,lon]
sst = np.reshape(sstn,(sstn.shape[0]//12,12,lat1.shape[0],lon1.shape[0]))

### Select month 
monthq = 8 - 1 
sstmonth = sst[:,monthq,:,:]

### Calculate trends over the 1982-2019 period
slope = np.empty((lat1.shape[0],lon1.shape[0]))
intercept = np.empty((lat1.shape[0],lon1.shape[0]))
rvalue = np.empty((lat1.shape[0],lon1.shape[0]))
pvalue = np.empty((lat1.shape[0],lon1.shape[0]))
std_err = np.empty((lat1.shape[0],lon1.shape[0]))
for i in range(lat1.shape[0]):
    for j in range(lon1.shape[0]):
        slope[i,j],intercept[i,j],rvalue[i,j],pvalue[i,j],std_err[i,j] = sts.linregress(time,sstmonth[:,i,j])
        
### Calculate significance
mask = pvalue.copy()
mask[np.where(mask <= 0.05)] = 1.
mask[np.where(mask != 1)] = 0.

### Trends 
trends = slope * mask
trends[np.where(trends == 0.0)] = np.nan

### Read in sea ice data
datasic = Dataset(directorydata + 'icec.mnmean.nc')
sicn = datasic.variables['icec'][1:] # 0-100
latsic = datasic.variables['lat'][:]
lonsic = datasic.variables['lon'][:]
datasic.close()

sic = np.reshape(sicn,(sicn.shape[0]//12,12,latsic.shape[0],lonsic.shape[0]))
siccopy = sic.copy()
siccopy[np.where(siccopy <= 15.)] = 0

yearq = np.where((years >= 1982) & (years <= 2010))[0]
medianice = np.nanmedian(siccopy[yearq,:,:,:],axis=0)

### Select month
varsic = sic[:,monthq,:,:]
climoice = medianice[monthq,:,:]

### Mask data
sic[np.where(sic<15.)]=np.nan
sic[np.where(sic>=15.)]=100.

###########################################################################
###########################################################################
###########################################################################
### Plot Variables
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
def setcolor(x, color):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)

### Set limits for contours and colorbars
limit = np.arange(-0.1,0.101,0.005)
barlim = np.round(np.arange(-0.1,0.11,0.05),2)
lonq,latq = np.meshgrid(lon1,lat1)

fig = plt.figure()
ax1 = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

varn, lons_cyclic = addcyclic(trends, lon1)
varsic , lons_cyclic = addcyclic(varsic,lonsic)
varsicmean , lons_cyclic = addcyclic(climoice,lonsic)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)
          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(x,y,varsic[-1],100,extend='both',alpha=1,zorder=2,
                 colors='w')
cs = m.contourf(x,y,varn,limit,extend='both')
cs2 = m.contour(x,y,varsicmean,
                np.arange(15,30,15),alpha=1,
                linewidths=2,colors='gold',zorder=3)

m.drawcoastlines(color='darkgray',linewidth=0.3,zorder=4)
m.fillcontinents(color='dimgrey',lake_color='darkgrey',zorder=3)

cmap = cmocean.cm.balance
cs.set_cmap(cmap)
cbar = m.colorbar(cs,drawedges=False,location='bottom',pad = 0.14,size=0.1)

parallels = np.arange(50,86,10)
meridians = np.arange(-180,180,90)
m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                color='k',zorder=10)
par=m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0.5,
                    fontsize=4,color='k',zorder=11)
setcolor(par,'k')

cbar.set_label(r'\textbf{SST Linear Trend [$^\circ$C/year]}',fontsize=11,color='dimgrey')  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelcolor='k')
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()

plt.savefig(directoryfigure + 'SST_trends_1982-2019.png',dpi=600)
print('Completed: Script done!')


        