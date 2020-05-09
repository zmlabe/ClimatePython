"""
Script calculates time series mean for oisstv2 data
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 20 January 2020
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
from netCDF4 import Dataset
import calc_Utilities as UT
import scipy.stats as sts
import calc_SeaIceConc_CDRv3 as ICE

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Figures/'
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
monthq = 7 
sstmonth = sst[:,monthq,:,:]
sstyear = sstmonth[-1,:,:]

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
sic,yearssic,latsic,lonsic = ICE.readSIC(7,True)
sicyear = sic[-1,:,:]

### Mask data
sic[np.where(sic > 100)] = np.nan
siccopy = sic.copy()
siccopy[np.where(siccopy < 15.)] = 0
sic[np.where(sic<15.)]=np.nan
sic[np.where(sic>=15.)]=100.

### Calculate climatology (1982-2010)
yearq = np.where((years >= 1982) & (years <= 2010))[0]
medianice = np.nanmedian(siccopy[yearq,:,:],axis=0)

climoice = medianice[:,:]

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
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap      
             
fig = plt.figure()
ax1 = plt.subplot(121)

###########################################################################
### Set limits for contours and colorbars
limit = np.arange(-2,10.1,0.5)
barlim = np.arange(-2,11,2)
lonq,latq = np.meshgrid(lon1,lat1)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

###########################################################################
varn, lons_cyclic = addcyclic(sstyear, lon1)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

###########################################################################          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

###########################################################################
cs1 = m.contourf(lonsic,latsic,sicyear,alpha=1,zorder=3,colors='w',
                 latlon=True)
cs = m.contourf(x,y,varn,limit,extend='both',zorder=2)
csl = m.contour(x,y,varn,np.arange(10,15,5),linestyles='-',
               linewidths=2,colors='k',zorder=4)

###########################################################################
m.drawcoastlines(color='darkgray',linewidth=0.3,zorder=6)
m.fillcontinents(color='dimgrey',lake_color='darkgrey',zorder=5)

###########################################################################
cmap = plt.get_cmap('CMRmap')
new_cmap = truncate_colormap(cmap,0,0.85)
cs.set_cmap(new_cmap)
cs.cmap.set_under('w')

cbar = m.colorbar(cs,drawedges=False,location='bottom',pad = 0.14,size=0.1)

###########################################################################
parallels = np.arange(50,86,10)
meridians = np.arange(-180,180,90)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.5,
                color='k',zorder=10)
par=m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.5,
                    fontsize=4,color='k',zorder=11)
setcolor(par,'k')

###########################################################################
cbar.set_label(r'\textbf{2019 SST [$^\circ$C]}',fontsize=8,color='dimgrey')  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelcolor='k')
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=6)

###########################################################################
plt.annotate(r'\textbf{Barents}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.60,0.325),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{GREENLAND}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.265,0.32),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Baffin}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.215,0.39),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Bay}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.225,0.37),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Beaufort}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.28,0.73),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Chukchi}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.4,0.78),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{East}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.535,0.76),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Siberian}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.50,0.74),
            fontsize=6,color='w',zorder=11)
plt.annotate(r'\textbf{Laptev}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.62,0.63),
            fontsize=6,color='w',zorder=11)

plt.tight_layout()

### Set limits for contours and colorbars
limit = np.arange(-0.1,0.101,0.005)
barlim = np.round(np.arange(-0.1,0.11,0.05),2)
lonq,latq = np.meshgrid(lon1,lat1)

ax1 = plt.subplot(122)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

varn, lons_cyclic = addcyclic(trends, lon1)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)
          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

cs1 = m.contourf(lonsic,latsic,sic[-1],100,extend='both',alpha=1,zorder=2,
                 colors='w',latlon=True)
cs = m.contourf(x,y,varn,limit,extend='both')
cs2 = m.contour(lonsic,latsic,climoice,
                np.arange(15,30,15),alpha=1,
                linewidths=1.2,colors='gold',zorder=3,latlon=True)

m.drawcoastlines(color='darkgray',linewidth=0.3,zorder=4)
m.fillcontinents(color='dimgrey',lake_color='darkgrey',zorder=3)

cmap = cmocean.cm.balance
cs.set_cmap(cmap)
cbar = m.colorbar(cs,drawedges=False,location='bottom',pad = 0.14,size=0.1)

parallels = np.arange(50,86,10)
meridians = np.arange(-180,180,90)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.5,
                color='k',zorder=10)
par=m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.5,
                    fontsize=4,color='k',zorder=11)
setcolor(par,'k')

cbar.set_label(r'\textbf{1982-2019 SST Linear Trend [$^\circ$C year$\bf{^{-1}}$]}',fontsize=8,color='dimgrey')  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelcolor='k')
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=6)

plt.annotate(r'\textbf{[a]}',textcoords='axes fraction',
            xy=(0,0), xytext=(-1.00,0.93),
            fontsize=11,color='dimgrey')
plt.annotate(r'\textbf{[b]}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.,0.93),
            fontsize=11,color='dimgrey')
#plt.annotate(r'\textbf{AUG. 2019}',textcoords='axes fraction',
#            xy=(0,0), xytext=(-0.57,1.04),
#            fontsize=13,color='k',ha='center',va='center')
#plt.annotate(r'\textbf{AUG. 1982-2019}',textcoords='axes fraction',
#            xy=(0,0), xytext=(0.49,1.04),
#            fontsize=13,color='k',ha='center',va='center')

plt.savefig(directoryfigure + 'Fig5.4_R2.ps',dpi=600,format='ps')
plt.savefig(directoryfigure + 'Composite-Trends_1982-2019_R2.png',dpi=600)
print('Completed: Script done!')


        