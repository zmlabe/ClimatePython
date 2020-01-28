"""
Script calculates composite month for oisstv2 data
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 27 January 2020
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
import palettable.cubehelix as cm
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

### Read in sea ice data
datasic = Dataset(directorydata + 'icec.mnmean.nc')
sicn = datasic.variables['icec'][1:] # 0-100
latsic = datasic.variables['lat'][:]
lonsic = datasic.variables['lon'][:]
datasic.close()

sic = np.reshape(sicn,(sicn.shape[0]//12,12,latsic.shape[0],lonsic.shape[0]))
sicmonth = sic[:,monthq,:,:]

### Mask data
sic[np.where(sic<15.)]=np.nan
sic[np.where(sic>=15.)]=100.

###############################################################################
###############################################################################
###############################################################################
### Select year
sstyear = sstmonth[-1,:,:] # 2019
sicyear = sicmonth[-1,:,:] # 2019

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

###########################################################################
### Set limits for contours and colorbars
limit = np.arange(-2,10.1,0.5)
barlim = np.arange(-2,11,2)
lonq,latq = np.meshgrid(lon1,lat1)

###########################################################################
fig = plt.figure()
ax1 = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

###########################################################################
varn, lons_cyclic = addcyclic(sstyear, lon1)
varsic , lons_cyclic = addcyclic(sicyear,lonsic)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
x, y = m(lon2d, lat2d)

###########################################################################          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

###########################################################################
cs1 = m.contourf(x,y,varsic,alpha=1,zorder=3,colors='w')
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

cbar = m.colorbar(cs,drawedges=False,location='right',pad = 0.5,size=0.1)

###########################################################################
parallels = np.arange(50,86,10)
meridians = np.arange(-180,180,90)
m.drawparallels(parallels,labels=[True,True,False,False],linewidth=0.5,
                color='k',zorder=10)
par=m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0.5,
                    fontsize=4,color='k',zorder=11)
setcolor(par,'k')

###########################################################################
cbar.set_label(r'\textbf{SST [$^\circ$C]}',fontsize=11,color='dimgrey')  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='y', size=.001,labelcolor='k')
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=6)

###########################################################################
plt.tight_layout()

###########################################################################
plt.annotate(r'\textbf{AUG. 2019}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.73,0.99),
            fontsize=20,color='k')

plt.annotate(r'\textbf{Barents}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.61,0.32),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{GREENLAND}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.275,0.32),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Baffin}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.22,0.39),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Bay}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.23,0.37),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Beaufort}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.28,0.73),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Chukchi}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.4,0.78),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{East}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.53,0.76),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Siberian}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.50,0.74),
            fontsize=6.5,color='k',zorder=11)
plt.annotate(r'\textbf{Laptev}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.62,0.63),
            fontsize=6.5,color='k',zorder=11)

###########################################################################
plt.savefig(directoryfigure + 'SST-SIC_%s_%s.png' % (monthq,2019),dpi=600)
print('Completed: Script done!')