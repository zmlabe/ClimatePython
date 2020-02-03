"""
Figures plots sea ice anomalies for August 2019
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 27 January 2020
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import cmocean
from netCDF4 import Dataset
import calc_SeaIceConc_CDRv3 as ICE

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Figures/'
years = np.arange(1982,2019+1,1)

### Read in sea ice data
sic,yearssic,latsic,lonsic = ICE.readSIC(7,True)

datasst = Dataset(directorydata + 'sst.mnmean.nc')
sstn = datasst.variables['sst'][1:] # degrees C
latsst = datasst.variables['lat'][:]
lonsst = datasst.variables['lon'][:]
datasst.close()

### Reshape for [year,month,lat,lon]
lonsst2,latsst2 = np.meshgrid(lonsst,latsst)
sst = np.reshape(sstn,(sstn.shape[0]//12,12,latsst.shape[0],lonsst.shape[0]))

### Mask data
sic[np.where(sic > 100)] = np.nan
siccopy = sic.copy()
siccopy[np.where(siccopy < 15.)] = 0
sic[np.where(sic<15.)]=np.nan
sic[np.where(sic>=15.)]=100.

### Calculate climatology (1982-2010)
yearq = np.where((years >= 1982) & (years <= 2010))[0]
meansst = np.nanmean(sst[yearq,:,:,:],axis=0)
medianice = np.nanmedian(siccopy[yearq,:,:],axis=0)

### Calculate anomalies
anomsst = sst - meansst
anomsic = sic - medianice

### Select month
monsst = anomsst[:,7,:,:]
monice = sic[:,:,:]
climoice = medianice[:,:]

###############################################################################
###############################################################################
###############################################################################
### Select year
sstyear = monsst[-1,:,:] # 2019
sstyear18 = monsst[-2,:,:] # 2018
sicyear = monice[-1,:,:] # 2019
sicyear18 = monice[-2,:,:] # 2018

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
def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon(list(xy),edgecolor='k',facecolor='none',alpha=1,
                   linewidth=1,linestyle='-',zorder=4,clip_on=True)
    plt.gca().add_patch(poly)

###########################################################################
### Set limits for contours and colorbars
limit = np.arange(-3,3.1,0.25)
barlim = np.arange(-3,4,1)

###########################################################################
fig = plt.figure()
ax1 = plt.subplot(121)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

###########################################################################
varn, lons_cyclic = addcyclic(sstyear, lonsst)
lon2d, lat2d = np.meshgrid(lons_cyclic, latsst)
x, y = m(lon2d, lat2d)

###########################################################################          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

###########################################################################
cs1 = m.contourf(lonsic,latsic,sicyear,alpha=1,zorder=3,colors='w',latlon=True)
cs = m.contourf(x,y,varn,limit,extend='both',zorder=2)
cs3 = m.contour(lonsic,latsic,climoice,
                    np.arange(15,30,15),alpha=1,
                    linewidths=1.2,colors='gold',zorder=3,latlon=True)

###########################################################################
m.drawcoastlines(color='darkgray',linewidth=0.3,zorder=6)
m.fillcontinents(color='dimgrey',lake_color='darkgrey',zorder=5)

latsb = [ 76.4, 79.4, 79.4, 76.5 ]
lonsb = [ 38, 38, 60, 60 ]
draw_screen_poly(latsb,lonsb,m)
latsc = [ 68, 74, 74, 68 ]
lonsc = [ 180, 180, 200, 200 ]
draw_screen_poly(latsc,lonsc,m)
latsc = [ 54, 64, 64, 54 ]
lonsc = [ 180, 180, 200, 200 ]
draw_screen_poly(latsc,lonsc,m)

###########################################################################
cmap = cmocean.cm.balance
cs.set_cmap(cmap)

###########################################################################
###########################################################################
###########################################################################
ax1 = plt.subplot(122)

m = Basemap(projection='npstere',boundinglat=60,lon_0=0,
            resolution='l',round =True)

###########################################################################
varn, lons_cyclic = addcyclic(sstyear18, lonsst)
lon2d, lat2d = np.meshgrid(lons_cyclic, latsst)
x, y = m(lon2d, lat2d)

###########################################################################          
circle = m.drawmapboundary(fill_color='darkgrey',
                           color='dimgrey',linewidth=0.7)
circle.set_clip_on(False)

###########################################################################
cs1 = m.contourf(lonsic,latsic,sicyear18,alpha=1,zorder=3,colors='w',latlon=True)
cs = m.contourf(x,y,varn,limit,extend='both',zorder=2)
cs3 = m.contour(lonsic,latsic,climoice,
                    np.arange(15,30,15),alpha=1,
                    linewidths=1.2,colors='gold',zorder=3,latlon=True)

###########################################################################
m.drawcoastlines(color='darkgray',linewidth=0.3,zorder=6)
m.fillcontinents(color='dimgrey',lake_color='darkgrey',zorder=5)

latsb = [ 76.4, 79.4, 79.4, 76.5 ]
lonsb = [ 38, 38, 60, 60 ]
draw_screen_poly(latsb,lonsb,m)
latsc = [ 68, 74, 74, 68 ]
lonsc = [ 180, 180, 200, 200 ]
draw_screen_poly(latsc,lonsc,m)
latsc = [ 54, 64, 64, 54 ]
lonsc = [ 180, 180, 200, 200 ]
draw_screen_poly(latsc,lonsc,m)

###########################################################################
cmap = cmocean.cm.balance
cs.set_cmap(cmap)

###########################################################################
plt.tight_layout()

###########################################################################
plt.annotate(r'\textbf{[a]}',textcoords='axes fraction',
            xy=(0,0), xytext=(-1.0,0.93),
            fontsize=11,color='dimgrey')
plt.annotate(r'\textbf{[b]}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.,0.93),
            fontsize=11,color='dimgrey')
plt.annotate(r'\textbf{AUGUST 2018}',textcoords='axes fraction',
            xy=(0,0), xytext=(0.5,1.05),
            fontsize=17,color='k',ha='center',va='center')
plt.annotate(r'\textbf{AUGUST 2019}',textcoords='axes fraction',
            xy=(0,0), xytext=(-0.55,1.05),
            fontsize=17,color='k',ha='center',va='center')
plt.annotate(r'\textbf{Relative to 1982-2010}',textcoords='axes fraction',
            xy=(0,0), xytext=(-0.015,0),
            fontsize=6,color='k',ha='center',va='center')
###########################################################################

cbar_ax = fig.add_axes([0.304,0.11,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)

############################################################################
cbar.set_label(r'\textbf{SST Anomaly [$^\circ$C]}',fontsize=11,color='dimgrey')  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelcolor='k')
cbar.outline.set_edgecolor('dimgrey')
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=6)

###########################################################################
plt.savefig(directoryfigure + 'SST-SIC_%s_%s_Anomalies.png' % (8,2019),
            dpi=600)
print('Completed: Script done!')