"""
Plots coupled SST-ICE graphic for introductory talk at 2020's ocean science
meeting for Mary-Louise Timmermans
 
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

### Define constants
directorydata = '/home/zlabe/Documents/Projects/Tests/oisstv2/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/Tests/oisstv2/Figures_AGU20/v3/'
years = np.arange(1982,2019+1,1)

datasic = Dataset(directorydata + 'icec.mnmean.nc')
sicn = datasic.variables['icec'][1:] # 0-100
latsic = datasic.variables['lat'][:]
lonsic = datasic.variables['lon'][:]
datasic.close()

datasst = Dataset(directorydata + 'sst.mnmean.nc')
sstn = datasst.variables['sst'][1:] # degrees C
latsst = datasst.variables['lat'][:]
lonsst = datasst.variables['lon'][:]
datasst.close()

### Reshape for [year,month,lat,lon]
lonsic2,latsic2 = np.meshgrid(lonsic,latsic)
lonsst2,latsst2 = np.meshgrid(lonsst,latsst)
sic = np.reshape(sicn,(sicn.shape[0]//12,12,latsic.shape[0],lonsic.shape[0]))
sst = np.reshape(sstn,(sstn.shape[0]//12,12,latsst.shape[0],lonsst.shape[0]))

siccopy = sic.copy()
siccopy[np.where(siccopy <= 15.)] = 0

### Calculate climatology (1982-2010)
yearq = np.where((years >= 1982) & (years <= 2010))[0]
meansst = np.nanmean(sst[yearq,:,:,:],axis=0)
medianice = np.nanmedian(siccopy[yearq,:,:,:],axis=0)

### Calculate anomalies
anomsst = sst - meansst
anomsic = sic - medianice

### Select month
monq = 8 # index for September
monsst = anomsst[:,monq,:,:]
monice = sic[:,monq,:,:]
climoice = medianice[monq,:,:]

### Mask data
sic[np.where(sic<15.)]=np.nan
sic[np.where(sic>=15.)]=100.

##############################################################################
##############################################################################
##############################################################################
### Define parameters (dark)
def setcolor(x, color):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')

## Plot global temperature anomalies
style = 'polar'

### Define figure
if style == 'ortho':
    m = Basemap(projection='ortho',lon_0=-90,
                lat_0=90,resolution='l',round=True,area_thresh=10000)
elif style == 'polar':
    m = Basemap(projection='npstere',boundinglat=52,lon_0=0,
                resolution='l',round =True)

for i in range(years.shape[0]):
    fig = plt.figure()
    ax = plt.subplot(111)
    for txt in fig.texts:
        txt.set_visible(False)
    
    varsst = monsst[i,:,:]
    varsic = monice[i,:,:]
    
    varsst , lons_cyclic = addcyclic(varsst,lonsst)
    varsic , lons_cyclic = addcyclic(varsic,lonsic)
    varsicmean , lons_cyclic = addcyclic(climoice,lonsic)
    lon2d, lat2d = np.meshgrid(lons_cyclic, latsst)
    x, y = m(lon2d, lat2d)
    
    m.drawmapboundary(fill_color='k')
    
    # Make the plot continuous
    barlim = np.arange(-3,4,1)
    
    cs1 = m.contourf(x,y,varsic,100,extend='both',alpha=1,zorder=2,
                     colors='gainsboro')
    cs2 = m.contourf(x,y,varsst,np.arange(-3,3.01,0.25),extend='both',alpha=1,
                     zorder=1)
    cs3 = m.contour(x,y,varsicmean,
                    np.arange(15,30,15),alpha=1,
                    linewidths=2,colors='gold',zorder=3)
    
    m.drawlsmask(land_color='dimgrey',ocean_color='k')
    m.drawcoastlines(color='k',linewidth=0.1,zorder=6)
    m.fillcontinents(color='dimgrey',lake_color='dimgrey',zorder=5)
    parallels = np.arange(50,86,5)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.0,
                    color='w')
    par=m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0.0,
                        fontsize=4,color='w')
    setcolor(par,'darkgrey')
    
    ### Set colormap
    cmap2 = cmocean.cm.balance                             
    cs2.set_cmap(cmap2)
                
    cbar = m.colorbar(cs2,drawedges=False,location='bottom',pad = 0.14,size=0.1)
    ticks = np.arange(0,8,1)
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.set_label(r'\textbf{SEA SURFACE TEMPERATURE ANOMALIES [\textbf{$\bf{^\circ}$C}]}',
                             fontsize=14,
                             color='w')
    cbar.ax.tick_params(axis='x', size=.0001)
    cbar.outline.set_edgecolor('darkgrey')
#    cbar.dividers.set_color('darkgrey')
#    cbar.dividers.set_linewidth(0.6)
    cbar.outline.set_linewidth(0.6)
    cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    
    if i == 37:
        ccc = 'w'
    else:
        ccc = 'darkgrey'
    t1 = plt.annotate(r'\textbf{%s}' % years[i],textcoords='axes fraction',
                xy=(0,0), xytext=(-0.45,0.83),
            fontsize=50,color=ccc)
    
    plt.annotate(r'\textbf{1982-2010 : Median Sea Ice Extent : Yellow}',textcoords='axes fraction',
                xy=(0,0), xytext=(0.93,0.905),
            fontsize=5,color='gold')
    
###########################################################################
###########################################################################
    if i < 10:        
        plt.savefig(directoryfigure + 'icy_00%s.png' % i,dpi=300)
    elif i >= 10 and i < 100:
        if i < 25: 
            plt.savefig(directoryfigure + 'icy_0%s.png' % i,dpi=300)
        if i == 25: # 2007
            plt.savefig(directoryfigure + 'icy_0%s.png' % i,dpi=300)
            plt.savefig(directoryfigure + 'icy_026.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_027.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_028.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_029.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_030.png',dpi=300)
        elif i >= 26 and i < 30:
            plt.savefig(directoryfigure + 'icy_0%s.png' % (i+5),dpi=300)
#        elif i == 30: # 2012
#            plt.savefig(directoryfigure + 'icy_0%s.png' % (i+5),dpi=300)
#            plt.savefig(directoryfigure + 'icy_036.png',dpi=300)
#            plt.savefig(directoryfigure + 'icy_037.png',dpi=300)
#            plt.savefig(directoryfigure + 'icy_038.png',dpi=300)
#            plt.savefig(directoryfigure + 'icy_039.png',dpi=300)
#            plt.savefig(directoryfigure + 'icy_040.png',dpi=300)
        elif i >= 30:
            plt.savefig(directoryfigure + 'icy_0%s.png' % (i+5),dpi=300)
    else:
        plt.savefig(directoryfigure + 'icy_%s.png' % i,dpi=300)
    if i == 37:
        plt.savefig(directoryfigure + 'icy_990.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_991.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_992.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_993.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_994.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_995.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_996.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_997.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_998.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_999.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9991.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9992.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9993.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9994.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9995.png',dpi=300)
        plt.savefig(directoryfigure + 'icy_9996.png',dpi=300)
