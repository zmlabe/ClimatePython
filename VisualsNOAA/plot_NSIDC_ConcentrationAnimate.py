"""
Script reads Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS
Passive Microwave Data, Version 1 binary files for select variables and 
regrids according to selected grid style (e.g., NSIDC EASE grid data). 
 
Source : https://nsidc.org/data/nsidc-0051#
Author : Zachary Labe
Date : 6 November 2018
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime
import cmocean

### Define constants
directorydata = '/seley/zlabe/seaice/nsidc/SepMinSIC/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/VisualsNOAA/Concentration/'
now = datetime.datetime.now()
month = now.month

years = np.arange(1979,2018+1,1)

data = Dataset(directorydata + 'NSIDC_09_1979-2018.nc')
lat2 = data.variables['lat'][:]
lon2 = data.variables['lon'][:]
sic = data.variables['sic'][:,:,:]
data.close()

data = Dataset(directorydata + 'NSIDC_09_Median.nc')
mean = data.variables['sic'][:,:]
data.close()

### Define parameters (dark)
def setcolor(x, color):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')

## Plot global temperature anomalies
style = 'polar'

### Define figure
if style == 'ortho':
    m = Basemap(projection='ortho',lon_0=-90,
                lat_0=90,resolution='l',round=True)
elif style == 'polar':
    m = Basemap(projection='npstere',boundinglat=55,lon_0=270,resolution='l',round =True)

for i in range(sic.shape[0]):
    fig = plt.figure()
    ax = plt.subplot(111)
    for txt in fig.texts:
        txt.set_visible(False)
    
    var = sic[i,:,:]

    m.drawmapboundary(fill_color='k')
    m.drawlsmask(land_color='k',ocean_color='k')
    m.drawcoastlines(color='k',linewidth=0.1)
    m.fillcontinents(color='dimgrey')
    
    # Make the plot continuous
    barlim = np.arange(0.1,1.1,1)
    
    cs = m.contourf(lon2,lat2,var,np.arange(0,101,5),
                    alpha=1,latlon=True)
    cs1 = m.contour(lon2,lat2,mean,
                    np.arange(15,30,15),alpha=1,latlon=True,
                    linewidths=2,colors='darkorange')
    
    parallels = np.arange(50,86,5)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.0,
                    color='w')
    par=m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0.0,
                        fontsize=4,color='w')
    setcolor(par,'darkgrey')
    
    cmap = cmocean.cm.ice  
    cs.set_cmap(cmap)
    if i >= 39:
        plt.annotate(r'\textbf{%s}' % years[i],textcoords='axes fraction',
                    xy=(0,0), xytext=(-0.51,0.92),
                fontsize=50,color='w')
    else:
        plt.annotate(r'\textbf{%s}' % years[i],textcoords='axes fraction',
                xy=(0,0), xytext=(-0.51,0.92),
            fontsize=50,color='darkgrey')
    plt.annotate(r'\textbf{1981-2010 Median : Orange}',textcoords='axes fraction',
                xy=(0,0), xytext=(-0.57,-0.22),
            fontsize=8,color='darkorange')
    
    cbar_ax = fig.add_axes([0.315,0.05,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        drawedges=False)
    barlim = np.arange(0,101,50)
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.set_label(r'\textbf{SEA ICE CONCENTRATION [\%]}',fontsize=13,
                         color='darkgrey',labelpad=-36)
    cbar.ax.tick_params(axis='x', size=.001)
    plt.subplots_adjust(bottom=0.18)

    if i < 10:        
        plt.savefig(directoryfigure + 'icy_0%s.png' % i,dpi=300)
    else:
        plt.savefig(directoryfigure + 'icy_%s.png' % i,dpi=300)
        if i == 39:
            plt.savefig(directoryfigure + 'icy_991.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_992.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_993.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_994.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_995.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_996.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_997.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_998.png',dpi=300)
            plt.savefig(directoryfigure + 'icy_999.png',dpi=300)