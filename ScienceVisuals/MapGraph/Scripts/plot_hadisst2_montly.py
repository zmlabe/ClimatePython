"""
HadISST2 plot of a particularly month time series

Website   : http://apps.ecmwf.int/data-catalogues/era5/?class=ea
Author    : Zachary M. Labe
Date      : 30 March 2018
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import numpy as np
import datetime
import cmocean

### Define constants
directorydata = '/home/zlabe/dataall/SST/HadISST2/ERA5/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/' \
                'ScienceVisuals/MapGraph/Figures/'
now = datetime.datetime.now()
month = now.month

filename = directorydata + 'SST_ERA5.nc'
data = Dataset(filename)
lats = data.variables['latitude'][:]
lons = data.variables['longitude'][:]
var = data.variables['sst'][5,:,:]
data.close()

years = np.arange(2010,2017+1,1)
yearsq = np.repeat(years,12)

#### Convert
var = var - 273.15

### Define parameters (dark)
def setcolor(x, color):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)

###############################################################################
###############################################################################
###############################################################################
### Define figure 1
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',labelcolor='white')
plt.rc('axes',facecolor='black')             
             
m = Basemap(projection='merc',llcrnrlat=-84,urcrnrlat=84,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='l')

fig = plt.figure()
ax = plt.subplot(111)
for txt in fig.texts:
    txt.set_visible(False)

m.drawmapboundary(fill_color='k')
m.drawlsmask(land_color='w',ocean_color='k')
m.drawcoastlines(color='k',linewidth=0.4)

# Make the plot continuous
barlim = np.arange(0,32,30)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var,np.arange(-1.8,31.5,2),
                extend='both',
                alpha=1)
       
cs.set_cmap('jet')
        
t5 = plt.annotate(r'Sea Surface Temperature',
     textcoords='axes fraction',
     xy=(0,0),xytext=(.5,-0.07),
            fontsize=16,color='w',alpha=1,rotation=0,
            va='center',ha='center')

#    t6 = plt.annotate(r'\textbf{$\bf{^\circ}$\textbf{C}}',
#             textcoords='axes fraction',
#             xy=(0,0), xytext=(.78,-0.164),
#                fontsize=15,color='w',alpha=1,rotation=0,va='center',ha='center')
            
m.fillcontinents(color='w')
            
cbar = plt.colorbar(cs,drawedges=False,orientation='horizontal',pad = 0.1,
                  fraction=0.035)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(axis='x', size=.01)
cbar.ax.tick_params(labelsize=6) 
       
plt.savefig(directoryfigure + 'sst_01.png',dpi=300)
             
###############################################################################
###############################################################################
###############################################################################
### Define figure 2
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',labelcolor='white')
plt.rc('axes',facecolor='black')             
             
m = Basemap(projection='robin',lon_0=0,resolution='l')

fig = plt.figure()
ax = plt.subplot(111)
for txt in fig.texts:
    txt.set_visible(False)

m.drawmapboundary(fill_color='k')
m.drawlsmask(land_color='w',ocean_color='k')
m.drawcoastlines(color='k',linewidth=0.4)

# Make the plot continuous
barlim = np.arange(0,32,30)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var,np.arange(-1.8,31.5,0.2),
                extend='both',
                alpha=1)
       
cs.set_cmap('jet')
        
t5 = plt.annotate(r'Sea Surface Temperature',
     textcoords='axes fraction',
     xy=(0,0),xytext=(.5,-0.07),
            fontsize=16,color='w',alpha=1,rotation=0,
            va='center',ha='center')

#    t6 = plt.annotate(r'\textbf{$\bf{^\circ}$\textbf{C}}',
#             textcoords='axes fraction',
#             xy=(0,0), xytext=(.78,-0.164),
#                fontsize=15,color='w',alpha=1,rotation=0,va='center',ha='center')
            
m.fillcontinents(color='k')
            
cbar = plt.colorbar(cs,drawedges=False,orientation='horizontal',pad = 0.1,
                  fraction=0.035)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(axis='x', size=.01)
cbar.ax.tick_params(labelsize=6) 
       
plt.savefig(directoryfigure + 'sst_02.png',dpi=300)
             
###############################################################################
###############################################################################
###############################################################################
### Define figure 3
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')             
             
m = Basemap(projection='robin',lon_0=0,resolution='l')

fig = plt.figure()
ax = plt.subplot(111)
for txt in fig.texts:
    txt.set_visible(False)

m.drawmapboundary(fill_color='k')
m.drawlsmask(land_color='w',ocean_color='k')
m.drawcoastlines(color='k',linewidth=0.4)

# Make the plot continuous
barlim = np.arange(0,32,30)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var,np.arange(-1.8,31.5,0.2),
                extend='both',
                alpha=1)
       
cs.set_cmap('jet')
        
t5 = plt.annotate(r'\textbf{SEA SURFACE TEMPERATURE}',
     textcoords='axes fraction',
     xy=(0,0),xytext=(.5,-0.07),
            fontsize=16,color='darkgrey',alpha=1,rotation=0,
            va='center',ha='center')
t6 = plt.annotate(r'\textbf{$\bf{^\circ}$\textbf{C}}',
         textcoords='axes fraction',
         xy=(0,0), xytext=(.78,-0.164),
            fontsize=15,color='darkgrey',alpha=1,rotation=0,va='center',ha='center')
            
m.fillcontinents(color='k')
            
cbar = plt.colorbar(cs,drawedges=False,orientation='horizontal',pad = 0.1,
                  fraction=0.035)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(axis='x', size=.01)
cbar.ax.tick_params(labelsize=11) 
       
plt.savefig(directoryfigure + 'sst_03.png',dpi=300)
             
###############################################################################
###############################################################################
###############################################################################
### Define figure 4
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')             
             
m = Basemap(projection='robin',lon_0=0,resolution='l')

fig = plt.figure()
ax = plt.subplot(111)
for txt in fig.texts:
    txt.set_visible(False)

m.drawmapboundary(fill_color='k')
m.drawlsmask(land_color='w',ocean_color='k')
m.drawcoastlines(color='k',linewidth=0.4)

# Make the plot continuous
barlim = np.arange(0,32,30)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var,np.arange(-1.8,31.5,0.2),
                extend='both',
                alpha=1)
       
cs.set_cmap(cmocean.cm.thermal)
        
t5 = plt.annotate(r'\textbf{SEA SURFACE TEMPERATURE}',
     textcoords='axes fraction',
     xy=(0,0),xytext=(.5,-0.07),
            fontsize=16,color='darkgrey',alpha=1,rotation=0,
            va='center',ha='center')
t6 = plt.annotate(r'\textbf{$\bf{^\circ}$\textbf{C}}',
         textcoords='axes fraction',
         xy=(0,0), xytext=(.78,-0.164),
            fontsize=15,color='darkgrey',alpha=1,rotation=0,va='center',ha='center')
            
m.fillcontinents(color='k')
            
cbar = plt.colorbar(cs,drawedges=False,orientation='horizontal',pad = 0.1,
                  fraction=0.035)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(axis='x', size=.01)
cbar.ax.tick_params(labelsize=11) 
       
plt.savefig(directoryfigure + 'sst_04.png',dpi=300)
             
###############################################################################
###############################################################################
###############################################################################
### Define figure 5
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')             
             
m = Basemap(projection='moll',lon_0=0,resolution='l',area_thresh=10000)

fig = plt.figure()
ax = plt.subplot(111)
for txt in fig.texts:
    txt.set_visible(False)

m.drawmapboundary(fill_color='k')
m.drawlsmask(land_color='k',ocean_color='k')
m.drawcoastlines(color='k',linewidth=0.4)

# Make the plot continuous
barlim = np.arange(0,32,30)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var,np.arange(-1.8,31.5,0.2),
                extend='both',
                alpha=1)
       
cs.set_cmap(cmocean.cm.thermal)
        
t5 = plt.annotate(r'\textbf{SEA SURFACE TEMPERATURE}',
     textcoords='axes fraction',
     xy=(0,0),xytext=(.5,-0.07),
            fontsize=16,color='darkgrey',alpha=1,rotation=0,
            va='center',ha='center')
t6 = plt.annotate(r'\textbf{$\bf{^\circ}$\textbf{C}}',
         textcoords='axes fraction',
         xy=(0,0), xytext=(.78,-0.164),
            fontsize=15,color='darkgrey',alpha=1,rotation=0,va='center',ha='center')
            
m.fillcontinents(color='k')
            
cbar = plt.colorbar(cs,drawedges=False,orientation='horizontal',pad = 0.1,
                  fraction=0.035)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))  
cbar.ax.tick_params(axis='x', size=.01)
cbar.ax.tick_params(labelsize=11) 
       
plt.savefig(directoryfigure + 'sst_05.png',dpi=300)