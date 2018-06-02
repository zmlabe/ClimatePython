"""
Script plots temperature anomaly over the Arctic for a set time period

Notes
-----
    Author : Zachary Labe
    Date   : 2 June 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import cmocean

### Define directories
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/' \
                'ScienceVisuals/MapGraph/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/' \
                'ScienceVisuals/MapGraph/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting polar temperature anomalies - %s----' % titletime)

### Alott time series
year1 = 1948
year2 = 2018
years = np.arange(year1,year2+1,1)

### Retrieve data
data = Dataset(directorydata + 'May_2018_R1_T925.nc')
time = data.variables['time'][:]
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
temp = data.variables['air'][:]
data.close()

### Define parameters (dark)
def setcolor(x, color):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)

###############################################################################
###############################################################################
###############################################################################
### Plot map of T anomalies 1
plt.rc('savefig',facecolor='k')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='w')
plt.rc('ytick',color='w')
plt.rc('axes',labelcolor='w')
plt.rc('axes',facecolor='k')        
             
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
            round =True)

ax = plt.subplot(111)

var = temp[0,:,:]

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=1.1)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-5,6,5)

var, lons_cyclic = addcyclic(var, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var[:,:],
                np.arange(-5,5.1,0.1),extend='both')

cs.set_cmap('jet')

plt.title(r'925 mb Temperature Anomaly - MAY 2018',
         color='w',size=17)

cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.05,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001)

fig.subplots_adjust(wspace=-0.65)
fig.subplots_adjust(hspace=-0.01)
    
#plt.savefig(directoryfigure + 'May_2018_R1_T925_01.png',dpi=300)
#
################################################################################
################################################################################
################################################################################
#### Plot map of T anomalies 2
#plt.rc('savefig',facecolor='k')
#plt.rc('axes',edgecolor='k')
#plt.rc('xtick',color='w')
#plt.rc('ytick',color='w')
#plt.rc('axes',labelcolor='w')
#plt.rc('axes',facecolor='k')        
#             
#fig = plt.figure()
#ax = plt.subplot(111)
#
#m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
#            round =True)
#
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='k',linewidth=1.1)
#m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#
## Make the plot continuous
#barlim = np.arange(-5,6,5)
#
#cs = m.contourf(x,y,var[:,:],
#                np.arange(-5,5.1,0.1),extend='both')
#
#cs.set_cmap('seismic')
#
#plt.title(r'925 mb Temperature Anomaly - MAY 2018',
#         color='w',size=17)
#
#cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='both',extendfrac=0.05,drawedges=False)
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(list(map(str,barlim)))
#cbar.ax.tick_params(axis='x', size=.001)
#
#fig.subplots_adjust(wspace=-0.65)
#fig.subplots_adjust(hspace=-0.01)
#    
#plt.savefig(directoryfigure + 'May_2018_R1_T925_02.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot map of T anomalies 3
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='k')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='k')        
             
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
            round=True)

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=1.1)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-5,6,5)

cs = m.contourf(x,y,var[:,:],
                np.arange(-5,5.1,0.1),extend='both')

cs.set_cmap('seismic')

plt.title(r'\textbf{925 mb TEMPERATURE ANOMALY - MAY 2018}',
         color='darkgrey',size=17)

cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.05,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001)

plt.text(-0.18,-0.093,r'\textbf{$^\circ$C}',
         color='darkgrey',size=15)

fig.subplots_adjust(wspace=-0.65)
fig.subplots_adjust(hspace=-0.01)
    
plt.savefig(directoryfigure + 'May_2018_R1_T925_03.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot map of T anomalies 4
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='k')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='k')        
             
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
            round=True,area_thresh=10000.)

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=1.1)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-5,6,5)

cs = m.contourf(x,y,var[:,:],
                np.arange(-5,5.1,0.1),extend='both')

cs.set_cmap('seismic')

plt.title(r'\textbf{925 mb TEMPERATURE ANOMALY - MAY 2018}',
         color='darkgrey',size=17)

cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.05,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001)

plt.text(-0.18,-0.093,r'\textbf{$^\circ$C}',
         color='darkgrey',size=15)

fig.subplots_adjust(wspace=-0.65)
fig.subplots_adjust(hspace=-0.01)
    
plt.savefig(directoryfigure + 'May_2018_R1_T925_04.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot map of T anomalies 5
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='k')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='k')        
             
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
            round=True,area_thresh=10000.)

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=1.1)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-5,6,5)

cs = m.contourf(x,y,var[:,:],
                np.arange(-5,5.1,0.1),extend='both')

cs.set_cmap(cmocean.cm.balance)

plt.title(r'\textbf{925 mb TEMPERATURE ANOMALY - MAY 2018}',
         color='darkgrey',size=17)

cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.05,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001)

plt.text(-0.18,-0.093,r'\textbf{$^\circ$C}',
         color='darkgrey',size=15)

fig.subplots_adjust(wspace=-0.65)
fig.subplots_adjust(hspace=-0.01)
    
plt.savefig(directoryfigure + 'May_2018_R1_T925_05.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot map of T anomalies 6
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='k')
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='k')        
             
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=54,lon_0=270,resolution='l',
            round=True,area_thresh=10000.)

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='dimgrey',linewidth=1.1)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
parallels = np.arange(50,86,5)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.001,color='darkgrey')
par=m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0.001,fontsize=4.5,color='darkgrey')
setcolor(par,'darkgrey')

# Make the plot continuous
barlim = np.arange(-5,6,5)

cs = m.contourf(x,y,var[:,:],
                np.arange(-5,5.1,0.1),extend='both')

cs.set_cmap(cmocean.cm.balance)

cbar_ax = fig.add_axes([0.315,0.05,0.4,0.035])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.05,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001)

plt.text(-.70,24.75,r'\textbf{925 mb TEMPERATURE ANOMALY : \underline{MAY 2018}}',
         color='darkgrey',size=17)
plt.text(-.18,-0.093,r'\textbf{$^\circ$C}',
         color='darkgrey',size=15)

fig.subplots_adjust(wspace=-0.65)
fig.subplots_adjust(hspace=-0.01)
    
plt.savefig(directoryfigure + 'May_2018_R1_T925_06.png',dpi=300)