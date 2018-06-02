"""
Script plots the global temperatures from various data sets

Notes
-----
    Author : Zachary Labe
    Date   : 2 June 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import datetime
import cmocean

### Define directories
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/' \
                'ScienceVisuals/LineGraph/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/' \
                'ScienceVisuals/LineGraph/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting global temperatures - %s----' % titletime)

### Alott time series
year1 = 1850
year2 = 2017
years = np.arange(year1,year2+1,1)

### Temperature data sets
datasets = [r'HadCRUT4',r'GISTEMP',r'NOAA',r'Berkeley',r'Cowtan \& Way']

### Read in data
filename = 'annual_globaltemps.txt'
data = np.genfromtxt(directorydata + filename,unpack=True,delimiter=',',
                     skip_header=1)

###############################################################################
###############################################################################
###############################################################################    
### Plot figure
### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])

################################################################################  
################################################################################
################################################################################
#### Basic Graph #1
#matplotlib.rc('savefig', facecolor='black')
#matplotlib.rc('axes', edgecolor='w')
#matplotlib.rc('xtick', color='w')
#matplotlib.rc('ytick', color='w')
#matplotlib.rc('axes', labelcolor='w')
#matplotlib.rc('axes', facecolor='black')        
#        
#fig = plt.figure()
#ax = plt.subplot(111) 
#
#color=iter(plt.cm.jet(np.linspace(0.,1,len(data))))
#for i in range(data.shape[0]):
#    cma=next(color)
#    plt.plot(years,data[i,:],alpha=1,linewidth=1,label=datasets[i],
#             color=cma)
#    
#plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
#plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
#plt.xlim([1850,2018])
#plt.ylim([-1,1])
#
#l = plt.legend(fontsize=6.5,loc='upper left')
#for text in l.get_texts():
#    text.set_color('w') 
#
#plt.ylabel('Anomaly (degrees C)')
#plt.title('Annual Mean Global Temperature',color='w')
#plt.savefig(directoryfigure + 'annual_globaltemps_01.png',dpi=300)
#
################################################################################  
################################################################################
################################################################################
#### Basic Graph #2
#matplotlib.rc('savefig', facecolor='black')
#matplotlib.rc('axes', edgecolor='w')
#matplotlib.rc('xtick', color='w')
#matplotlib.rc('ytick', color='w')
#matplotlib.rc('axes', labelcolor='w')
#matplotlib.rc('axes', facecolor='black')        
#        
#fig = plt.figure()
#ax = plt.subplot(111) 
#
#color=iter(plt.cm.viridis(np.linspace(0.,1,len(data))))
#for i in range(data.shape[0]):
#    cma=next(color)
#    plt.plot(years,data[i,:],alpha=1,linewidth=1.5,label=datasets[i],
#             color=cma)
#    
#plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
#plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
#plt.xlim([1850,2018])
#plt.ylim([-1,1])
#
#l = plt.legend(fontsize=6.5,loc='upper left')
#for text in l.get_texts():
#    text.set_color('w') 
#
#plt.ylabel('Anomaly (degrees C)')
#plt.title('Annual Mean Global Temperature',color='w')
#plt.savefig(directoryfigure + 'annual_globaltemps_02.png',dpi=300)
#
################################################################################  
################################################################################
################################################################################
#### Basic Graph #3
#matplotlib.rc('savefig', facecolor='black')
#matplotlib.rc('axes', edgecolor='w')
#matplotlib.rc('xtick', color='w')
#matplotlib.rc('ytick', color='w')
#matplotlib.rc('axes', labelcolor='w')
#matplotlib.rc('axes', facecolor='black')    
#        
#fig = plt.figure()
#ax = plt.subplot(111) 
#
#adjust_spines(ax, ['left', 'bottom'])            
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_linewidth(2)
#ax.spines['left'].set_linewidth(2)
#ax.tick_params('both',length=5.5,width=2,which='major',direction='out')    
#
#color=iter(cmocean.cm.thermal(np.linspace(0.15,1,len(data))))
#for i in range(data.shape[0]):
#    cma=next(color)
#    plt.plot(years,data[i,:],alpha=1,linewidth=1.5,label=datasets[i],
#             color=cma)
#    
#plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
#plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
#plt.xlim([1850,2018])
#plt.ylim([-1,1])
#
#l = plt.legend(fontsize=6.5,loc='upper left')
#for text in l.get_texts():
#    text.set_color('w') 
#
#plt.ylabel('Anomaly (degrees C)')
#plt.title('Annual Mean Global Temperature',color='w',fontsize=15)
#plt.savefig(directoryfigure + 'annual_globaltemps_03.png',dpi=300)

###############################################################################  
###############################################################################
###############################################################################
### Basic Graph #4
matplotlib.rc('savefig', facecolor='black')
matplotlib.rc('axes', edgecolor='w')
matplotlib.rc('xtick', color='w')
matplotlib.rc('ytick', color='w')
matplotlib.rc('axes', labelcolor='w')
matplotlib.rc('axes', facecolor='black')  
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})  
        
fig = plt.figure()
ax = plt.subplot(111) 

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=5.5,width=2,which='major',direction='out')    

color=iter(cmocean.cm.thermal(np.linspace(0.15,1,len(data))))
for i in range(data.shape[0]):
    cma=next(color)
    plt.plot(years,data[i,:],alpha=1,linewidth=1.5,label=datasets[i],
             color=cma)
    
plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
plt.xlim([1850,2018])
plt.ylim([-1,1])

l = plt.legend(fontsize=6.5,loc='upper left')
for text in l.get_texts():
    text.set_color('w') 

plt.title('Annual Mean Global Temperature',color='w',fontsize=15)
plt.ylabel('Anomaly ($^\circ$C)')
plt.savefig(directoryfigure + 'annual_globaltemps_04.png',dpi=300)

###############################################################################  
###############################################################################
###############################################################################
### Basic Graph #5
matplotlib.rc('savefig', facecolor='black')
matplotlib.rc('axes', edgecolor='darkgrey')
matplotlib.rc('xtick', color='darkgrey')
matplotlib.rc('ytick', color='darkgrey')
matplotlib.rc('axes', labelcolor='darkgrey')
matplotlib.rc('axes', facecolor='black')  
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})  
        
fig = plt.figure()
ax = plt.subplot(111) 

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=5.5,width=2,which='major',direction='out')    

color=iter(cmocean.cm.thermal(np.linspace(0.15,1,len(data))))
for i in range(data.shape[0]):
    cma=next(color)
    plt.plot(years,data[i,:],alpha=1,linewidth=1.5,label=datasets[i],
             color=cma)
    
plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
plt.xlim([1850,2018])
plt.ylim([-1,1])

l = plt.legend(fontsize=6.5,loc='upper left')
for text in l.get_texts():
    text.set_color('darkgrey') 

plt.title('ANNUAL MEAN GLOBAL TEMPERATURE',color='darkgrey',fontsize=19)
plt.ylabel('Anomaly ($^\circ$C)',color='darkgrey')
plt.savefig(directoryfigure + 'annual_globaltemps_05.png',dpi=300)

###############################################################################  
###############################################################################
###############################################################################
### Basic Graph #6
matplotlib.rc('savefig', facecolor='black')
matplotlib.rc('axes', edgecolor='darkgrey')
matplotlib.rc('xtick', color='darkgrey')
matplotlib.rc('ytick', color='darkgrey')
matplotlib.rc('axes', labelcolor='darkgrey')
matplotlib.rc('axes', facecolor='black')  
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})  
        
fig = plt.figure()
ax = plt.subplot(111) 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})

adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=5.5,width=2,which='major',direction='out')

color=iter(cmocean.cm.thermal(np.linspace(0.15,1,len(data))))
for i in range(data.shape[0]):
    cma=next(color)
    plt.plot(years,data[i,:],color=cma,alpha=1,linewidth=1.5,
             label=datasets[i])
    
plt.xticks(np.arange(1850,2025,20),np.arange(1850,2025,20),rotation=0,fontsize=9)
plt.yticks(np.arange(-1,1.1,0.5),np.arange(-1,1.1,0.5),rotation=0,fontsize=9)
plt.xlim([1850,2018])
plt.ylim([-1,1])

plt.title(r'\textbf{ANNUAL MEAN GLOBAL TEMPERATURE}',color='darkgrey',fontsize=19)
plt.ylabel(r'Anomaly ($^\circ$C)',color='darkgrey',fontsize=11)
l = plt.legend(shadow=False,fontsize=6.5,loc='lower center',
           bbox_to_anchor=(0.5, -0.025),fancybox=True,ncol=5,
            frameon=False)
for text in l.get_texts():
    text.set_color('darkgrey') 

plt.savefig(directoryfigure + 'annual_globaltemps_06.png',dpi=300)
