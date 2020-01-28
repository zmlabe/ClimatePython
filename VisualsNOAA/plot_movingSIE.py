"""
Script reads Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS
Passive Microwave Data, Version 1 binary files for select variables and 
regrids according to selected grid style (e.g., NSIDC EASE grid data). 
 
Source : https://nsidc.org/data/nsidc-0051#
Author : Zachary Labe
Date : 19 August 2018
"""

### Import Modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import datetime
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/seaice/sic/sic_rename/'    
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/VisualsNOAA/' 

yearmin = 1979               # first time includes 12 months
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
months = np.arange(1,13,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print('\n' 'Satellite SIC Read & Regrid - %s' '\n' % titletime) 

### Read in NSIDC official SIE data
sie15 = np.genfromtxt('sie_tresh_15.txt',unpack=True,delimiter=',')
sie30 = np.genfromtxt('sie_tresh_30.txt',unpack=True,delimiter=',')
sie60 = np.genfromtxt('sie_tresh_60.txt',unpack=True,delimiter=',')
sie75 = np.genfromtxt('sie_tresh_75.txt',unpack=True,delimiter=',')
sie90 = np.genfromtxt('sie_tresh_90.txt',unpack=True,delimiter=',')

### Plot figure
fig = plt.figure()
ax = plt.subplot(111) 

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})
matplotlib.rc('savefig', facecolor='black')
matplotlib.rc('axes', edgecolor='darkgrey')
matplotlib.rc('xtick', color='darkgrey')
matplotlib.rc('ytick', color='darkgrey')
matplotlib.rc('axes', labelcolor='darkgrey')
matplotlib.rc('axes', facecolor='black')

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
        
adjust_spines(ax, ['left', 'bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=5.5,width=2,which='major')

bar1, = plt.plot(years,sie15,color=cmocean.cm.ice(0.99))
bar2, = plt.plot(years,sie30,color=cmocean.cm.ice(0.8))
bar3, = plt.plot(years,sie60,color=cmocean.cm.ice(0.65))
bar4, = plt.plot(years,sie75,color=cmocean.cm.ice(0.45))
bar5, = plt.plot(years,sie90,color=cmocean.cm.ice(0.2))
bar11 = plt.fill_between(years,sie15,0,color=cmocean.cm.ice(0.99))
bar22 = plt.fill_between(years,sie30,0,color=cmocean.cm.ice(0.8))
bar33 = plt.fill_between(years,sie60,0,color=cmocean.cm.ice(0.65))
bar44 = plt.fill_between(years,sie75,0,color=cmocean.cm.ice(0.45))
bar55 = plt.fill_between(years,sie90,0,color=cmocean.cm.ice(0.2))

plt.text(2019.4,sie15[-1]-0.15,r'\textbf{15\%}',color=cmocean.cm.ice(0.99),fontsize=15)
plt.text(2019.4,sie30[-1]-0.15,r'\textbf{30\%}',color=cmocean.cm.ice(0.8),fontsize=15)
plt.text(2019.4,sie60[-1]-0.15,r'\textbf{60\%}',color=cmocean.cm.ice(0.65),fontsize=15)
plt.text(2019.4,sie75[-1]-0.15,r'\textbf{75\%}',color=cmocean.cm.ice(0.45),fontsize=15)
plt.text(2019.4,sie90[-1]-0.15,r'\textbf{90\%}',color=cmocean.cm.ice(0.2),fontsize=15)

plt.text(2019,5.5,r'\textbf{Sea Ice}',color='darkgrey',fontsize=8,
          ha='center',va='center')
plt.text(2019.,5.2,r'\textbf{Concentration}',color='darkgrey',fontsize=8,
          ha='center',va='center')

plt.ylabel(r'\textbf{Extent [$\bf{\times}$10$^{6}$\ \textbf{km}$^2$]}',
           fontsize=13,color='darkgrey',labelpad=10)
plt.title(r'\textbf{SEPTEMBER -- ARCTIC SEA ICE EXTENT}',fontsize=20,
          color='w')

plt.xticks(np.arange(1979,2020,4),map(str,np.arange(1979,2020,4)),rotation=0,fontsize=10)
ylabels = map(str,np.arange(0,10,1))
plt.yticks(np.arange(0,10,1),ylabels,fontsize=10)
plt.ylim([0,8])
plt.xlim([1979,2019])

#def animate(i):
#    x = years[0:(i)]
#    y15 = sie15[0:(i)]
#    y30 = sie30[0:(i)]
#    y60 = sie60[0:(i)]
#    y75 = sie75[0:(i)]
#    y90 = sie90[0:(i)]
#    p15 = plt.fill_between(x,y15,0,color=cmocean.cm.ice(0.99))
#    p30 = plt.fill_between(x,y30,0,color=cmocean.cm.ice(0.8))
#    p60 = plt.fill_between(x,y60,0,color=cmocean.cm.ice(0.65))
#    p75 = plt.fill_between(x,y75,0,color=cmocean.cm.ice(0.45))
#    p90 = plt.fill_between(x,y90,0,color=cmocean.cm.ice(0.2))
#    return p15,p30,p60,p75,p90
#    
#ani = animation.FuncAnimation(fig,animate,frames=70,interval=180,repeat=True,
#                              blit=True)
#ani.save('NSIDC_SIC_thresh_moving.gif',writer='imagemagick',dpi=300)

plt.savefig('NSIDC_SIC_thresh_moving.png',dpi=300)