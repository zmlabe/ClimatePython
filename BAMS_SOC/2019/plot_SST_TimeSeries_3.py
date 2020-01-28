"""
Script plot time series for polar cap, barents, and chukchi data
 
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
month = 'aug'

### Read in data
polar = np.genfromtxt(directorydata + 'SST_Arctic_%s_1982-2019.txt' % month)
ch = np.genfromtxt(directorydata + 'SST_Chukchi_%s_1982-2019.txt' % month)
ba = np.genfromtxt(directorydata + 'SST_Barents_%s_1982-2019.txt' % month)

### Calculate trend lines
slopep,interceptp,rvaluep,pvaluep,std_errp = sts.linregress(years,polar)
linep = slopep*years + interceptp
confp = std_errp*1.96

slopec,interceptc,rvaluec,pvaluec,std_errc = sts.linregress(years,ch)
linec = slopec*years + interceptc
confc = std_errc*1.96

slopeb,interceptb,rvalueb,pvalueb,std_errb = sts.linregress(years,ba)
lineb = slopeb*years + interceptb
confb = std_errb*1.96

###############################################################################
###############################################################################
############################################midnightblue###################################
### Create subplots SST data
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

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
        
###############################################################################
fig = plt.figure(figsize=(9,5))
gs = fig.add_gridspec(2,2)

ax1 = fig.add_subplot(gs[:,0])
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('w')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(0)
ax1.spines['bottom'].set_linewidth(2)
ax1.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.plot(years,polar,linewidth=3.8,color='k',marker='o',
         markersize=6,clip_on=False,
         label = r'\textbf{Arctic [$\bf{>}$67$\bf{^{\circ}}$N]: %s $\bf{\pm}$ %s $\bf{^{\circ}}$C/yr}' % \
         (np.round(slopep,2),np.round(confp,2)))
plt.plot(years,linep,linewidth=1.5,color='k',linestyle='--',
         dashes=(1,0.4),clip_on=False)

l = plt.legend(fontsize=12,loc='upper center',bbox_to_anchor=(0.5,0.09),
               ncol=1,frameon=False,fancybox=True,shadow=False)

ax1.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35)
plt.xticks(np.arange(1982,2020,3),map(str,np.arange(1982,2020,3)),
           size=6)
plt.yticks(np.arange(-10,11,0.5),map(str,np.arange(-10,11,0.5)),
           size=6)
plt.xlim([1982,2019])
plt.ylim([-1.5,1.5])
plt.text(1982,1.39,r'\textbf{[a]}',fontsize=12,color='dimgrey')

plt.ylabel(r'\textbf{Mean SST Anomaly [$^\circ$C]}',fontsize=12,color='k')

###############################################################################
ax2 = fig.add_subplot(gs[0,1])
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('w')
ax2.spines['left'].set_color('dimgrey')
ax2.spines['bottom'].set_color('dimgrey')
ax2.spines['left'].set_linewidth(2)
ax2.spines['right'].set_linewidth(0)
ax2.spines['bottom'].set_linewidth(2)
ax2.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.plot(years,ba,linewidth=2.5,color='deepskyblue',marker='o',
         markersize=4,clip_on=False,
         label = r'\textbf{Northern Barents: %s $\bf{\pm}$ %s $\bf{^{\circ}}$C/yr}' % \
         (np.round(slopeb,2),np.round(confb,2)))
plt.plot(years,lineb,linewidth=1.5,color='midnightblue',linestyle='--',
         dashes=(1,0.4),clip_on=False)

l = plt.legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,0.09),
               ncol=1,frameon=False,fancybox=True,shadow=False)

ax2.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35)
plt.xticks(np.arange(1982,2020,3),map(str,np.arange(1982,2020,3)),
           size=6)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),
           size=6)
plt.xlim([1982,2019])
plt.ylim([-3.5,3.5])
plt.text(1982,3.1,r'\textbf{[b]}',fontsize=12,color='dimgrey')

###############################################################################
ax3 = fig.add_subplot(gs[1,1])
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('w')
ax3.spines['left'].set_color('dimgrey')
ax3.spines['bottom'].set_color('dimgrey')
ax3.spines['left'].set_linewidth(2)
ax3.spines['right'].set_linewidth(0)
ax3.spines['bottom'].set_linewidth(2)
ax3.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.plot(years,ch,linewidth=2.5,color='crimson',marker='o',
         markersize=4,clip_on=False,
         label = r'\textbf{Chukchi: %s $\bf{\pm}$ %s $\bf{^{\circ}}$C/yr}' % \
         (np.round(slopec,2),np.round(confc,2)))
plt.plot(years,linec,linewidth=1.5,color='maroon',linestyle='--',
         dashes=(1,0.4),clip_on=False)

l = plt.legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,0.09),
               ncol=1,frameon=False,fancybox=True,shadow=False)

ax3.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35)
plt.xticks(np.arange(1982,2020,3),map(str,np.arange(1982,2020,3)),
           size=6)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),
           size=6)
plt.xlim([1982,2019])
plt.ylim([-3.5,3.5])
plt.text(1982,3.1,r'\textbf{[c]}',fontsize=12,color='dimgrey')

plt.tight_layout()

plt.savefig(directoryfigure + 'SST_TimeSeries_%s_3.png' % month,dpi = 600)