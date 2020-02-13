"""
Script plot time series for Barents sea ice and SST
 
Source : Reynolds et al. 2002 [Journal of Climate]
Data: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
Author : Zachary Labe
Date : 4 February 2019
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sts
import calc_SIE as ICEALL
import calc_RegionalSIE as ICEREG

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Data/'
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Figures/'
years = np.arange(1982,2019+1,1)
month = 'aug'
monq = 8

### Read in SST data
ba = np.genfromtxt(directorydata + 'SST_Barents_%s_1982-2019.txt' % month)

### Calculate trend lines
slopeb,interceptb,rvalueb,pvalueb,std_errb = sts.linregress(years,ba)
lineb = slopeb*years + interceptb
confb = std_errb*1.96

### Read in sea ice data
#extb,nsidc = ICEALL.calcSIE(15,'north_barents',monq)
extb = ICEREG.read_RegionalSIE('Barents',monq)

###############################################################################
###############################################################################
###############################################################################
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
        
fig = plt.figure()
ax = plt.subplot(111)
        
adjust_spines(ax, ['left','bottom'])            
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none') 
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2) 
ax.tick_params('both',length=5.5,width=2,which='major',color='dimgrey')     

plt.plot(years,ba,linewidth=4,color='steelblue',marker='o',
         markersize=6,clip_on=False,
         label = r'\textbf{Northern Barents: %s $\bf{\pm}$ %s $\bf{^{\circ}}$C year$\bf{^{-1}}$}' % \
         (np.round(slopeb,2),np.round(confb,2)),zorder=11)
plt.plot(years,lineb,linewidth=2,color='k',linestyle='--',
         dashes=(1,0.4),clip_on=False,zorder=12)

l = plt.legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,0.06),
               ncol=1,frameon=False,fancybox=True,shadow=False)

plt.xticks(np.arange(1982,2020,3),map(str,np.arange(1982,2020,3)),
           size=7)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),
           size=7)
plt.xlim([1982,2019])
plt.ylim([-3,3])
plt.ylabel(r'\textbf{Mean SST Anomaly [$^\circ$C]}',fontsize=11,color='k')
plt.text(1980.3,3.20,r'\textbf{[a]}',fontsize=12,color='dimgrey')
        
################################################################################
a = plt.axes([.48, .67, .46, .2])

#adjust_spines(a, ['left', 'bottom'])
a.spines['top'].set_color('none')
a.spines['right'].set_color('w')
a.spines['left'].set_color('dimgrey')
a.spines['bottom'].set_color('dimgrey')
a.spines['left'].set_linewidth(2)
a.spines['right'].set_linewidth(0)
a.spines['bottom'].set_linewidth(2)
a.tick_params('both',length=4,width=2,which='major',color='dimgrey')

rects = plt.bar(np.arange(0.5,38.5,1),extb,color='navy',
                edgecolor='navy',width=0.6,clip_on=False)

ylabels = map(str,np.arange(0,2.3,0.2))
plt.xticks(np.arange(0.5,46.5,3),map(str,np.arange(1982,2025,3)),
           size=5)
plt.yticks(np.arange(0,2.1,0.05),map(str,np.round(np.arange(0,2.05,0.05),2)),
           size=5)
plt.ylim([0,0.2])
plt.xlim([0,37.7])
a.yaxis.grid(zorder=1,color='dimgrey',alpha=0.35)
a.tick_params(axis='x', which='major', pad=1)
a.tick_params(axis='y', which='major', pad=1)
plt.title(r'\textbf{Barents: Mean Sea Ice Extent [$\bf{\times}$10$\bf{^{6}}$\ \textbf{km}$\bf{^2}$]}',
           fontsize=7,alpha=1,color='k') 
plt.text(-1.7,0.23,r'\textbf{[b]}',fontsize=12,color='dimgrey')

plt.savefig(directoryfigure + 'SST_TimeSeries_%s_Barents.png' % month,dpi = 600)