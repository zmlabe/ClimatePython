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
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
from scipy.interpolate import griddata as g
import datetime
import numpy.ma as ma

### Define directories
directorydata = '/surtsey/zlabe/seaice/sic/sic_rename/'    
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/VisualsNOAA/' 

yearmin = 1979               # first time includes 12 months
yearmax = 2019
years = np.arange(yearmin,yearmax+1,1)
years2 = np.arange(yearmin,yearmax,1)
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
nsidc = np.genfromtxt('SeaIceIndex_09.txt',unpack=True,usecols=4,
                      delimiter=',',skip_header=1)

def sicRead(directory,years):
    """
    Reads binary sic data
    """
    
    months = [9]
    
    ### Read binary lat x lon arrays
    lons = np.fromfile(directory + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directory + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    ### Read binary sea ice concentration
    ice = np.empty((len(years),136192))
    for i in range(len(years)):    
        if months[0] < 10:
            filename = 'sic_%s_0%s.bin' % (years[i],months[0])
        else:
            filename = 'sic_%s_%s.bin' % (years[i],months[0])
        infile = directory + filename
        
        with open(infile, 'rb') as fr:
            hdr = fr.read(300)
            ice[i,:] = np.fromfile(fr, dtype=np.uint8)
    
    ice = np.reshape(ice,(ice.shape[0],448,304))
    ice = ice/250               # Scale Factor
    
    ### Assign mask
    mask = np.fromfile(directory + 'gsfc_25n.msk',dtype='int8')
    mask = np.reshape(mask,(448,304))
    mask = mask.astype(float)
    mask[np.where(mask == 1.0)] = np.nan
    mask[np.where(mask == 0.0)] = 1.0
    ice = ice*mask
    
    print('Completed: Read SIC data!' )
    return lats,lons,ice,hdr,mask

def read2019(directory):
    directory19 = '/seley/zlabe/seaice/nsidc/sep_2019/'
    days = np.arange(1,30+1,1)
    
    ### Read binary lat x lon arrays
    lons = np.fromfile(directory + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directory + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    ### Read binary sea ice concentration
    ice = np.empty((len(days),136192))
    for i in range(len(days)):    
        if days[i] < 10:
            filename = 'nt_2019090%s_f18_nrt_n.bin' % (days[i])
        else:
            filename = 'nt_201909%s_f18_nrt_n.bin' % (days[i])
        infile = directory19 + filename
        
        with open(infile, 'rb') as fr:
            hdr = fr.read(300)
            ice[i,:] = np.fromfile(fr, dtype=np.uint8)
    
    ice = np.reshape(ice,(ice.shape[0],448,304))
    ice = ice/250               # Scale Factor
    
    ### Take monthly average
    icem = np.mean(ice,axis=0)
    
    ### Assign mask
    mask = np.fromfile(directory + 'gsfc_25n.msk',dtype='int8')
    mask = np.reshape(mask,(448,304))
    mask = mask.astype(float)
    mask[np.where(mask == 1.0)] = np.nan
    mask[np.where(mask == 0.0)] = 1.0
    icem = icem*mask
    
    ### Add axis to append
    icem = icem[np.newaxis,:,:]
    
    return lats,lons,icem,hdr,mask

### Call data
lat2,lon2,sieall,hdr,mask = sicRead(directorydata,years2)
lat19,lon19,sie19,hdr,mask = read2019(directorydata)

### Append data
sie = np.append(sieall,sie19,axis=0)

def calculateThresh(sie,years,thresh):
    """
    Calculates extents depending on SIC threshold
    """
    
    ### Function to calculate SIE for various SIC
    ### Extent is a binary 0 or 1 for SIC threshold
    sieq = sie.copy()
    sieq[np.where(sieq<thresh)]=np.nan
    sieq[np.where(sieq>=thresh)]=1
    
    ### Calculate sea ice extent
    ext = np.zeros((years.shape[0]))
    valyr = np.zeros((sieq.shape))
    for yr in range(years.shape[0]):
        for i in range(lat2.shape[0]):
            for j in range(lon2.shape[1]):
                if sieq[yr,i,j] == 1.0:
                   ### Area of 25 km grid cell; 625 = 25 x 25]
                   valyr[yr,i,j] = 625 # km
        ext[yr] = np.nansum(valyr[yr,:,:])/1e6
        
    return ext

### Call thresholds
ext15 = calculateThresh(sie,years,0.15)
ext30 = calculateThresh(sie,years,0.30)
ext60 = calculateThresh(sie,years,0.60)
ext75 = calculateThresh(sie,years,0.75)
ext90 = calculateThresh(sie,years,0.90)

np.savetxt(directoryfigure + 'sie_tresh_15.txt',ext15,delimiter=',')
np.savetxt(directoryfigure + 'sie_tresh_30.txt',ext30,delimiter=',')
np.savetxt(directoryfigure + 'sie_tresh_60.txt',ext60,delimiter=',')
np.savetxt(directoryfigure + 'sie_tresh_75.txt',ext75,delimiter=',')
np.savetxt(directoryfigure + 'sie_tresh_90.txt',ext90,delimiter=',')

### Testing data (do not use)
fig = plt.figure()
plt.plot(ext15)
plt.plot(nsidc)
plt.savefig(directoryfigure + 'test.png',dpi=300)
fig = plt.figure()
plt.plot(ext15)
plt.plot(ext30)
plt.plot(ext60)
plt.plot(ext75)
plt.plot(ext90)
plt.savefig(directoryfigure + 'test2.png',dpi=300)
