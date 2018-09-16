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
directorydata = '/seley/zlabe/seaice/nsidc/SepMinSIC/'    
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/VisualsNOAA/' 

yearmin = 1979               # first time includes 12 months
yearmax = 2018
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

def sicRead(directory,years):
    """
    Reads binary sic data
    """
    
    ### Read binary lat x lon arrays
    directorydata1 = '/surtsey/zlabe/seaice/sic/sic_rename/' 
    lons = np.fromfile(directorydata1 + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directorydata1 + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    ### Read binary sea ice concentration
    ice = np.empty((len(years),136192))
    for i in range(len(years)):    
        if years[i] <= 1986:
            filename = 'nt_%s0915_n07_v1.1_n.bin' % (years[i])
        elif years[i]>1986 and years[i]<=1991:
            filename = 'nt_%s0915_f08_v1.1_n.bin' % (years[i])
        elif years[i]>1991 and years[i]<=1995:
            filename = 'nt_%s0915_f11_v1.1_n.bin' % (years[i])
        elif years[i]>1995 and years[i]<=2007:
            filename = 'nt_%s0915_f13_v1.1_n.bin' % (years[i])
        elif years[i]>2007 and years[i]<=2018:
            filename = 'nt_%s0915_f17_v1.1_n.bin' % (years[i])
        infile = directory + filename
        
        with open(infile, 'rb') as fr:
            hdr = fr.read(300)
            ice[i,:] = np.fromfile(fr, dtype=np.uint8)
    
    ice = np.reshape(ice,(ice.shape[0],448,304))
    ice = ice/250               # Scale Factor
    
    ### Assign mask
    mask = np.fromfile(directorydata1 + 'gsfc_25n.msk',dtype='int8')
    mask = np.reshape(mask,(448,304))
    mask = mask.astype(float)
    mask[np.where(mask == 1.0)] = np.nan
    mask[np.where(mask == 0.0)] = 1.0
    ice = ice*mask
    
    ice = ice*100.
    
    print('Completed: Read SIC data!' )
    return lats,lons,ice,hdr,mask

### Call data
lat2,lon2,sic,hdr,mask = sicRead(directorydata,years)

### Calculate climatology
yearmean = np.where((years>=1981) & (years<=2010))[0]
mean = np.nanmedian(sic[yearmean],axis=0)

def netcdfFile(lats,lons,var,directory):
    print('\n>>> Using netcdfFile function!')
    
    name = 'NSIDC_09_1979-2018.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = '15 September 1979-2018 SIC data'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable('sic','f4',('years','lat','lon'))
    
    ### Units
    varns.units = '%'
    ncfile.title = 'Daily NSIDC Data'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NSIDC Data'
    ncfile.references = 'SSMIS (DMSP)]'
    
    ### Data
    years[:] = np.arange(1979,2018+1,1)
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print('*Completed: Created netCDF4 File!')
    
def netcdfFileMean(lats,lons,var,directory):
    print('\n>>> Using netcdfFile function!')
    
    name = 'NSIDC_09_Median.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = '15 September Median (1981-2010) SIC data'
    
    ### Dimensions
    ncfile.createDimension('lat',var.shape[0])
    ncfile.createDimension('lon',var.shape[1])
    
    ### Variables
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable('sic','f4',('lat','lon'))
    
    ### Units
    varns.units = '%'
    ncfile.title = 'Daily NSIDC Data Climatology'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NSIDC Data'
    ncfile.references = 'SSMIS (DMSP)]'
    
    ### Data
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print('*Completed: Created netCDF4 File!')
    
netcdfFile(lat2,lon2,sic,directorydata)
netcdfFileMean(lat2,lon2,mean,directorydata)