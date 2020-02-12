"""
Read in regional mask data (https://nsidc.org/sites/nsidc.org/files/files/data/noaa/g02135/Sea-Ice-Analysis-Spreadsheets-Overview.pdf)
"""

import matplotlib.pyplot as plt
import numpy as np

### Define constants
directorydata = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Data/'
directorydata2 = '/surtsey/zlabe/seaice/sic/sic_rename/'    
directoryfigure = '/home/zlabe/Documents/Projects/ClimatePython/BAMS_SOTC/2019/Figures/'

### Read in mask 
filename = 'Arctic_region_mask_Meier_AnnGlaciol2007.msk'
mask = np.fromfile(directorydata + filename,
                            dtype=np.uint8).reshape((448,304))

### Read in lat/lon
lons = np.fromfile(directorydata2  + 'psn25lons_v3.dat',dtype='<i4')
lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
lats = np.fromfile(directorydata2  + 'psn25lats_v3.dat',dtype='<i4')
lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor

### Find Barents (#8)
arg = np.where(mask.ravel() == 8)[0]
latq = lats.ravel()[arg]
lonq = lons.ravel()[arg]

### Save polygon points
np.savetxt(directorydata + 'Lats_Barents.txt',latq)
np.savetxt(directorydata + 'Lons_Barents.txt',lonq)