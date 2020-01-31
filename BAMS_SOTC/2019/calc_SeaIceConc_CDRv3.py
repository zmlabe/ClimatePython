
def readSIC(monthindex,oisst):
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    import math
    
    ### Define directories
    directorydata = '/surtsey/zlabe/seaice/CDRv3/monthly/' 
    directorydata2 = '/home/zlabe/Documents/Projects/Tests/oisstv2/Data/'
    
    ### Define attributes
    years = np.arange(1979,2018+1,1)
    yearsall = np.arange(1979,2019+1,1)
    months = np.arange(1,12+1,1)
    monthq = [str(item).zfill(2) for item in months]
    satellite1 = np.repeat(['n07'],103)
    satellite2 = np.repeat(['f08'],53)
    satellite3 = np.repeat(['f11'],45)
    satellite4 = np.repeat(['f13'],147)
    satellite5 = np.repeat(['f17'],132)
    satall = (satellite1,satellite2,satellite3,satellite4,satellite5)
    sat = np.concatenate(satall)
    
    ### Read in all data into [years,months,lat,lon]
    count = 0
    sic = np.empty((years.shape[0],months.shape[0],448,304))
    for i in range(sic.shape[0]):
        for j in range(sic.shape[1]):
            count += 1
            time = '%s%s' % (years[i],monthq[j])
            if time == '198801':
                none = np.empty((448,304))
                none.fill(np.nan)
                sic[i,j,:,:] = none
            else:
                filename = 'seaice_conc_monthly_nh_%s_%s_v03r01.nc' % (sat[count-1],
                                                                         time)
                data = Dataset(directorydata + filename)
                sic[i,j,:,:] = data.variables['goddard_merged_seaice_conc_monthly'][:]
                lat = data.variables['latitude'][:]
                lon = data.variables['longitude'][:]
                data.close()
                
    ### Mask missing data
    sic[np.where(sic == 251.)] = 1
    sic[np.where(sic > 1)] = np.nan
    sicq = sic*100.
    
    ### Select month
    sicm = sicq[:,monthindex,:,:]
    
    ### Read in 2019 data
    datasic = Dataset(directorydata2 + 'seaice_conc_monthly_icdr_nh_f18_20190%s_v01r00.nc' % (monthindex+1))
    sic19 = datasic.variables['seaice_conc_monthly_cdr'][:]
    datasic.close()

    ### Convert to %
    sic19[np.where(np.asarray(sic19) == 251.)] = 1
    sicq19 = sic19[:,:,:] * 100.
    
    ### Compile data sets for CDR
    sicall = np.append(sicm,sicq19,axis=0)
    
    ### oisstv2 data length
    if oisst == True:
        yearsnew = np.arange(1982,2019+1,1)
        yearq = np.where((yearsall >= 1982) & (yearsall <= 2019))[0]
        sicalloi = sicall[yearq,:,:]
    else:
        yearsnew = yearsall
        sicalloi = sicall 
        
    
    return sicalloi,yearsnew,lat,lon

### Test function
import numpy as np
import matplotlib.pyplot as plt
sic,years,lat,lon = readSIC(7,True)