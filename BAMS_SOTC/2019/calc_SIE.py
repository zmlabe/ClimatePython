def readNSIDCData(currentmnq):
    ### Import modules
    import numpy as np
    import urllib.request
    
    url = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/monthly/data/'
    
    if currentmnq < 10:
        monthz = 'N_0%s_extent_v3.0.csv' % currentmnq
    else:
        monthz = 'N_%s_extent_v3.0.csv' % currentmnq
    urls = url + monthz
    
    ### Read files
    raw_data = urllib.request.urlopen(urls)
    dataset = np.genfromtxt(raw_data, skip_header=1,delimiter=',',
                            usecols=[0,1,4,5])
    
    print('\nCompleted: Read NSDIC Sea Ice Index v3 data!')                        
    
    ### Set missing data to nan
    dataset[np.where(dataset==-9999)] = np.nan
    
    ### Variables
    extent = dataset[3:,2]
    return extent

def calcSIE(thresh,region,monq):
    ### Import modules
    import numpy as np
    import calc_SeaIceConc_CDRv3 as ICE
    
    ### Read in sea ice concentration data
    sic,years,lat,lon = ICE.readSIC(monq-1,True)
    lon36 = np.mod(lon,360)
    
    ## Read in land mask
    directorymask = '/surtsey/zlabe/seaice/sic/sic_rename/'
    mask = np.fromfile(directorymask + 'gsfc_25n.msk',dtype='int8')
    mask = np.reshape(mask,(448,304))
    mask = mask.astype(float)
    mask[np.where(mask == 1.0)] = np.nan
    mask[np.where(mask == 0.0)] = 1.0
    sicmask = sic * mask
    
    ### Check for missing data
    sicmask[np.where(sicmask > 100)] = np.nan
    sicmask[np.where(sicmask < thresh)] = np.nan
    
    ### Create binary
    sie = sicmask.copy()
    sie[np.where((sie >= 15.) & (sie <= 100.))] = 1.
    
    if region == 'arctic':
        sieregion = np.where(lat > 50,sie,np.nan) 
    elif region == 'north_barents':
        sieregion1 = np.where(((lat >= 76.4) & (lat <= 79.4)),sie,np.nan)
        sieregion = np.where(((lon36 >= 38) & (lon36 <= 60)),sieregion1,np.nan)
    elif region == 'chukchi':
        sieregion1 = np.where(((lat >= 68) & (lat <= 74)),sie,np.nan)
        sieregion = np.where(((lon36 >= 180) & (lon36 <= 200)),sieregion1,np.nan)
    elif region == 'bering':
        sieregion1 = np.where(((lat >= 54) & (lat <= 64)),sie,np.nan)
        sieregion = np.where(((lon36 >= 180) & (lon36 <= 200)),sieregion1,np.nan)
    elif region == 'global':
        sieregion = sie[:,:,:]
        
    ## Calculate sea ice extent
    area = np.fromfile(directorymask + 'psn25area_v3.dat',dtype='<i4')
    area = np.reshape(area,(448,304))/1000 # scale factor
    area = area.astype(float)
        
    #ext = np.zeros((years.shape[0]))
    #valyr = np.zeros((sieregion.shape))
    #for yr in range(years.shape[0]):
    #    for i in range(lat.shape[0]):
    #        for j in range(lon.shape[1]):
    #            if sieregion[yr,i,j] == 1.0:
    #               ### Area of 25 km grid cell; 625 = 25 x 25]
    #               valyr[yr,i,j] = 625 # km
    #    ext[yr] = np.nansum(valyr[yr,:,:])/1e6
        
    ext = np.zeros((years.shape[0]))
    valyr = np.zeros((sieregion.shape))
    for yr in range(years.shape[0]):
        for i in range(lat.shape[0]):
            for j in range(lon.shape[1]):
                if sieregion[yr,i,j] == 1.0:
                   ### Area of 25 km grid cell; 625 = 25 x 25]
                   valyr[yr,i,j] = area[i,j]
        ext[yr] = np.nansum(valyr[yr,:,:])/1e6
    
    ### Compare to nsidc data
    if region == 'arctic':
        nsidc = readNSIDCData(monq)
        corr = np.corrcoef(ext,nsidc)[0][1]
        print('\n---Correlation is %s!---\n' % np.round(corr,3))
    else:
        nsidc = np.array([np.nan]*years.shape[0])
        
    return ext,nsidc

### Call functions
#import numpy as np
#import matplotlib.pyplot as plt
#ext,nsidc = calcSIE(15,'north_barents',8)
#plt.plot(ext,label='weight')
#plt.plot(nsidc,label='NSIDC')
#plt.legend()
#plt.savefig('/home/zlabe/Desktop/test.png',dpi=300)