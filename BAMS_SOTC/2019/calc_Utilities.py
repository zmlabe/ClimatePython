"""
Functions are useful untilities for SITperturb experiments
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    [1] calcDecJan(varx,vary,lat,lon,level,levsq)
    [2] calcDecJanFeb(varx,vary,lat,lon,level,levsq)
    [3] calc_FDR_ttest(varx,vary,alpha_f)
    [4] calc_indttest(varx,vary)
    [5] calc_weightedAve(var,lats)
    [6] calc_spatialCorr(varx,vary,lats,lons,weight)
    [7] calc_RMSE(varx,vary,lats,lons,weight)
    [8] calc_spatialCorrHeight(varx,vary,lats,lons,weight)
    [9] calc_spatialCorrHeightLev(varx,vary,lats,lons,weight,levelq)
    [10] detrendData(datavar,years,level,yearmn,yearmx)
    [11] detrendDataR(datavar,years,level,yearmn,yearmx)
    [12] mk_test(x, alpha)
"""

def calcDecJan(varx,vary,lat,lon,level,levsq):
    """
    Function calculates average for December-January

    Parameters
    ----------
    varx : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    vary : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    level : string
        Height of variable (surface or profile)
    levsq : integer
        number of levels
        
    Returns
    -------
    varx_dj : 3d array or 4d array
        [year,lat,lon] or [year,lev,lat,lon]
    vary_dj : 3d array
        [year,lat,lon] or [year,lev,lat,lon]

    Usage
    -----
    varx_dj,vary_dj = calcDecJan(varx,vary,lat,lon,level,levsq)
    """
    print('\n>>> Using calcDecJan function!')
    
    ### Import modules
    import numpy as np
    
    ### Reshape for 3d variables
    if level == 'surface':    
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12),
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12),
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_dj = np.empty((varx.shape[0]-1,lat.shape[0],lon.shape[0]))
        vary_dj = np.empty((vary.shape[0]-1,lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djappendh = np.append(varxravel[11+i,:,:],varxravel[12+i,:,:])
            djappendf = np.append(varyravel[11+i,:,:],varyravel[12+i,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)
    ### Reshape for 4d variables
    elif level == 'profile':
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12.),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12.),levsq,
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_dj = np.empty((int(varx.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        vary_dj = np.empty((int(vary.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djappendh = np.append(varxravel[11+i,:,:,:],
                                  varxravel[12+i,:,:,:])
            djappendf = np.append(varyravel[11+i,:,:,:],
                                  varyravel[12+i,:,:,:])    
            varx_dj[counter,:,:] = np.nanmean(np.reshape(djappendh,
                                    (2,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                   
            vary_dj[counter,:,:] = np.nanmean(np.reshape(djappendf,
                                    (2,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                               
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
                                
    print('Completed: Organized data by months (ON,DJ,FM)!')

    print('*Completed: Finished calcDecJan function!')
    return varx_dj,vary_dj

###############################################################################
###############################################################################
###############################################################################

def calcDecJanFeb(varx,vary,lat,lon,level,levsq):
    """
    Function calculates average for December-January-February

    Parameters
    ----------
    varx : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    vary : 4d array or 5d array
        [year,month,lat,lon] or [year,month,lev,lat,lon]
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    level : string
        Height of variable (surface or profile)
    levsq : integer
        number of levels
        
    Returns
    -------
    varx_djf : 3d array or 4d array
        [year,lat,lon] or [year,lev,lat,lon]
    vary_djf : 3d array
        [year,lat,lon] or [year,lev,lat,lon]

    Usage
    -----
    varx_djf = calcDecJanFeb(varx,vary,lat,lon,level,levsq)
    """
    print('\n>>> Using calcDecJanFeb function!')
    
    ### Import modules
    import numpy as np
    
    ### Reshape for 3d variables
    if level == 'surface':    
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12),
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12),
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_djf = np.empty((varx.shape[0]-1,lat.shape[0],lon.shape[0]))
        vary_djf = np.empty((vary.shape[0]-1,lat.shape[0],lon.shape[0]) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djfappendh1 = np.append(varxravel[11+i,:,:],varxravel[12+i,:,:])
            djfappendf1 = np.append(varyravel[11+i,:,:],varyravel[12+i,:,:])  
            djfappendh = np.append(djfappendh1,varxravel[13+i,:,:])
            djfappendf = np.append(djfappendf1,varyravel[13+i,:,:]) 
            varx_djf[counter,:,:] = np.nanmean(np.reshape(djfappendh,
                                    (3,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)                   
            vary_djf[counter,:,:] = np.nanmean(np.reshape(djfappendf,
                                    (3,int(lat.shape[0]),int(lon.shape[0]))),
                                    axis=0)
    ### Reshape for 4d variables
    elif level == 'profile':
        varxravel = np.reshape(varx.copy(),
                           (int(varx.shape[0]*12.),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        varyravel = np.reshape(vary.copy(),
                               (int(vary.shape[0]*12.),levsq,
                                int(lat.shape[0]),int(lon.shape[0]))) 
                               
        varx_djf = np.empty((int(varx.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])))
        vary_djf = np.empty((int(vary.shape[0]-1),levsq,
                            int(lat.shape[0]),int(lon.shape[0])) )                 
        for i in range(0,varxravel.shape[0]-12,12):
            counter = 0
            if i >= 12:
                counter = i//12
            djfappendh1 = np.append(varxravel[11+i,:,:,:],
                                  varxravel[12+i,:,:,:])
            djfappendf1 = np.append(varyravel[11+i,:,:,:],
                                  varyravel[12+i,:,:,:]) 
            djfappendh = np.append(djfappendh1,
                                  varxravel[13+i,:,:,:])
            djfappendf = np.append(djfappendf1,
                                  varyravel[13+i,:,:,:])  
            varx_djf[counter,:,:] = np.nanmean(np.reshape(djfappendh,
                                    (3,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                   
            vary_djf[counter,:,:] = np.nanmean(np.reshape(djfappendf,
                                    (3,levsq,int(lat.shape[0]),
                                     int(lon.shape[0]))),axis=0)                               
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
                                
    print('Completed: Organized data by months (DJF)!')

    print('*Completed: Finished calcDecJanFeb function!')
    return varx_djf

###############################################################################
###############################################################################
###############################################################################
    
def calc_FDR_ttest(varx,vary,alpha_f):
    """
    Function first calculates statistical difference for 2 independent
    sample t-test and then adjusts using a false discovery rate (FDR)
    where alpha_o = alpha_FDR

    Parameters
    ----------
    varx : 2d or 3d array
    vary : 2d or 3d array
    alpha_f : float (alpha_o = alpha_FDR)
    
    Returns
    -------
    pruns : 1d or 2d array of adjusted p values

    Usage
    -----
    calc_FDR_ttest(varx,vary,alpha_f)
    """
    print('\n>>> Using calc_FDR_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    import statsmodels.stats.multitest as fdr
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_ind(varx,vary,nan_policy='omit')
    
    ### Ravel all 2d pvalues
    if pvalue.ndim == 2:
        pvalall = np.reshape(pvalue,(pvalue.shape[0]*pvalue.shape[1]))
    elif pvalue.ndim == 3:
        pvalall = np.reshape(pvalue,(pvalue.shape[0]*pvalue.shape[1]*pvalue.shape[2]))
    else:
        pvalall = pvalue
    
    ### Calculate false discovery rate
    prunsq = np.empty((pvalall.shape))
    score = np.empty((pvalall.shape))
    prunsq.fill(np.nan)
    score.fill(np.nan)
    
    ### Check for nans before correction!!
    mask = np.isfinite(pvalall[:])
    score[mask],prunsq[mask] = fdr.fdrcorrection(pvalall[mask],alpha=alpha_f,
                                      method='indep')
        
    ### Reshape into correct dimensions
    pruns = np.reshape(prunsq,(pvalue.shape))
    
    ### Mask variables by their adjusted p-values
    pruns[np.where(pruns >= alpha_f)] = np.nan
    pruns[np.where(pruns < alpha_f)] = 1.
    
    print('*Completed: Finished calc_FDR_ttest function!')
    return pruns

###############################################################################
###############################################################################
###############################################################################
    
def calc_indttest(varx,vary):
    """
    Function calculates statistical difference for 2 independent
    sample t-test

    Parameters
    ----------
    varx : 3d array
    vary : 3d array
    
    Returns
    -------
    stat = calculated t-statistic
    pvalue = two-tailed p-value

    Usage
    -----
    stat,pvalue = calc_ttest(varx,vary)
    """
    print('\n>>> Using calc_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_ind(varx,vary,nan_policy='omit')
    
    ### Significant at 95% confidence level
    pvalue[np.where(pvalue >= 0.1)] = np.nan
    pvalue[np.where(pvalue < 0.1)] = 1.
    
    print('*Completed: Finished calc_ttest function!')
    return stat,pvalue

###############################################################################
###############################################################################
###############################################################################

def calc_weightedAve(var,lats):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    
    Parameters
    ----------
    var : 5d,4d,3d,2d array of a gridded variable
    lats : 2d array of latitudes
    
    Returns
    -------
    meanvar : weighted average for 3d,2d,1d array

    Usage
    -----
    meanvar = calc_weightedAve(var,lats)
    """
    print('\n>>> Using calc_weightedAve function!')
    
    ### Import modules
    import numpy as np
    
    ### Calculate weighted average for various dimensional arrays
    if var.ndim == 5:
        meanvar = np.empty((var.shape[0],var.shape[1],var.shape[2]))
        for ens in range(var.shape[0]):
            for i in range(var.shape[1]):
                for j in range(var.shape[2]):
                    varq = var[ens,i,j,:,:]
                    mask = np.isfinite(varq) & np.isfinite(lats)
                    varmask = varq[mask]
                    areamask = np.cos(np.deg2rad(lats[mask]))
                    meanvar[ens,i,j] = np.nansum(varmask*areamask) \
                                        /np.sum(areamask)  
    elif var.ndim == 4:
        meanvar = np.empty((var.shape[0],var.shape[1]))
        for i in range(var.shape[0]):
            for j in range(var.shape[1]):
                varq = var[i,j,:,:]
                mask = np.isfinite(varq) & np.isfinite(lats)
                varmask = varq[mask]
                areamask = np.cos(np.deg2rad(lats[mask]))
                meanvar[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
    elif var.ndim == 3:
        meanvar = np.empty((var.shape[0]))
        for i in range(var.shape[0]):
            varq = var[i,:,:]
            mask = np.isfinite(varq) & np.isfinite(lats)
            varmask = varq[mask]
            areamask = np.cos(np.deg2rad(lats[mask]))
            meanvar[i] = np.nansum(varmask*areamask)/np.sum(areamask)
    elif var.ndim == 2:
        meanvar = np.empty((var.shape[0]))
        varq = var[:,:]
        mask = np.isfinite(varq) & np.isfinite(lats)
        varmask = varq[mask]
        areamask = np.cos(np.deg2rad(lats[mask]))
        meanvar = np.nansum(varmask*areamask)/np.sum(areamask)
    else:
        print(ValueError('Variable has the wrong dimensions!'))
     
    print('Completed: Weighted variable average!')
    
    print('*Completed: Finished calc_weightedAve function!')
    return meanvar

###############################################################################
###############################################################################
###############################################################################
    
def calc_spatialCorr(varx,vary,lats,lons,weight):
    """
    Calculates spatial correlation from pearson correlation coefficient
    
    Parameters
    ----------
    varx : 2d array
    vary : 2d array
    lons : 1d array of latitude
    weight : string (yes or no)
    
    Returns
    -------
    corrcoef : 1d array of correlation coefficient (pearson r)
    
    Usage
    -----
    corrcoef = calc_spatialCorr(varx,vary,lats,lons)
    """
    
    print('\n>>> Using calc_spatialCorr function!')
    ### Import modules
    import numpy as np
    
    if weight == 'yes': # Computed weighted correlation coefficient   
        ### mask 
        mask = 'yes'
        if mask == 'yes':
            latq = np.where(lats > 30)[0]
            lats = lats[latq]
            varx = varx[latq,:]
            vary = vary[latq,:]
            print('MASKING LATITUDES!')
        
        ### Create 2d meshgrid for weights 
        lon2,lat2 = np.meshgrid(lons,lats)
        
        ### Create 2d array of weights based on latitude
        gw = np.cos(np.deg2rad(lat2))
        
        def m(x, w):
            """Weighted Mean"""
    
            wave = np.sum(x * w) / np.sum(w)
            print('Completed: Computed weighted average!') 
            return wave
        
        def cov(x, y, w):
            """Weighted Covariance"""
            
            wcov = np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)
            print('Completed: Computed weighted covariance!')
            return wcov
        
        def corr(x, y, w):
            """Weighted Correlation"""
            
            wcor = cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))
            print('Completed: Computed weighted correlation!')
            return wcor
        
        corrcoef = corr(varx,vary,gw)
        
    elif weight == 'no':   
        ### Correlation coefficient from numpy function (not weighted)
        corrcoef= np.corrcoef(varx.ravel(),vary.ravel())[0][1]
        print('Completed: Computed NON-weighted correlation!')
        
    else:
        ValueError('Wrong weighted arguement in function!')
    
    print('*Completed: Finished calc_SpatialCorr function!')
    return corrcoef

###############################################################################
###############################################################################
###############################################################################

def calc_RMSE(varx,vary,lats,lons,weight):
        """
        Calculates root mean square weighted average
        
        Parameters
        ----------
        varx : 2d array
        vary : 2d array
        lons : 1d array of latitude
        weight : string (yes or no)
        
        Returns
        -------
        rmse : 1d array
        
        Usage
        -----
        rmse = calc_RMSE(varx,vary,lats,lons)
        """
        
        print('\n>>> Using calc_RMSE function!')
        ### Import modules
        import numpy as np
        from sklearn.metrics import mean_squared_error
        
        if weight == 'yes': # Computed weighted correlation coefficient   
            ### mask
            mask = 'yes'
            if mask == 'yes':
                latq = np.where(lats > 30)[0]
                lats = lats[latq]
                varx = varx[latq,:]
                vary = vary[latq,:]
                print('MASKING LATITUDES!')
            
            ### Create 2d meshgrid for weights 
            lon2,lat2 = np.meshgrid(lons,lats)
            
            ### Create 2d array of weights based on latitude
            gw = np.cos(np.deg2rad(lat2))
            
            ### Calculate rmse 
            sq_err = (varx - vary)**2
            rmse = np.sqrt((np.sum(sq_err*gw))/np.sum(gw))
     
        elif weight == 'no':   
            ### Root mean square error from sklearn (not weighted)
            rmse = np.sqrt(mean_squared_error(varx.ravel(),vary.ravel()))
            print('Completed: Computed NON-weighted correlation!')
            
        else:
            ValueError('Wrong weighted arguement in function!')
        
        print('*Completed: Finished calc_RMSE function!')
        return rmse
    
###############################################################################
###############################################################################
###############################################################################
    
def calc_spatialCorrHeight(varx,vary,levs,lons,weight):
    """
    Calculates spatial correlation from pearson correlation coefficient for
    grids over vertical height (17 pressure coordinate levels)
    
    Parameters
    ----------
    varx : 2d array
    vary : 2d array
    levs : 1d array of levels
    lons : 1d array of latitude
    weight : string (yes or no)
    
    Returns
    -------
    corrcoef : 1d array of correlation coefficient (pearson r)
    
    Usage
    -----
    corrcoef = calc_spatialCorrHeight(varx,vary,lats,lons)
    """
    
    print('\n>>> Using calc_spatialCorrHeight function!')
    ### Import modules
    import numpy as np
    
    if weight == 'yes': # Computed weighted correlation coefficient   
        
        ### Create 2d meshgrid for weights 
        lon2,lev2 = np.meshgrid(lons,levs)
        
        ### Create 2d array of weights based on latitude
        gwq = np.array([0.25,0.25,0.25,0.25,0.25,0.25,0.4,0.5,0.5,0.5,
              0.5,0.5,0.5,0.7,0.7,0.7,1.])
        gw,gw2 = np.meshgrid(lons,gwq)
        
        def m(x, w):
            """Weighted Mean"""
    
            wave = np.sum(x * w) / np.sum(w)
            print('Completed: Computed weighted average (17 P Levels)!') 
            return wave
        
        def cov(x, y, w):
            """Weighted Covariance"""
            
            wcov = np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)
            print('Completed: Computed weighted covariance (17 P Levels)!')
            return wcov
        
        def corr(x, y, w):
            """Weighted Correlation"""
            
            wcor = cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))
            print('Completed: Computed weighted correlation (17 P Levels)!')
            return wcor
        
        corrcoef = corr(varx,vary,gw2)
        
    elif weight == 'no':   
        ### Correlation coefficient from numpy function (not weighted)
        corrcoef= np.corrcoef(varx.ravel(),vary.ravel())[0][1]
        print('Completed: Computed NON-weighted correlation!')
        
    else:
        ValueError('Wrong weighted argument in function!')
    
    print('*Completed: Finished calc_SpatialCorrHeight function!')
    return corrcoef

###############################################################################
###############################################################################
###############################################################################
    
def calc_spatialCorrHeightLev(varx,vary,levs,lons,weight,levelq):
    """
    Calculates spatial correlation from pearson correlation coefficient for
    grids over vertical height (17 pressure coordinate levels). Change the 
    weighting for different level correlations
    
    Parameters
    ----------
    varx : 2d array
    vary : 2d array
    levs : 1d array of levels
    lons : 1d array of latitude
    weight : string (yes or no)
    levelq : string (all, tropo, strato)
    
    Returns
    -------
    corrcoef : 1d array of correlation coefficient (pearson r)
    
    Usage
    -----
    corrcoef = calc_spatialCorrHeight(varx,vary,lats,lons,levels)
    """
    
    print('\n>>> Using calc_spatialCorrHeightLev function!')
    ### Import modules
    import numpy as np
    
    if weight == 'yes': # Computed weighted correlation coefficient   
        
        ### Create 2d meshgrid for weights 
        lon2,lev2 = np.meshgrid(lons,levs)
        
        if levelq == 'col':
            ### Create 2d array of weights based on latitude
            gwq = np.array([0.25,0.25,0.25,0.25,0.25,0.25,0.4,0.5,0.5,0.5,
                  0.5,0.5,0.5,0.7,0.7,0.7,1.])
            gw,gw2 = np.meshgrid(lons,gwq)
        elif levelq == 'tropo':
           gwq = np.array([1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.2,0.2,0.,0.,0.,
                           0.,0.,0.,0.,0.])
           gw,gw2 = np.meshgrid(lons,gwq)
        elif levelq == 'strato':
           gwq = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.5,1.,1.,1.,1.
                           ,1.,1.])
           gw,gw2 = np.meshgrid(lons,gwq)
        
        def m(x, w):
            """Weighted Mean"""
    
            wave = np.sum(x * w) / np.sum(w)
            print('Completed: Computed weighted average (17 P Levels)!') 
            return wave
        
        def cov(x, y, w):
            """Weighted Covariance"""
            
            wcov = np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)
            print('Completed: Computed weighted covariance (17 P Levels)!')
            return wcov
        
        def corr(x, y, w):
            """Weighted Correlation"""
            
            wcor = cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))
            print('Completed: Computed weighted correlation (17 P Levels)!')
            return wcor
        
        corrcoef = corr(varx,vary,gw2)
        
    elif weight == 'no':   
        ### Correlation coefficient from numpy function (not weighted)
        corrcoef= np.corrcoef(varx.ravel(),vary.ravel())[0][1]
        print('Completed: Computed NON-weighted correlation!')
        
    else:
        ValueError('Wrong weighted argument in function!')
    
    print('*Completed: Finished calc_SpatialCorrHeightLev function!')
    return corrcoef

###############################################################################
###############################################################################
###############################################################################

def detrendData(datavar,years,level,yearmn,yearmx):
    """
    Function removes linear trend

    Parameters
    ----------
    datavar : 4d numpy array or 5d numpy array 
        [ensemble,year,lat,lon] or [ensemble,year,level,lat,lon]
    years : 1d numpy array
        [years]
    level : string
        Height of variable (surface or profile)
    yearmn : integer
        First year
    yearmx : integer
        Last year
    
    Returns
    -------
    datavardt : 4d numpy array or 5d numpy array 
        [ensemble,year,lat,lon] or [ensemble,year,level,lat,lon]
        

    Usage
    -----
    datavardt = detrendData(datavar,years,level,yearmn,yearmx)
    """
    print('\n>>> Using detrendData function! \n')
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Import modules
    import numpy as np
    import scipy.stats as sts

    ### Slice time period
    sliceq = np.where((years >= yearmn) & (years <= yearmx))[0]
    datavar = datavar[:,sliceq,:,:]
    
    ### Detrend data array
    if level == 'surface':
        x = np.arange(datavar.shape[1])
        
        slopes = np.empty((datavar.shape[0],datavar.shape[2],datavar.shape[3]))
        intercepts = np.empty((datavar.shape[0],datavar.shape[2],
                               datavar.shape[3]))
        for ens in range(datavar.shape[0]):
            print('-- Detrended data for ensemble member -- #%s!' % (ens+1))
            for i in range(datavar.shape[2]):
                for j in range(datavar.shape[3]):
                    mask = np.isfinite(datavar[ens,:,i,j])
                    y = datavar[ens,:,i,j]
                    
                    if np.sum(mask) == y.shape[0]:
                        xx = x
                        yy = y
                    else:
                        xx = x[mask]
                        yy = y[mask]
                    
                    if np.isfinite(np.nanmean(yy)):
                        slopes[ens,i,j],intercepts[ens,i,j], \
                        r_value,p_value,std_err = sts.linregress(xx,yy)
                    else:
                        slopes[ens,i,j] = np.nan
                        intercepts[ens,i,j] = np.nan
        print('Completed: Detrended data for each grid point!')
                                
    print('\n>>> Completed: Finished detrendData function!')
    return slopes

###############################################################################
###############################################################################
###############################################################################

def detrendDataR(datavar,years,level,yearmn,yearmx):
    """
    Function removes linear trend from reanalysis data

    Parameters
    ----------
    datavar : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]
    years : 1d numpy array
        [years]
    level : string
        Height of variable (surface or profile)
    yearmn : integer
        First year
    yearmx : integer
        Last year
    
    Returns
    -------
    datavardt : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]
        
    Usage
    -----
    datavardt = detrendDataR(datavar,years,level,yearmn,yearmx)
    """
    print('\n>>> Using detrendData function! \n')
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### Slice time period
    sliceq = np.where((years >= yearmn) & (years <= yearmx))[0]
    datavar = datavar[sliceq,:,:]
    
    ### Detrend data array
    if level == 'surface':
        x = np.arange(datavar.shape[0])
        
        slopes = np.empty((datavar.shape[1],datavar.shape[2]))
        intercepts = np.empty((datavar.shape[1],datavar.shape[2]))
        std_err = np.empty((datavar.shape[1],datavar.shape[2]))
        for i in range(datavar.shape[1]):
            for j in range(datavar.shape[2]):
                mask = np.isfinite(datavar[:,i,j])
                y = datavar[:,i,j]
                
                if np.sum(mask) == y.shape[0]:
                    xx = x
                    yy = y
                else:
                    xx = x[mask]
                    yy = y[mask]
                
                if np.isfinite(np.nanmean(yy)):
                    slopes[i,j],intercepts[i,j], \
                    r_value,p_value,std_err[i,j] = sts.linregress(xx,yy)
                else:
                    slopes[i,j] = np.nan
                    intercepts[i,j] = np.nan
        print('Completed: Detrended data for each grid point!')

    print('\n>>> Completed: Finished detrendDataR function!')
    return slopes,std_err

###############################################################################
###############################################################################
###############################################################################

def mk_test(x, alpha):
    """
    This function is derived from code originally posted by Sat Kumar Tomer
    (satkumartomer@gmail.com)
    See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm
    The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
    1987) is to statistically assess if there is a monotonic upward or downward
    trend of the variable of interest over time. A monotonic upward (downward)
    trend means that the variable consistently increases (decreases) through
    time, but the trend may or may not be linear. The MK test can be used in
    place of a parametric linear regression analysis, which can be used to test
    if the slope of the estimated linear regression line is different from
    zero. The regression analysis requires that the residuals from the fitted
    regression line be normally distributed; an assumption not required by the
    MK test, that is, the MK test is a non-parametric (distribution-free) test.
    Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
    viewed as an exploratory analysis and is most appropriately used to
    identify stations where changes are significant or of large magnitude and
    to quantify these findings.
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)
    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics
    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05)
    """
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Import modules
    import numpy as np
    from scipy.stats import norm
    
    n = len(x)

    # calculate S
    s = 0
    for k in range(n-1):
        for j in range(k+1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n*(n-1)*(2*n+5) - np.sum(tp*(tp-1)*(2*tp+5)))/18

    if s > 0:
        z = (s - 1)/np.sqrt(var_s)
    elif s < 0:
        z = (s + 1)/np.sqrt(var_s)
    else: # s == 0:
        z = 0

    # calculate the p_value
    p = 2*(1-norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1-alpha/2)

    if (z < 0) and h:
        trend = 'decreasing'
    elif (z > 0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'
        
    ### Significant at 95% confidence level
    if p >= alpha:
        p = np.nan
    elif p < alpha:
        p = 1.

    return trend, h, p, z