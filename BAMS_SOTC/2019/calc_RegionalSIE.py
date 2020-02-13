def read_RegionalSIE(regionsea,mon):
    """
    Retrieve data from NSIDC on regional sea ice extent (daily)
    """
    
    ### Import modules
    import numpy as np
    import pandas as pd

    ### Select parameters
    years = np.arange(1982,2019+1,1)
    
    ### Load url
    url = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/seaice_analysis/' \
            'Sea_Ice_Index_Regional_Monthly_Data_G02135_v3.0.xlsx'
    
    ### Read file
    time = ((mon-1)*2)+1
    df_reg = pd.read_excel(url,sheet_name='%s-Extent-km^2' % (regionsea),
                               header=3,
                               usecols=range(time,time+1,1))
    reg = df_reg.values
    
    ### Convert
    regsieq = reg/1e6
    
    ### Plot through 2019
    regsie = regsieq.squeeze()[3:-1]
    
    print('\nCompleted: Read sea ice data for %s!' % (regionsea))           
    return regsie

#### Call functions
regsie = read_RegionalSIE('Barents',8)
print(regsie)