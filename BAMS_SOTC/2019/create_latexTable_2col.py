"""
Script creates a latex table for manuscript (two-column)

Notes
-----
    Author : Zachary Labe
    Date   : 22 February 2018
"""

### Import modules
import numpy as np
import datetime
from tabulate import tabulate

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Create LaTeX Table (2 col)- %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

def readData(vari):
    """
    Read in data as numpy array for generating latex table
    
    Parameters
    ----------
    vari : string 
        variable to retrieve data file       
        
    Returns
    -------
    datay : 2d array
        [variables,month] ---> should be 6 months (Oct-Mar)

    Usage
    -----
    datay = readData(vari)
    
    """
    print('\n>>> Using readData function!')
    
    ### Select file
    if vari == 'RMSE':
        filename = directorydata2 + 'rmse_DJF.txt'
    elif vari == 'patterncorr':
        filename = directorydata2 + 'patterncorr_DJF.txt'
    elif vari == 'ratio':
        filename = directorydata2 + 'sicsitratio_DJF.txt'
    else:
        ValueError('Variable statistics not calculated!')
    
    ### Read in data    
    data = np.genfromtxt(filename,skip_header=2,delimiter=',')
    print('Completed: Read in data for %s!' % vari)
    
    ### Transpose table so it is [variables,months]. This allows for months
    ### as the column headers
    if any([data.shape[0] == 6, data.shape[0] == 3]):
        data = data.transpose()
    
    print('*Completed: Finished readData function!')
    return data

### Select variable for script
var1 = 'ratio'
var2 = 'patterncorr'

### Call function read in variable for table
data1 = readData(var1)
data2 = readData(var2)

datay = np.array([data1,data2]).transpose()

#################################################################################
#################################################################################
#################################################################################
#### Create LaTeX table

### Enter row titles
rowIDs = ['U10','Z30','U300','Z500','SLP','T2M','RNET']

### Enter column titles

if datay.shape[1] == 6:
    cols = ['October','November','December','January','February','March']
elif datay.shape[1] == 3:
    cols = ['OCT-NOV','DEC-JAN','FEB-MAR']
elif datay.shape[1] == 2:
    cols = ['Response Contribution (\%)','Pattern Correlations (r)']
else:
    print(ValueError('Something is wrong with the data!'))

### Build table
table = tabulate(datay,tablefmt='latex',showindex=rowIDs,headers=cols,
                 numalign='center')

### Save table
file = open(directorydata2 + 'DJF_latextable.txt','w')
file.write(table)
file.close()
