# ClimatePython [![DOI](https://zenodo.org/badge/135844507.svg)](https://zenodo.org/badge/latestdoi/135844507)
This repository contains a variety of Python tutorials and lectures I've used in climate/atmospheric sciences. Check back soon!

###### Under construction... ```[Python 3.6]```

## Contact
Zachary Labe - [Research Website](http://sites.uci.edu/zlabe/) - [@ZLabe](https://twitter.com/ZLabe)

## Publications
+ Timmermans, M.-L., **Z.M. Labe**, and C. Ladd (2020). [The Arctic] Sea surface temperature [in “State of the Climate in 2019”], *Bull. Amer. Meteor. Soc.*, DOI:10.1175/BAMS-D-20-0086.1 [[HTML]](https://journals.ametsoc.org/bams/article/101/8/S239/353884/The-Arctic)[[BibTeX]](https://sites.uci.edu/zlabe/files/2020/08/TheArctic_BAMS_SOTC2019_BibTeX.pdf)

## Description
+ ```ScienceVisuals/```: Lectures on improving science figures
    + ```LineGraph/```: Script on improving line graphs 
    + ```MapGraph/``` : Script on improving map plots
    + ```Presentations/```: Powerpoints on discussing science visualization [[SlideShare]](https://www.slideshare.net/ZacharyLabe)
+ ```Data/```: Additional data files not provided by Python URL functions
+ ```Figures/```: Arbitrary figures as examples from listed scripts
+ ```Scripts/```: Main [Python](https://www.python.org/) scripts/functions used in data analysis and plotting. 
+ ```requirements.txt```: List of environments and modules associated with the most recent version of this project. A Python [Anaconda3 Distribution](https://docs.continuum.io/anaconda/) was used for the analysis. Tools including [NCL 6.4.0](https://www.ncl.ucar.edu/), [CDO](https://code.mpimet.mpg.de/projects/cdo), and [NCO](http://nco.sourceforge.net/) were also used for initial data manipulation. [ImageMagick](https://www.imagemagick.org/script/index.php) is used for most of the animations (GIF). All code has been tested with Python ```3.6```.

## Data
###### Global Surface Temperature Data Sets (station-based)
+ Berkeley Earth Surface Temperature Project (BEST) : [[DATA]](http://berkeleyearth.org/data/)
    + Rohde, R., et al. "A new estimate of the average Earth surface land temperature spanning 1753 to 2011. Geoinfor Geostat Overview 1: 1." of 7 (2013): 2. [[Publication]](https://www.scitechnol.com/new-estimate-of-the-average-earth-surface-land-temperature-spanning-to-1eCc.php?article_id=450)
+ Cowtan and Way surface temperature record (C&W) : [[DATA]](http://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html)
    + Cowtan, Kevin, and Robert G. Way. "Coverage bias in the HadCRUT4 temperature series and its impact on recent temperature trends." Quarterly Journal of the Royal Meteorological Society 140.683 (2014): 1935-1944. [[Publication]](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.2297)
+ Hadley Centre/CRU gridded surface temperature data set version 4 (HadCRUT4) : [[DATA]](https://crudata.uea.ac.uk/cru/data/temperature/)
    + Morice, Colin P., et al. "Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: The HadCRUT4 data set." Journal of Geophysical Research: Atmospheres 117.D8 (2012). [[Publication]](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011JD017187)
+ NASA GISS Surface Temperature Analysis (GISTEMP) : [[DATA]](https://data.giss.nasa.gov/gistemp/)
    + Hansen, James, et al. "Global surface temperature change." Reviews of Geophysics 48.4 (2010). [[Publication]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010RG000345)
+ NOAA Merged Land Ocean Global Surface Temperature (NOAAGlobalTemp) : [[DATA]](https://www.ncdc.noaa.gov/data-access/marineocean-data/noaa-global-surface-temperature-noaaglobaltemp)
    + Vose, Russell S., et al. "NOAA's merged land–ocean surface temperature analysis." Bulletin of the American Meteorological Society 93.11 (2012): 1677-1685. [[Publication]](https://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00241.1)
###### Reanalysis Data 
+ ERA5 : [[DATA]](http://apps.ecmwf.int/data-catalogues/era5/?class=ea)
+ NCEP/NCAR Reanalysis 1 (R1): [[DATA]](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html)
    + Kalnay, E., and co-authors, 1996: The NCEP/NCAR 40-year reanalysis project. Bulletin of the American meteorological Society, 77(3), 437-471 [[Publication]](http://journals.ametsoc.org/doi/abs/10.1175/1520-0477(1996)077%3C0437:TNYRP%3E2.0.CO;2)
###### Sea Surface Temperature Data 
+ NOAA Optimum Interpolation Sea Surface Temperature V2 (OISSTv2): [[DATA]](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)
    + Reynolds, R. W., Rayner, N. A., Smith, T. M., Stokes, D. C., & Wang, W. (2002). An improved in situ and satellite SST analysis for climate. Journal of climate, 15(13), 1609-1625. [[Publication]](https://journals.ametsoc.org/doi/full/10.1175/1520-0442%282002%29015%3C1609%3AAIISAS%3E2.0.CO%3B2)
