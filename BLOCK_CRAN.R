# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# devtools::install_github("dkahle/ggmap")
# install.packages("~/Downloads/spatial.tools_1.6.2.tar.gz", repos=NULL, type="source")
# devtools::install_github("SantanderMetGroup/climate4R.datasets")
# devtools::install_github(c("SantanderMetGroup/loadeR.java",
#"SantanderMetGroup/climate4R.UDG",
#"SantanderMetGroup/loadeR",
#"SantanderMetGroup/transformeR",
#"SantanderMetGroup/visualizeR",
#"SantanderMetGroup/downscaleR"))
# install.packages("~/Downloads/rcompanion_2.3.26.tar.gz", repos=NULL, type="source", dependencies = T)
# devtools::install_github("pacificclimate/ClimDown")
# install.packages("~/Downloads/climdex.pcic_1.1-11.tar.gz", repos=NULL, type="source")

# Download IDFtool from Github
# https://github.com/dazamora/IDFtool/

# Library devtools MUST be previously installed
# devtools::install_github("dazamora/IDFtool")

# Library investr MUST be downloaded from CRAN and installed as source
# directly in RStudio Console as:
# install.packages("~/Downloads/investr_1.4.0.tar.gz", repos=NULL, type="source")
# devtools::install_github("pacificclimate/climdex.pcic")
# install_github("ClimDesign/fixIDF",ref="main")
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(qmap)
require(hydroGOF)
require(tdr)
require(dplyr)
require(ggplot2)
require(hyfo)
require(lubridate)
require(plyr)
require(transformeR)
require(pastecs)
require(climdex.pcic)
require(IDFtool)
require(openxlsx)
require(evd)
require(nsRFA)
# require(raster)
# require(sf)
# require(car)
# require(chron)
# require(compare)
# require(DescTools)
# require(tidyr)
# library(purrr)
# require(forecast)
# require(ggmap)
# require(intamap)
# require(KScorrect)
# require(lattice)
# require(maps)
# require(maptools)
# require(MASS)
# require(matrixStats)
# require(ncdf4)
# require(psych)
# require(rcompanion)
# require(readr)
# require(readxl)
# require(reshape)
# require(reshape2)
# require(rgdal)
# require(rgeos)
# require(RNetCDF)
# require(scales)
# require(mmap)
# require(spatial.tools)
# require(spatialEco)
# require(SpatialTools)
# require(tibble)
# require(tidyr)
# require(viridis)
# require(weathermetrics)
# require(rstudioapi)
# require(eurocordexr)
# require(imputeTS)
# require(MBC)
# require(EnvStats)
# require(investr)
# require(stringr)
# require(precintcon)
# require(trend)

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------
