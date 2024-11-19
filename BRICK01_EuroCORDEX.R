# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Construction Engineering School
# MSc.Eng. Maikel Mendez Morales
# https://www.tec.ac.cr
# Email: maikel.mendez@gmail.com; mamendez@itcr.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# Skype: maikel.mendez
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INFO: This script is intended for the extraction of Daily CORDEX-Datasets based on specific Lat/Long position
# for specific weather-stations. Its structure will change as a function of the NetCDF temporal array
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
# List_CAM22_HadGEM2_RegCM47_hist.txt: NetCDF list of historical CAM22 for specific GCM-RCM member
# List_CAM22_HadGEM2_RegCM47_sim.txt: NetCDF list of simulated future CAM22 for specific GCM-RCM member and RCP
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# BRICK01_df_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing dataset vectors of cell-to-point time-series extraction
# date: date as month/day/year
# prec: precipitation [mm/day]
# grid_lon: Weather station longitude [degrees]
# grid_lat: Weather station latitude [degrees]
# sta_id: Weather station ID
# rcm_id: CORDEX member ID
# crtm05_x: CRTM05 X [m]
# crtm05_y: CRTM05 X [y]
#-------------------------------------------------------------------------------------------------------------------

# Workspace is cleared
gc(); rm(list = ls())

# Scientific notation is disabled
options(scipen=999)

# Start time is recorded
start.time <- Sys.time()

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/HadGEM2_RegCM47_AIJS")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(gstat)
require(sp)
require(tidyverse)
require(lubridate)
require(DescTools)
require(ggplot2)
require(matrixStats)
require(pastecs)
require(ncdf4)
require(RNetCDF)
require(rgdal)
require(eurocordexr)
require(openxlsx)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX + Weather Station Metadata + Input Data
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Weather station ID
#-------------------------------------------------------------------------------------------------------------

# Reading Weather station attributes values and parameters
attri <- read.delim("attrifile.txt", header = TRUE, sep = "\t")

# Weather station ID is defined
sta.id <- as.character(attri[1, 2]) # Aeropuerto Internacional Juan Santamaria

# RCM ID is defined
rcm.id <- as.character(attri[2, 2])

# RCP ID is defined
rcp.id <- as.character(attri[3, 2])

# Weather station coordinates are defined
sta.long <- as.numeric(attri[4, 2]) # Aeropuerto Internacional Juan Santamaria
sta.lat <- as.numeric(attri[5, 2])  # Aeropuerto Internacional Juan Santamaria

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: CORDEX Historical
#-------------------------------------------------------------------------------------------------------------

# A path to Historical CORDEX  data is defined
cordex.path <- as.character(attri[6, 2])

# Input data is loaded and a data.frame is created for Historical CORDEX
df.base.cordex <- read.table(as.character(attri[7, 2]), header = FALSE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX:CAM22 Data Historical Extraction
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A check on data directory is performed
path <- cordex.path 
dat_inv <- get_inventory(path)
inv_check <- check_inventory(dat_inv)
inv_check

# A list of data files is requested
path <- cordex.path 
dat <- get_inventory(path)
print(dat)
dat_file <- get_inventory(path, add_files = TRUE)
print(dat_file)

# Pixel counter is defined (years)
rcm.counter <- length(df.base.cordex[[1]])

# A list container is defined
df.list <- NULL
df.list <- list()

# ==============================
# Outermost loop is initialized
# ==============================

# RCP model selector is defined
for (j in 1:rcm.counter) {
  
  # Individual CORDEX files are accessed
  fn1 <- df.base.cordex[[1]][j]
  
  # rotpole_nc_point_to_dt {eurocordexr} is used to extract time series
  df.rot.extract <- rotpole_nc_point_to_dt(filename = fn1,
                                           variable = "pr", # this can vary as a function of CORDEX member
                                           point_lon = sta.long,
                                           point_lat = sta.lat,
                                           verbose = TRUE,
                                           interpolate_to_standard_calendar = TRUE,
                                           add_grid_coord = TRUE)
  # Last data.frame row is isolated
  last_row <- tail(df.rot.extract, n=1)
  
  # December 31 is isolated
  last_row[1,1] <- last_row[1,1] + 1
  
  # A non-precipitation values is assigned
  last_row[1,2] <- 0
  
  # Leap year is isolated
  leap.p <- last_row[1,1]
  
  # Leap year is converted to Date class
  leap.p <- as.Date(leap.p$date)
  
  # Leap year is tested
  lubridate::leap_year(leap.p-1)
  
    # December 31 is added if needed
  if(length(df.rot.extract$date) == 364){
    df.rot.extract <- rbind(df.rot.extract, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract$date) == 365 & lubridate::leap_year(leap.p-1) == TRUE){
    df.rot.extract <- rbind(df.rot.extract, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract$date) == 1825){
    df.rot.extract <- rbind(df.rot.extract, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract$date) == 1826 & lubridate::leap_year(leap.p-1) == TRUE){
    df.rot.extract <- rbind(df.rot.extract, last_row)}
  
  # Units are transformed from [kg m-2 s-1] to [mm/day]
  df.rot.extract$pr <- df.rot.extract$pr*86400
  
  # Data containers are filled
  df.list[[j]] <- df.rot.extract
  
  # ============================
} # End of outermost-loop
  # ============================

# An external list of data.frames is created
df.export.rcm01 <- plyr::ldply(df.list, data.frame) # Changes with RCP !!!

# Relevant values are rounded to 6 decimals
df.export.rcm01$pr <- round(df.export.rcm01$pr, 6)

# Weather station ID is added
df.export.rcm01$sta_id <- sta.id

# RCM ID is added
df.export.rcm01$rcm_id <- rcm.id

# A ggplot object is created
ggraph01 <- ggplot() +
  geom_point(aes(x = date,y = pr,colour = rcm_id), data=df.export.rcm01, size = 0.50) +
  geom_path(aes(x = date,y = pr,colour = sta_id), data=df.export.rcm01, linewidth = 0.50) +
  geom_hline(data = df.export.rcm01, yintercept = 200.0,linewidth = 0.50,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  ggtitle(paste(rcm.id, "GCM-RCM Raw daily historical precipition for Weather Station")) +
  xlab("Day of the year") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(text=element_text(size=14, family="serif"), legend.position="bottom")

# A ggplot object is requested
ggraph01

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: CORDEX RCP Future Simulated
#-------------------------------------------------------------------------------------------------------------

# A path to Future RCP CORDEX data is defined
cordex.path02 <- as.character(attri[8, 2])

# Input data is loaded and a data.frame is created for Future CORDEX
df.base.cordex02 <- read.table(as.character(attri[9, 2]), header = FALSE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX:CAM22 Data Simulated Extraction
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A check on data directory is performed
path <- cordex.path02 
dat_inv <- get_inventory(path)
inv_check <- check_inventory(dat_inv)
inv_check

# A list of data files is requested
path <- cordex.path02 
dat <- get_inventory(path)
print(dat)
dat_file <- get_inventory(path, add_files = TRUE)
print(dat_file)

# Pixel counter is defined
rcm.counter02 <- length(df.base.cordex02[[1]])

# A list container is defined
df.list02 <- NULL
df.list02 <- list()

# ==============================
# Outermost loop is initialized
# ==============================

# RCP model selector is defined
for (j in 1:rcm.counter02) {
  
  # Individual CORDEX files are accessed
  fn2 <- df.base.cordex02[[1]][j]
  
  # rotpole_nc_point_to_dt {eurocordexr} is used to extract time series
  df.rot.extract02 <- rotpole_nc_point_to_dt(filename = fn2,
                                             variable = "pr", # this can vary as a function of CORDEX member
                                             point_lon = sta.long,
                                             point_lat = sta.lat,
                                             verbose = TRUE,
                                             interpolate_to_standard_calendar = TRUE,
                                             add_grid_coord = TRUE)
  
  # Last data.frame row is isolated
  last_row <- tail(df.rot.extract02, n=1)
  
  # December 31 is isolated
  last_row[1,1] <- last_row[1,1] + 1
  
  # A non-precipitation values is assigned
  last_row[1,2] <- 0
  
  # Leap year is isolated
  leap.p <- last_row[1,1]
  
  # Leap year is converted to Date class
  leap.p <- as.Date(leap.p$date)
  
  # Leap year is tested
  lubridate::leap_year(leap.p-1)
  
  # December 31 is added if needed
  if(length(df.rot.extract02$date) == 365 & lubridate::leap_year(leap.p-1) == TRUE){
    df.rot.extract02 <- rbind(df.rot.extract02, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract02$date) == 364){
    df.rot.extract02 <- rbind(df.rot.extract02, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract02$date) == 1825){
    df.rot.extract02 <- rbind(df.rot.extract02, last_row)}
  
  # December 31 is added if needed
  if(length(df.rot.extract02$date) == 1826 & lubridate::leap_year(leap.p-1) == TRUE){
    df.rot.extract02 <- rbind(df.rot.extract02, last_row)}
  
  # Units are transformed from [kg m-2 s-1] to [mm/day] 
  df.rot.extract02$pr <- df.rot.extract02$pr*86400
  
  # Data containers are filled
  df.list02[[j]] <- df.rot.extract02
  
  # ==============================
} # End of outermost-loop
  # ==============================

# An external list of data.frames is created
df.export.rcm02 <- plyr::ldply(df.list02, data.frame) # Changes with RCP !!!

# Relevant values are rounded to 6 decimals
df.export.rcm02$pr <- round(df.export.rcm02$pr, 6)

# Weather station ID is added
df.export.rcm02$sta_id <- sta.id

# RCM ID is added
df.export.rcm02$rcm_id <- rcm.id

# A ggplot object is created
ggraph02 <- ggplot() +
  geom_point(aes(x = date,y = pr,colour = rcm_id), data=df.export.rcm02, size = 0.50) +
  geom_path(aes(x = date,y = pr,colour = sta_id), data=df.export.rcm02, linewidth = 0.50) +
  geom_hline(data = df.export.rcm02, yintercept = 200.0,linewidth = 0.50,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  ggtitle(paste(rcm.id, "Raw daily historical precipition for Weather Station")) +
  xlab("Day of the year") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(text=element_text(size=14, family="serif"), legend.position="bottom")

# A ggplot object is requested
ggraph02

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Stationarity test of time series
#-------------------------------------------------------------------------------------------------------------

# A time-series object is created
t.sta.pop <- ts(df.export.rcm01$pr, frequency = 365)
t.sta.pop2 <- ts(df.export.rcm02$pr, frequency = 365)

# Simple time-series object are plotted for verification purposes
plot(t.sta.pop, main = "Simple adf.test plot for historical time series")
plot(t.sta.pop2, main = "Simple adf.test plot for simulated time series")

# Augmented Dickey-Fuller Test is requested
# If p-value > 0.05 = nonstationary; If p-value < 0.05 = stationary
# Since the p-value < 0.05, we conclude that there is enough evidence to reject 
# the Null hypothesis, meaning that the time series is indeed stationary
tseries::adf.test(t.sta.pop)
tseries::adf.test(t.sta.pop2)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: epsg coordinates transformation
#-------------------------------------------------------------------------------------------------------------

# Weather station coordinates are defined
df.coord <- data.frame(x=c(sta.long),y=c(sta.lat))

# Weather station coordinates are assigned
coordinates(df.coord) <- ~ x + y

# WGS epsg is defined
proj4string(df.coord) <- CRS("+init=epsg:4326")

# WGS coordinates are transformed to CRTM05
temp.coord <- spTransform(df.coord, CRS("+init=epsg:5367"))

# CRTM05 coordinates are transformed to vector
temp.coord <- as.vector(temp.coord@coords)

# CRTM05 coordinates are discretized 
crtm05.x <- round(temp.coord[1], 3)
crtm05.y <- round(temp.coord[2], 3)

# -------------------------------------------------------
# CRTM05 Coordinates are added to relevant objects
# -------------------------------------------------------

# crtm05.x is added
df.export.rcm01$crtm05_x <- crtm05.x
df.export.rcm02$crtm05_x <- crtm05.x

# crtm05.y is added
df.export.rcm01$crtm05_y <- crtm05.y
df.export.rcm02$crtm05_y <- crtm05.y

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# IDF historical intensities are requested
View(df.export.rcm01) # Hist. 13149 obs 1970-2005
View(df.export.rcm02) # Fut.  34394 obs 2005-2099

# data.frame variables are renamed
names(df.export.rcm01)[names(df.export.rcm01) == "pr"] <- c("prec")
names(df.export.rcm02)[names(df.export.rcm02) == "pr"] <- c("prec")

# An output file name is concatenated
t.output01 <- paste("BRICK01_df_",rcm.id,"_",rcp.id,"_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK01_df_",rcm.id,"_",rcp.id,"_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(df.export.rcm01,
     df.export.rcm02,
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("df.export.rcm01" = df.export.rcm01, # Hist. 13149 obs 1970-2005
                   "df.export.rcm02" = df.export.rcm02) # Fut.  34394 obs 2005-2099
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
