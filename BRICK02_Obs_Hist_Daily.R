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
# INFO: This script is intended for the extraction of Daily Observed Historical datasets
# for specific weather-stations
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:

# Lluvia_Horaria_84_169.txt: hourly historical records
# date: date as day/month/year
# time: hour of day from 100 to 2400
# prec: precipitation [mm/hour]

# Lluvia_diaria_84_169.txt: daily historical records
# date: date as day/month/year
# prec: precipitation [mm/day]
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# BRICK02_df_obs_AIJS: R-object containing dataset vectors of daily historical records with the following 
# objects: df_obs_1970_2005, df_obs_1970_1990, df_obs_1991_2005, df_obs_1970_2022,
# date: date as month/day/year
# prec: precipitation [mm/day]
# grid_lon: Weather station longitude [degrees]
# grid_lat: Weather station latitude [degrees]
# sta_id: Weather station ID
# imn_id: IMN weather-sations included
# crtm05_x: CRTM05 X [m]
# crtm05_y: CRTM05 X [y]

# objects: df.base.daily.parsing
# date: date as month/day/year
# prec: precipitation [mm/day]
# max_hr: max.precipitation [mm/hour]
# min_hr: min.precipitation [mm/hour]
# mean_hr: mean.precipitation [mm/hour]
# grid_lon: Weather station longitude [degrees]
# grid_lat: Weather station latitude [degrees]
# sta_id: Weather station ID
# imn_id: IMN weather-sations included
# crtm05_x: CRTM05 X [m]
# crtm05_y: CRTM05 X [y]
#-------------------------------------------------------------------------------------------------------------------

# Workspace is cleared
# gc(); rm(list = ls())

# Scientific notation is disabled
options(scipen=999)

# Start time is recorded
start.time <- Sys.time()

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/HadGEM2_RegCM47_AIJS")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(dplyr)
require(gstat)
require(sp)
require(rgdal)
require(tidyverse)
require(lubridate)
require(DescTools)
require(ggplot2)
require(matrixStats)
require(anytime)
require(imputeTS)
require(tsbox)
require(openxlsx)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Weather Station Metadata + Input Data
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
# SUB-BLOCK: Weather Matadata
#-------------------------------------------------------------------------------------------------------------

# Hourly observed data are loaded
df.base.hour <- read.table("Lluvia_Horaria_84_169.txt",header=T,sep="\t",quote="")

# Daily observed data are loaded
df.base.daily01 <- read.table("Lluvia_diaria_84_169.txt",header=T,sep="\t",quote="")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Hourly data is aggregated into 24-hr (daily) series and NAs are imputated
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Time is transformed to 0-23 hr range
df.base.hour$time <- ((df.base.hour$time)/100)-1

# Variables data + time are converted to as.POSIXct objects
df.base.hour$date <- as.POSIXct(paste(df.base.hour$date, df.base.hour$time), format="%d/%m/%Y %H")

# Incomplete days are removed
df.base.hour <- df.base.hour[16:202517, ]

# data.frame head is checked, it MUST return 00:00:00
head(df.base.hour, n = 1)

# data.frame tail is checked, it MUST return 23:00:00
tail(df.base.hour, n = 1)

# Dimensions of date fields are checked for consistency
# They MUST exhibit the same length
min(str_length(df.base.hour$date))
max(str_length(df.base.hour$date))

# NAs presence is evaluated
summary(is.na(df.base.hour))

# Hourly data is aggregated into 24-hr (daily) series
df.base.daily.parsing <- df.base.hour %>%
  mutate(date_time = as.POSIXct(date,format="%Y-%m-%d %H:%M:%S")) %>%
  group_by(date(date)) %>%
  summarise(sum = sum(prec,na.rm=TRUE),
            max = max(prec,na.rm=TRUE),
            min = min(prec,na.rm=TRUE),
            mean = round(mean(prec,na.rm=TRUE), 5))

# NAs presence is evaluated
summary(is.na(df.base.daily.parsing))

# data.frame variables are renamed
names(df.base.daily.parsing) <- c("date", "prec", "max_hr", "min_hr", "mean_hr")

# A prec. ONLY subset is created
df.base.daily02 <- df.base.daily.parsing[, c(1,2)]

# A date summary is requested for field date
summary(df.base.daily02$date)

# Dimensions of date field are checked
# They MUST exhibit the same length
min(str_length(df.base.daily02$date))
max(str_length(df.base.daily02$date))

# A complete synthetic timeseries is created to check for NAs
df.dates.base01 <- as.data.frame(seq(as.Date(min(df.base.daily02$date)),
                                     as.Date(max(df.base.daily02$date)),
                                     by = "1 day"))

# data.frame variables are renamed
names(df.dates.base01) <- c("date")

# data.frame head is requested
head(df.dates.base01$date)

# data.frame tail is requested
tail(df.dates.base01$date)

# A dummy precipitation variables is created
df.dates.base01$prec <- 0

# mutate-joins {dplyr} is applied to pinpoint NAs
df.blended.01 <- full_join(df.dates.base01, df.base.daily02, by = c("date"), keep = FALSE)

# A blended prec field is created
df.blended.01$prec <- (df.blended.01$prec.x + df.blended.01$prec.y)

# Blended precipitation + date fields are isolated
df.isolate01 <- df.blended.01[, c(1,4)]

# A time-series object is created
ts.isolate01 <- ts_ts(ts_long(df.isolate01))

# Object class is requested
class(ts.isolate01)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Isolation of missing values
#-------------------------------------------------------------------------------------------------------------

# BEWARE !!! A summary on missing values (NAs) is requested
print(statsNA(ts.isolate01))

# ggplot_na_distribution is used to Visualize the Distribution of Missing Values
ggplot_na_distribution(ts.isolate01)

# ggplot_na_gapsize is used to Visualize Occurrences of NA gap sizes
ggplot_na_gapsize(ts.isolate01)

# A missing values-only data.frame is created
df.isolate02 <- df.isolate01[is.na(df.isolate01$prec), ]

# An export/import XLS object is created
write.xlsx(df.isolate02, file = "missing.xlsx", rowNames = TRUE)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Data auto-imputation by {imputeTS}
#-------------------------------------------------------------------------------------------------------------

# na_kalman {imputeTS} is applied to imputate NAs by Kalman Smoothing and State Space Models

# Option 01
ts.kalman01 <- na_seasplit(ts.isolate01, algorithm = "interpolation")

# Option 02
#ts.kalman01 <- na_kalman(ts.isolate01, model = "auto.arima", smooth = TRUE)#, type = "trend")

# A ggplot2 object is requested
ggplot_na_imputations(ts.isolate01, ts.kalman01)

# time-series object is transformed to vector
ts.kalman01 <- as.vector(ts.kalman01)

# Imputated vector is added to data.frame
df.blended.01$prec_def <- ts.kalman01

# Variables are rounded to 6 decimal
df.blended.01$prec_def <- round(df.blended.01$prec_def, 6)

# Blended precipitation + date fields are isolated
df.blended.export.01 <- df.blended.01[, c(1,5)]

# data.frame variables are renamed
names(df.blended.export.01) <- c("date", "prec")

# A ggplot object is created
ggraph01 <- ggplot() +
  geom_path(aes(x = date,y = prec), colour = "blue", data=df.blended.export.01, linewidth = 0.50) +
  geom_point(aes(x = date,y = prec), colour = "black", data=df.blended.export.01, size = 1.50) +
  geom_hline(yintercept = 100.0,linewidth = 0.50,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, "Daily historical precipition for Weather Station")) +
  xlab("Day of the year") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(text=element_text(size=14, family="serif"), legend.position="bottom")

# A ggplot object is requested
ggraph01

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Daily dataset is analized and NAs are imputated
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Dimensions of date field are checked
# They MUST exhibit the same length
min(str_length(df.base.daily01$date))
max(str_length(df.base.daily01$date))

# as_date {lubridate} is used to assign dates to variables
df.base.daily01$date <- as.Date(df.base.daily01$date, format = "%d/%m/%Y")

# A complete synthetic time series is created to check for NAs
df.dates.base02 <- as.data.frame(seq(as.Date(min(df.base.daily01$date)),
                                     as.Date(max(df.base.daily01$date)),
                                     by = "1 day"))

# data.frame variables are renamed
names(df.dates.base02) <- c("date")

# data.frame head is requested
head(df.dates.base02$date)

# data.frame tail is requested
tail(df.dates.base02$date)

# A dummy precipitation variables is created
df.dates.base02$prec <- 0

# mutate-joins {dplyr} is applied to pinpoint NAs
df.blended.02 <- full_join(df.dates.base02, df.base.daily01, by = c("date"), keep = FALSE)

# A blended prec field is created
df.blended.02$prec <- (df.blended.02$prec.x + df.blended.02$prec.y)

# Blended precipitation + date fields are isolated
df.isolate02 <- df.blended.02[, c(1,4)]

# A time-series object is created
ts.isolate02 <- ts_ts(ts_long(df.isolate02))

# Object class is requested
class(ts.isolate02)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Data imputation by {imputeTS}
#-------------------------------------------------------------------------------------------------------------

# A summary is requested
print(statsNA(ts.isolate02))

# ggplot_na_distribution is used to Visualize the Distribution of Missing Values
ggplot_na_distribution(ts.isolate02)

# ggplot_na_gapsize is used to Visualize Occurrences of NA gap sizes
ggplot_na_gapsize(ts.isolate02)

# na_kalman {imputeTS} is applied to imputate NAs by Kalman Smoothing and State Space Models

# Option 01
ts.kalman02 <- na_seasplit(ts.isolate02, algorithm = "interpolation")

# Option 02
#ts.kalman02 <- na_kalman(ts.isolate02, model = "auto.arima", smooth = TRUE)#, type = "trend")

# A ggplot2 object is requested
ggplot_na_imputations(ts.isolate02, ts.kalman02)

# time-series object is transformed to vector
ts.kalman02 <- as.vector(ts.kalman02)

# Imputated vector is added to data.frame
df.blended.02$prec_def <- ts.kalman02

# Variables are rounded to 1 decimal
df.blended.02$prec_def <- round(df.blended.02$prec_def, 3)

# Blended precipitation + date fields are isolated
df.blended.export.02 <- df.blended.02[, c(1,5)]

# data.frame variables are renamed
names(df.blended.export.02) <- c("date", "prec")

# A ggplot object is created
ggraph02 <- ggplot() +
  geom_path(aes(x = date,y = prec), colour = "black", data=df.blended.export.02, linewidth = 0.50) +
  geom_point(aes(x = date,y = prec), colour = "blue", data=df.blended.export.02, size = 0.50) +
  geom_hline(yintercept = 100.0,linewidth = 0.50,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  ggtitle(paste(sta.id, "Daily historical precipition for Weather Station")) +
  xlab("Day of the year") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(text=element_text(size=14, family="serif"), legend.position="bottom")

# A ggplot object is requested
ggraph02

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Daily is blended into 24-hr (daily) series ONLY
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Both data.frames are joint to create a complete data.frame 
df.base.daily.union <- full_join(df.blended.export.02,
                                 df.blended.export.01,
                                 by = c("date"), keep = FALSE)

# A dummy variable is created
df.base.daily.union$prec <- -999

# NAs are replaced by by ifelse counterpart
df.base.daily.union <- df.base.daily.union %>% mutate(prec = ifelse(!is.na(prec.x),prec.x,prec.y))

# Blended precipitation + date fields are isolated
df.base.daily.union <- df.base.daily.union[, c(1,4)]

# A ggplot object is created
ggraph03 <- ggplot() +
  geom_path(aes(x = date,y = prec),colour = "black", data=df.base.daily.union, linewidth = 0.50) +
  geom_point(aes(x = date,y = prec),colour = "blue", data=df.base.daily.union, size = 0.50) +
  geom_hline(data = df.base.daily.union, yintercept = 100.0,linewidth = 0.50,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  ggtitle(paste(sta.id, "Daily historical precipition for Weather Station")) +
  xlab("Day of the year") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(text=element_text(size=14, family="serif"), legend.position="bottom")

# A ggplot object is requested
ggraph03

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
df.base.daily.union$grid_lon <- sta.long
df.base.daily.union$grid_lat <- sta.lat
df.base.daily.union$sta_id <- sta.id
df.base.daily.union$imn_id <- c("84-169/84-21")
df.base.daily.union$crtm05_x <- crtm05.x
df.base.daily.union$crtm05_y <- crtm05.y
df.base.daily.parsing$grid_lon <- sta.long
df.base.daily.parsing$grid_lat <- sta.lat
df.base.daily.parsing$sta_id <- sta.id
df.base.daily.parsing$imn_id <- c("84-169/84-21")
df.base.daily.parsing$crtm05_x <- crtm05.x
df.base.daily.parsing$crtm05_y <- crtm05.y

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX/Observational datasets time-slicing
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# cordex.output.CAM-22.ICTP.MOHC-HadGEM2-ES.historical.r1i1p1.RegCM4-7.v0.day.pr
# datetime_start = 1970-01-01T12:00:00Z
# datetime_stop = 2005-12-30T12:00:00Z 

# 36-year observation dataset is selected
df_obs_1970_2005 <- df.base.daily.union[df.base.daily.union$date >= "1970-01-01" & df.base.daily.union$date <= "2005-12-31", ]

# 21-year calibration observation dataset is selected
df_obs_1970_1990 <- df.base.daily.union[df.base.daily.union$date >= "1970-01-01" & df.base.daily.union$date <= "1990-12-31", ]

# 15-year validation observation dataset is selected
df_obs_1991_2005 <- df.base.daily.union[df.base.daily.union$date >= "1991-01-01" & df.base.daily.union$date <= "2005-12-31", ]

# 52-year IDF observation dataset is selected
df_obs_1970_2022 <- df.base.daily.union[df.base.daily.union$date >= "1970-01-01" & df.base.daily.union$date <= "2022-12-31", ]

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Stationarity test of time series
#-------------------------------------------------------------------------------------------------------------

# A data.frame is cloned
df_time_analysis <- df_obs_1970_2022[, c(1,2)]

# Year variable is added to data.frame
df_time_analysis$year <- floor_date(df_time_analysis$date, "year")

# Precipitation sum is aggregated by year [mm/year]
df_time_analysis_annual <- df_time_analysis %>%
  group_by(year) %>%
  summarize(mean = sum(prec)) # BEWARE!!! max could be used to analyze daily AMPs

# A time-series object is created
ts_obs_annual <- ts(df_time_analysis_annual$mean, frequency=365)

# A time-series object is plotted
plot(ts_obs_annual, main = "Simple adf.test plot for historical time series")

# Augmented Dickey-Fuller Test is requested
# If p-value > 0.05 = nonstationary; If p-value < 0.05 = stationary
# Since the p-value < 0.05, we conclude that there is enough evidence to reject 
# the Null hypothesis, meaning that the time series is indeed stationary
tseries::adf.test(ts_obs_annual)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Data Export
#-------------------------------------------------------------------------------------------------------------

# An output file name is concatenated
t.output01 <- paste("BRICK02_df_obs","_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK02_df_obs","_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(df_obs_1970_2005,
     df_obs_1970_1990,
     df_obs_1991_2005,
     df_obs_1970_2022,
     df.base.daily.parsing,
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("df_obs_1970_2005" = df_obs_1970_2005,
                   "df_obs_1970_1990" = df_obs_1970_1990,
                   "df_obs_1991_2005" = df_obs_1991_2005,
                   "df_obs_1970_2022" = df_obs_1970_2022,
                   "df_base_daily_parsing" = df.base.daily.parsing)
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
