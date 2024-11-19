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
# INFO: This script is intended for the generation of historical IDFs curves on specific Lat/Long position
# of weather-stations for durations [5, 10, 15, 30, 60, 120, 180, 360, 720, 1440] minutes
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
#
# EST_AIJS_SUB01.txt: Daily AMPs for 5,10,15,30,60,120,360 and 720 min for 1970-1998
# EST_AIJS_SUB02.txt: Daily AMPs for 5,10,15,30,60,120,360 and 720 min for 1980-1998
# EST_AIJS_SUB03.txt: Daily AMPs for 5,10,15 and 30 min for 1997-2023
# Lluvia_Horaria_84_169.txt: Daily hourly precipitation for 1998-2023

# BRICK02_df_obs_AIJS.RData: R-object containing dataset vectors of daily historical records with the following 
# objects: df_obs_1970_2005, df_obs_1970_1990, df_obs_1991_2005, df_obs_1970_2022,
# date: date as month/day/year
# prec: precipitation [mm/day]
# grid_lon: Weather station longitude [degrees]
# grid_lat: Weather station latitude [degrees]
# sta_id: Weather station ID
# imn_id: IMN weather-sations included
# crtm05_x: CRTM05 X [m]
# crtm05_y: CRTM05 X [y]

#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
#
# BRICK06_matrix_bc_AMP_par_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing:
# df.amp.compiled.int             #1: AMPs expressed as int.[mm/day] for 1970-2022
# df.amp.compiled.vol             #2: AMPs expressed as vol.[mm] for 1970-2022
# df.stat.desc.int                #3: Stat. par. as int.[mm/day] for 1970-2022
# df.merge.int.sel                #4: AMPs KS/AD/CvM/AIC/BIC for 1970-2022
# df.para.vol.gev                 #5: GEV AMPs geo.par. Location/Scale/Shape for 1970-2022
# IDF.Auto.int.gev.coefficients   #6: GEV A/B/C reg.mol.parameters for 1970-2022
# IDF.Auto.int.gev.intensities    #7: PDF mod.adjusted int.[mm/day] for 1970-2022
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
require(devtools)
require(investr)
require(IDFtool)
require(DescTools)
require(dplyr)
require(ggplot2)
require(lubridate)
require(matrixStats)
require(pastecs)
require(tidyr)
require(fixIDF)
require(stringr)
require(nsRFA)
require(openxlsx)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Preprocessed binary datasets are loaded
#-------------------------------------------------------------------------------------------------------------

# RData objects are loaded from BRICK02
load("BRICK02_df_obs_AIJS.RData")

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

# Length of validation period is defined [years]
val.length <- as.numeric(attri[10, 2])

# Precipitation threshold is defined
# Marra et.al 2023; Mengel et.al 2021; The Environment and Climate Change Canada (ECCC)
umbral.rain <- as.numeric(attri[11, 2])

# Generalized Pareto Distribution (GPD) theta threshold is defined
p.theta <- as.numeric(attri[12, 2])

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Probability Distribution Fitting (PDF)
#-------------------------------------------------------------------------------------------------------------

# A duration vector is created [min]
v.duration <- c(5, 10, 15, 30, 60, 120, 360, 720, 1440)

# A periods vector is created [years]
v.periods <- c(2, 3, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)

# PDF fitting method is selected
pdf.fit <- "lmoments"

# -------------------------------------------------------------------------------------------------
# L-moments (Lmoments) 
# Argument-label: lmoments

# Probability-Weighted Moments (PWD)
# Argument-label: pwd

# Maximum Likelihood (MLE) 
# Argument-label: mle

# Moments (MME)
# Argument-label: mme
# -------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Observed AMPs
#-------------------------------------------------------------------------------------------------------------

# Sub-hourly observed data are loaded
df.base.sub01 <- read.table("EST_AIJS_SUB01.txt",header=T,sep="\t",quote="")         # 1970-1998

# Sub-hourly observed data are loaded
df.base.sub02 <- read.table("EST_AIJS_SUB02.txt",header=T,sep="\t",quote="")         # 1980-1998

# Sub-hourly observed data are loaded
df.base.sub03 <- read.table("EST_AIJS_SUB03.txt",header=T,sep="\t",quote="")         # 1997-2023

# Hourly observed data are loaded
df.base.hour04 <- read.table("Lluvia_Horaria_84_169.txt",header=T,sep="\t",quote="") # 1998-2023

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Observed dataset cleaning
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1970-1998 sub-hourly dataset
#-------------------------------------------------------------------------------------------------------------------

# date character is converted to class-date and added as a new column
df.base.sub01$date <- as.Date(df.base.sub01$date, format = "%d/%m/%Y")

# Dimensions of date field are checked
min(str_length(df.base.sub01$date))
max(str_length(df.base.sub01$date))

# Only complete days and complete years are included within hourly-data
# This MUST be done manually by inspection !!!
df.base.sub01 <- df.base.sub01[ , ]

# A head/tail (length = 5) is requested for verification
head(df.base.sub01) # complete year !!!
tail(df.base.sub01) # complete year !!!

# lubridate Library functions are applied 
df.base.sub01$YEAR <- lubridate::year(df.base.sub01$date) # Years component of a date-time

# If present, negative values are replaced with 0.0's
df.base.sub01[3][df.base.sub01[3] < 0] <- 0.0

# If present, NAs are replace with 0.0's
df.base.sub01$prec[is.na(df.base.sub01$prec)] <- 0.0

# Sub-hourly observed data are aggregated by year
df.aggre.sub01 <- df.base.sub01 %>%
  group_by(YEAR, instance) %>%
  summarise(prec= max(prec))

# variables are transformed from long to wide format
df.wide.sub01 <- reshape2::dcast(data = df.aggre.sub01,
                       formula = YEAR ~ instance,
                       fun.aggregate = sum,
                       value.var = "prec")

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1980-1998 sub-hourly dataset
#-------------------------------------------------------------------------------------------------------------------

# date character is converted to class-date and added as a new column
df.base.sub02$date <- as.Date(df.base.sub02$date, format = "%d/%m/%Y")

# Dimensions of date field are checked
min(str_length(df.base.sub02$date))
max(str_length(df.base.sub02$date))

# Only complete days and complete years are included within hourly-data
# This MUST be done manually by inspection !!!
df.base.sub02 <- df.base.sub02[ , ]

# A head/tail (length = 5) is requested for verification
head(df.base.sub02) # complete year !!!
tail(df.base.sub02) # complete year !!!

# lubridate Library functions are applied 
df.base.sub02$YEAR <- lubridate::year(df.base.sub02$date) # Years component of a date-time

# If present, negative values are replaced with 0.0's
df.base.sub02[2][df.base.sub02[2] < 0] <- 0.0
df.base.sub02[3][df.base.sub02[3] < 0] <- 0.0
df.base.sub02[4][df.base.sub02[4] < 0] <- 0.0
df.base.sub02[5][df.base.sub02[5] < 0] <- 0.0
df.base.sub02[6][df.base.sub02[6] < 0] <- 0.0
df.base.sub02[7][df.base.sub02[7] < 0] <- 0.0
df.base.sub02[8][df.base.sub02[8] < 0] <- 0.0
df.base.sub02[9][df.base.sub02[9] < 0] <- 0.0

# If present, NAs are replace with 0.0's
df.base.sub02$X5[is.na(df.base.sub02$X5)] <- 0.0
df.base.sub02$X10[is.na(df.base.sub02$X10)] <- 0.0
df.base.sub02$X15[is.na(df.base.sub02$X15)] <- 0.0
df.base.sub02$X30[is.na(df.base.sub02$X30)] <- 0.0
df.base.sub02$X60[is.na(df.base.sub02$X60)] <- 0.0
df.base.sub02$X120[is.na(df.base.sub02$X120)] <- 0.0
df.base.sub02$X360[is.na(df.base.sub02$X360)] <- 0.0
df.base.sub02$X720[is.na(df.base.sub02$X720)] <- 0.0

# Sub-hourly observed data are aggregated by year
df.aggre.02 <- df.base.sub02 %>%
  group_by(YEAR) %>%
  summarise(X5= max(X5),
            X10= max(X10),
            X15= max(X15),
            X30= max(X30),
            X60= max(X60),
            X120= max(X120),
            X360= max(X360),
            X720= max(X720))

# Variables are transformed from long to wide format
df.wide.sub02 <- df.aggre.02

# Precipitation depths [mm] are transformed to intensities [mm/hour]
df.wide.sub02[2] <- df.wide.sub02[2]/(5/60)
df.wide.sub02[3] <- df.wide.sub02[3]/(10/60)
df.wide.sub02[4] <- df.wide.sub02[4]/(15/60)
df.wide.sub02[5] <- df.wide.sub02[5]/(30/60)
df.wide.sub02[6] <- df.wide.sub02[6]/(60/60)
df.wide.sub02[7] <- df.wide.sub02[7]/(120/60)
df.wide.sub02[8] <- df.wide.sub02[8]/(360/60)
df.wide.sub02[9] <- df.wide.sub02[9]/(720/60)

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1970-1998 sub-hourly dataset is completed
#-------------------------------------------------------------------------------------------------------------------

# Names are assigned to matrix columns
colnames(df.wide.sub01) <- c("year", "min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720")

# Names are assigned to matrix columns
colnames(df.wide.sub02) <- c("year", "min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720")

# A data.frame subset is requested given differences between dates, volumes an intensities
df.wide.sub01 <- df.wide.sub01[1:10, ]

# data.frames are rbinded
df.wide.temp <- rbind(df.wide.sub01, df.wide.sub02)

# 0 values are replaced with the mean per column in data.frame
df.wide.temp <-  df.wide.temp %>% mutate(across(.cols = min_5:min_720,
                                                .fns = ~ ifelse(.x == 0, mean(.x), .x)))

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1997-2023 sub-hourly dataset
#-------------------------------------------------------------------------------------------------------------------

# date character is converted to class-date and added as a new column
df.base.sub03$date <- as.Date(df.base.sub03$date, format = "%d/%m/%Y")

# Dimensions of date field are checked
min(str_length(df.base.sub03$date))
max(str_length(df.base.sub03$date))

# Only complete days and complete years are included within hourly-data
# This MUST be done manually by inspection !!!
df.base.sub03 <- df.base.sub03[66:8418 , ]

# A head/tail (length = 5) is requested for verification
head(df.base.sub03) # complete year !!!
tail(df.base.sub03) # complete year !!!

# lubridate Library functions are applied 
df.base.sub03$YEAR <- lubridate::year(df.base.sub03$date) # Years component of a date-time

# If present, negative values are replaced with 0.0's
df.base.sub03[2][df.base.sub03[2] < 0] <- 0.0
df.base.sub03[3][df.base.sub03[3] < 0] <- 0.0
df.base.sub03[4][df.base.sub03[4] < 0] <- 0.0
df.base.sub03[5][df.base.sub03[5] < 0] <- 0.0

# If present, NAs are replace with with 0.0's
df.base.sub03$X5[is.na(df.base.sub03$X5)] <- 0.0
df.base.sub03$X10[is.na(df.base.sub03$X10)] <- 0.0
df.base.sub03$X15[is.na(df.base.sub03$X15)] <- 0.0
df.base.sub03$X30[is.na(df.base.sub03$X30)] <- 0.0

# Sub-hourly observed data are aggregated by year
df.aggre.03 <- df.base.sub03 %>%
  group_by(YEAR) %>%
  summarise(X5= max(X5),
            X10= max(X10),
            X15= max(X15),
            X30= max(X30))

# Variables are transformed from long to wide format
df.wide.sub03 <- df.aggre.03

# Names are assigned to matrix columns
colnames(df.wide.sub03) <- c("year", "min_5",	"min_10",	"min_15",	"min_30")

# Precipitation depths [mm] are transformed to intensities [mm/hour]
df.wide.sub03[2] <- df.wide.sub03[2]/(5/60)
df.wide.sub03[3] <- df.wide.sub03[3]/(10/60)
df.wide.sub03[4] <- df.wide.sub03[4]/(15/60)
df.wide.sub03[5] <- df.wide.sub03[5]/(30/60)

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1998-2023 hourly dataset
#-------------------------------------------------------------------------------------------------------------------

# Time is transformed to 0-23 hr range
df.base.hour04$time <- ((df.base.hour04$time)/100)-1

# Variables data + time are converted to as.POSIXct objects
df.base.hour04$date <- as.POSIXct(paste(df.base.hour04$date,
                                        df.base.hour04$time), format="%d/%m/%Y %H")

# Dimensions of date field are checked
min(str_length(df.base.hour04$date))
max(str_length(df.base.hour04$date))

# Blended precipitation + date fields are isolated
df.base.hour04 <- df.base.hour04[, c(1,3)]

# lubridate Library functions are applied 
df.base.hour04$YEAR <- lubridate::year(df.base.hour04$date) # Years component of a date-time

# Hourly data are aggregated into 2-hour [120-min] timestep
df.base.hour04.two <- as.data.frame(rowsum(df.base.hour04[,c(2,3)],
                                           as.integer(gl(nrow(df.base.hour04), 2, nrow(df.base.hour04)))))
df.base.hour04.two$YEAR <- as.integer((df.base.hour04.two$YEAR)/2)

# Hourly data are aggregated into 3-hour timestep [180-min]
df.base.hour04.three <- as.data.frame(rowsum(df.base.hour04[,c(2,3)],
                                             as.integer(gl(nrow(df.base.hour04), 3, nrow(df.base.hour04)))))
df.base.hour04.three$YEAR <- as.integer((df.base.hour04.three$YEAR)/3)

# Hourly data are aggregated into 6-hour timestep [360-min]
df.base.hour04.six <- as.data.frame(rowsum(df.base.hour04[,c(2,3)],
                                           as.integer(gl(nrow(df.base.hour04), 6, nrow(df.base.hour04)))))
df.base.hour04.six$YEAR <- as.integer((df.base.hour04.six$YEAR)/6)

# Hourly data are aggregated into 12-hour timestep [720-min]
df.base.hour04.twelve <- as.data.frame(rowsum(df.base.hour04[,c(2,3)],
                                              as.integer(gl(nrow(df.base.hour04), 12, nrow(df.base.hour04)))))
df.base.hour04.twelve$YEAR <- as.integer((df.base.hour04.twelve$YEAR)/12)

# Hourly data are aggregated into 24-hour timestep [1440-min]
df.base.hour04.daily <- as.data.frame(rowsum(df.base.hour04[,c(2,3)],
                                             as.integer(gl(nrow(df.base.hour04), 24, nrow(df.base.hour04)))))
df.base.hour04.daily$YEAR <- as.integer((df.base.hour04.daily$YEAR)/24)

# Precipitation depths [mm] are transformed to intensities [mm/hour]
df.base.hour04[2] <- df.base.hour04[2]/(60/60)
df.base.hour04.two[1] <- df.base.hour04.two[1]/(120/60)
df.base.hour04.three[1] <- df.base.hour04.three[1]/(180/60)
df.base.hour04.six[1] <- df.base.hour04.six[1]/(360/60)
df.base.hour04.twelve[1] <- df.base.hour04.twelve[1]/(720/60)
df.base.hour04.daily[1] <- df.base.hour04.daily[1]/(1440/60)

# Hourly observed data are aggregated by year
df.aggr.one <- df.base.hour04 %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Hourly observed data are aggregated into 2-hour timestep by year
df.aggr.two <- df.base.hour04.two %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Hourly observed data are aggregated into 3-hour timestep by year
df.aggr.three <- df.base.hour04.three %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Hourly observed data are aggregated into 6-hour timestep by year
df.aggr.six <- df.base.hour04.six %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Hourly observed data are aggregated into 12-hour timestep by year
df.aggr.twelve <- df.base.hour04.twelve %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Hourly observed data are aggregated into 24-hour timestep by year
df.aggr.daily <- df.base.hour04.daily %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Irrelevant rows are deleted
df.aggr.daily <- df.aggr.daily [-c(1), ]

# Sub-hourly and hourly data are c-binded
df.aggr.total <- df.aggr.one
df.aggr.total$x120 <- df.aggr.two$prec
df.aggr.total$x180 <- df.aggr.three$prec
df.aggr.total$x360 <- df.aggr.six$prec
df.aggr.total$x720 <- df.aggr.twelve$prec
df.aggr.total$x1440 <- df.aggr.daily$prec

# ONLY complete years are selected
df.aggr.total <- df.aggr.total[-c(1,26), ]

# Names are assigned to matrix columns
names(df.aggr.total) <- c("year", "min_60", "min_120", "min_180", "min_360", "min_720", "min_1440")

#-------------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: 1998-2023 dataset is completed
#-------------------------------------------------------------------------------------------------------------------

# data.frames are rbinded
df.wide.sub04 <- cbind(df.wide.sub03, df.aggr.total)

# Irrelevant variables are deleted
df.wide.sub04 <- df.wide.sub04[, -c(6, 9, 12)]

# lubridate Library functions are applied 
df_obs_1970_2022$YEAR <- lubridate::year(df_obs_1970_2022$date) # Years component of a date-time

# Hourly observed data are aggregated into 24-hour timestep by year
df.aggr.daily24 <- df_obs_1970_2022 %>%
  group_by(YEAR) %>%
  summarise(prec = max(prec))

# Precipitation depths [mm] are transformed to intensities [mm/hour]
df.aggr.daily24[2] <- df.aggr.daily24[2]/(1440/60)

# data.frames are rbinded
df.wide.sub05 <- rbind(df.wide.sub01, df.wide.sub02, df.wide.sub04)

# Daily totals are added to data.frame
df.wide.sub05$min_1440 <- df.aggr.daily24$prec

# 0 values are replaced with the mean per column in data.frame
df.wide.sub05 <-  df.wide.sub05 %>% mutate(across(.cols = min_5:min_1440,
                                                  .fns = ~ifelse(.x == 0, mean(.x), .x)))

# A AMPs compiled data.frame expressed as intensity [mm/hour] is created
df.amp.compiled.int <- df.wide.sub05

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: 1970-2023 matrices are prepared as volume [mm] and intensity [mm/hour]
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A data.frame is transformed to matrix object
matrix.base.int <- as.matrix(df.wide.sub05)

# Names are assigned to matrix columns
colnames(matrix.base.int) <- c("year", "5",	"10",	"15",	"30", "60", "120", "360", "720", "1440")

# A data.frame is cloned
df.wide.sub06 <- df.wide.sub05

# Intensities [mm/hour] are transformed to Precipitation depths [mm] 
df.wide.sub06[2] <- df.wide.sub06[2]*(5/60)
df.wide.sub06[3] <- df.wide.sub06[3]*(10/60)
df.wide.sub06[4] <- df.wide.sub06[4]*(15/60)
df.wide.sub06[5] <- df.wide.sub06[5]*(30/60)
df.wide.sub06[6] <- df.wide.sub06[6]*(60/60)
df.wide.sub06[7] <- df.wide.sub06[7]*(120/60)
df.wide.sub06[8] <- df.wide.sub06[8]*(360/60)
df.wide.sub06[9] <- df.wide.sub06[9]*(720/60)
df.wide.sub06[10] <- df.wide.sub06[10]*(1440/60)

# A AMPs compiled data.frame expressed as volume [mm] is created
df.amp.compiled.vol <- df.wide.sub06

# A data.frame is transformed to matrix object
matrix.base.vol <- as.matrix(df.wide.sub06)

# Names are assigned to matrix columns
colnames(matrix.base.vol) <- c("year", "5",	"10",	"15",	"30", "60", "120", "360", "720", "1440")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Automatic rain-gauge IDFs development
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# -------------------------------------------------------------------------------------------------
# LOG-PEARSON TYPE 3 DISTRIBUTION [lp3]
# -------------------------------------------------------------------------------------------------
# Generic                LOCATION   SCALE     SHAPE
# IDFtool                mu         sigma     gamma
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# GUMBEL DISTRIBUTION [EV1]
# -------------------------------------------------------------------------------------------------
# Generic                LOCATION   SCALE     SHAPE
# IDFtool                xi         alpha     -----
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# GENERALIZED EXTREME VALUE (GEV) DISTRIBUTION [GEV]
# -------------------------------------------------------------------------------------------------
# Generic                LOCATION   SCALE     SHAPE
# IDFtool                xi         alpha     kappa
# -------------------------------------------------------------------------------------------------

# stat.desc {pastecs} function is requested
df.stat.desc.int <- as.data.frame(round(stat.desc(matrix.base.int[ , ], norm = TRUE), 3))

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: IDF Curves expressed as Intensities [mm/hour]
#-------------------------------------------------------------------------------------------------------------

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.int.lp3 <- IDFCurve(Data = matrix.base.int, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "log.pearson3", M.fit = pdf.fit, 
                             Plot = FALSE, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.int.lp3 <- as.data.frame(IDF.Auto.int.lp3$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.int.lp3$Distribution$`5 min`$Parameters$para

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.int.ev1 <- IDFCurve(Data = matrix.base.int, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "gumbel", M.fit = pdf.fit, 
                             Plot = FALSE, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.int.ev1 <- as.data.frame(IDF.Auto.int.ev1$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.int.ev1$Distribution$`1440 min`$Parameters$para

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.int.gev <- IDFCurve(Data = matrix.base.int, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "gev", M.fit = pdf.fit, 
                             Plot = 3, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.int.gev <- as.data.frame(IDF.Auto.int.gev$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.int.gev$Distribution$`1440 min`$Parameters$para
para.int.gev <- as.data.frame(IDF.Auto.int.gev$Distribution$`1440 min`$Parameters$para)

# A factor variable is created
fit.int.lp3$model <- c("lp3")
fit.int.ev1$model <- c("ev1")
fit.int.gev$model <- c("gev")

# data.frames objects are merged
df.merge.int <- rbind(fit.int.lp3, fit.int.ev1, fit.int.gev)

# Rownames are added as variables
df.merge.int <- tibble::rownames_to_column(df.merge.int, var = "dur")

#-------------------------------------------------------------------------------------------------------------------
# A for-in loop is defined to evaluate KS, AD and CVM TRUE or FALSE
# Kolmogorov-Smirnov (KS) [critical value = 0.20517]
# Anderson-Darling (AD)   [critical value = 2.5018]
# Cramer-von Mises (CVM)  [critical value = 0.22101]
#-------------------------------------------------------------------------------------------------------------------

# List containers are created
df.merge.int$ruleKS <- NULL
df.merge.int$ruleAD <- NULL
df.merge.int$ruleOmega2 <- NULL

# List containers are defined
df.merge.int$ruleKS <- NA
df.merge.int$ruleAD <- NA
df.merge.int$ruleOmega2 <- NA

# Internal counter is defined
counter <- length(df.merge.int$KS)

# ============================
# Inner loop is initialized
# ============================
for(i in 1:counter) {
  
  # Kolmogorov-Smirnov (KS) [critical value = 0.20517]
  if (df.merge.int$KS[i] < 0.20517) {
    df.merge.int$ruleKS[i] = "OK"
  }
  else {
    df.merge.int$ruleKS[i] = "FAILED"
  }
  
  # Anderson-Darling (AD) [critical value = 2.5018]
  if (df.merge.int$AD[i] < 2.5018) {
    df.merge.int$ruleAD[i] = "OK"
  }
  else {
    df.merge.int$ruleAD[i] = "FAILED"
  }
  
  # Cramer-von Mises (CVM) [critical value = 0.22101]
  if (df.merge.int$Omega2[i] < 0.22101) {
    df.merge.int$ruleOmega2[i] = "OK"
  }
  else {
    df.merge.int$ruleOmega2[i] = "FAILED"
  }  
  # ==============================
} # Inner loop is closed
  # ==============================

#-------------------------------------------------------------------------------------------------------------------
# A for-in loop is defined to generate AIC and BIC for comparison purposes
#-------------------------------------------------------------------------------------------------------------------

# data.frame dimensions are requested
dim.array <- dim(df.wide.sub05)

# data.frame dimensions are isolated
dim.array[2]

# List containers are created
list.AIC.lp3 <- NULL
list.BIC.lp3 <- NULL
list.AIC.ev1 <- NULL
list.BIC.ev1 <- NULL
list.AIC.gev <- NULL
list.BIC.gev <- NULL

# List containers are defined
list.AIC.lp3 <- list()
list.BIC.lp3 <- list()
list.AIC.ev1 <- list()
list.BIC.ev1 <- list()
list.AIC.gev <- list()
list.BIC.gev <- list()

# ============================
# Inner loop is initialized
# ============================
for(j in 2:dim.array[2]) {
  
  # MSClaio2008 {nsRFA} is used for Model Selection Criteria
  temp.metric <- MSClaio2008(sample = df.wide.sub05[ , j],
                             dist = c("LP3", "GUMBEL", "GEV"),
                             crit = c("AIC", "BIC"))
  
  list.AIC.lp3[[j]] <- (temp.metric$AIC[1])
  list.BIC.lp3[[j]] <- (temp.metric$BIC[1])
  
  list.AIC.ev1[[j]] <- (temp.metric$AIC[2])
  list.BIC.ev1[[j]] <- (temp.metric$BIC[2])
  
  list.AIC.gev[[j]] <- (temp.metric$AIC[3])
  list.BIC.gev[[j]] <- (temp.metric$BIC[3])
  
  # ============================
} # Inner loop is closed
  # ============================

# An external list of data.frames is created
vec.AIC.lp3 <- as.vector(do.call(rbind, list.AIC.lp3))
vec.BIC.lp3 <- as.vector(do.call(rbind, list.BIC.lp3))

# An external list of data.frames is created
vec.AIC.ev1 <- as.vector(do.call(rbind, list.AIC.ev1))
vec.BIC.ev1 <- as.vector(do.call(rbind, list.BIC.ev1))

# An external list of data.frames is created
vec.AIC.gev <- as.vector(do.call(rbind, list.AIC.gev))
vec.BIC.gev <- as.vector(do.call(rbind, list.BIC.gev))

# Vectors are concatenated
vec.AIC.def <- c(vec.AIC.lp3, vec.AIC.ev1, vec.AIC.gev)
vec.BIC.def <- c(vec.BIC.lp3, vec.BIC.ev1, vec.BIC.gev)

# New variables are added to data.frame
df.merge.int$AIC <- vec.AIC.def #AIC02
df.merge.int$BIC <- vec.BIC.def #BIC02

# Irrelevant variables are deleted
df.merge.int.sel <- df.merge.int[, -c(10,11)]

# ID variables are incorporated
df.merge.int.sel$sta_id <- sta.id 
df.merge.int.sel$rcp_id <- c("hist")

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: IDF Curves expressed as volume [mm]
#-------------------------------------------------------------------------------------------------------------

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.vol.lp3 <- IDFCurve(Data = matrix.base.vol, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "log.pearson3", M.fit = pdf.fit, 
                             Plot = FALSE, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.vol.lp3 <- as.data.frame(IDF.Auto.vol.lp3$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.vol.lp3$Distribution$`1440 min`$Parameters$para

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.vol.ev1 <- IDFCurve(Data = matrix.base.vol, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "gumbel", M.fit = pdf.fit, 
                             Plot = FALSE, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.vol.ev1 <- as.data.frame(IDF.Auto.vol.ev1$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.vol.ev1$Distribution$`1440 min`$Parameters$para

# IDFCurve {IDFtool} function is called and an new object is created
IDF.Auto.vol.gev <- IDFCurve(Data = matrix.base.vol, Station = sta.id,
                             Duration = v.duration, Periods = v.periods, 
                             Type = "gev", M.fit = pdf.fit, 
                             Plot = FALSE, Strategy = 1, logaxe = "", 
                             CI = TRUE, CIpdf = TRUE, iter = 100,
                             goodtest = TRUE, Resolution = 600, 
                             SAVE = FALSE, name = TRUE)

# Goodness-of-fit tests are requested
fit.vol.gev <- as.data.frame(IDF.Auto.vol.gev$Test.fit)

# Distribution parameters are requested for AMP24 ONLY
IDF.Auto.vol.gev$Distribution$`1440 min`$Parameters$para
para.vol.gev.daily <- as.data.frame(IDF.Auto.vol.gev$Distribution$`1440 min`$Parameters$para)

# A factor variable is created
fit.vol.lp3$model <- c("lp3")
fit.vol.ev1$model <- c("ev1")
fit.vol.gev$model <- c("gev")

# data.frames objects are merged
df.merge.vol <- rbind(fit.vol.lp3, fit.vol.ev1, fit.vol.gev)

# Rownames are added as variables
df.merge.vol <- tibble::rownames_to_column(df.merge.vol, var = "dur")

# Distribution parameters are requested
para.vol.gev <- IDF.Auto.vol.gev$Distribution

# A summary of parameters data.frame is created
df.para.vol.gev <- as.data.frame(sapply(para.vol.gev, unlist))

# Names are assigned to data.frame
names(df.para.vol.gev) <- c("min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720", "min_1440")

# ID variables are incorporated
df.para.vol.gev$sta_id <- sta.id 
df.para.vol.gev$rcp_id <- c("hist")

# A summary of parameters data.frame is created
df.para.vol.gev.iso <- df.para.vol.gev[c(1:4), ]

# A summary IDF coefficients data.frame is created
df.para.int.gev.coeff <- as.data.frame(IDF.Auto.int.gev$Models$HIDFUN$Coefficients)

# Rownames are added as variables
df.para.int.gev.coeff <- tibble::rownames_to_column(df.para.int.gev.coeff, var = "period")

# ID variables are incorporated
df.para.int.gev.coeff$sta_id <- sta.id 
df.para.int.gev.coeff$rcp_id <- c("hist")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# AMPs compiled data.frames expressed as intensity [mm/hour] 
# and volume [mm] are requested
View(df.amp.compiled.int)
View(df.amp.compiled.vol)
View(df.stat.desc.int)

# PDFs parameter data.frames are requested
View(df.merge.int.sel)

# GEV geometric parameters as volume [mm] are requested
View(df.para.vol.gev)
View(df.para.vol.gev.iso)

# GEV geometric coefficients as volume [mm] are requested
View(df.para.int.gev.coeff)

# An output file name is concatenated
t.output01 <- paste("BRICK06_matrix_bc_AMP_par_",rcm.id,"_",rcp.id,"_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK06_matrix_bc_AMP_par_",rcm.id,"_",rcp.id,"_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(df.amp.compiled.int,    #1
     df.amp.compiled.vol,    #2
     df.stat.desc.int,       #3
     df.merge.int.sel,       #4
     df.para.vol.gev,        #5
     df.para.vol.gev.iso,    #6
     df.para.int.gev.coeff,  #7
     IDF.Auto.int.gev,       #8
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("df.amp.compiled.int" = df.amp.compiled.int,
                   "df.amp.compiled.vol" = df.amp.compiled.vol,
                   "df.stat.desc.int" = df.stat.desc.int,
                   "df.merge.int.sel" = df.merge.int.sel,
                   "df.para.vol.gev" = df.para.vol.gev,
                   "IDF.Auto.int.gev.coefficients" = df.para.int.gev.coeff,
                   "IDF.Auto.int.gev.intensities" = IDF.Auto.int.gev$Intensities)
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
