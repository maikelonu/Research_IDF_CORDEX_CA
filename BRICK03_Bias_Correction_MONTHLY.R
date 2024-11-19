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
# INFO: This script is intended for the bias correction (BC) of Daily CORDEX-Datasets based on specific 
# Lat/Long position of weather stations
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:

# BRICK01_df_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing dataset vectors of cell-to-point time-series extraction
# date: date as month/day/year
# prec: precipitation [mm/day]
# grid_lon: Weather station longitude [degrees]
# grid_lat: Weather station latitude [degrees]
# sta_id: Weather station ID
# rcm_id: CORDEX member ID
# crtm05_x: CRTM05 X [m]
# crtm05_y: CRTM05 X [y]

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

# BRICK03_matrix_bc_24h_par_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing:
# amp_24hr_sim_comp            #1: AMP24-BC prec.[mm/day] for 2006-2099
# amp_24hr_hist_comp           #2  AMP24-BC prec.[mm/day] for 1991-2005
# amp_24hr_arch_comp           #3  AMP24-BC prec.[mm/day] for 1970-2005
# amp_24hrAMP_hist_comp        #4  Monthly AMP24-BC prec.[mm/day] for 1991-2005
# df_amp_24hr_sim_duplicated   #5  AMP24-BC Number of duplicates for 2006-2099
# df_amp_24hr_hist_duplicated  #6  AMP24-BC Number of duplicates for 1991-2005
# df_amp_24hr_arch_duplicated  #7  AMP24-BC Number of duplicates for 1970-2005
# df.base.desc                 #8  Stat. parameters for sim hist and arch
# df.perf.desc.24hrAMP         #9  AMP24-BC ModelMetrics perf. for 1991-2005
# df.perf.desc.daily           #10 Daily ModelMetrics perf. for 1991-2005
# df.perf.desc.daily.baseline  #11 Daily ModelMetrics perf. for baseline 1971-2000
# df.perf.climdex              #12 Daily ClimDex perf. for 1991-2005
# df.perf.climdex.baseline     #13 Daily ClimDex perf. for baseline 1971-2000
# df.rx1day_YEAR.long          #14 AMP24-BC rx1day for 1991-2005
# df.rx1day_YEAR_hist_fut_comp #15 AMP24-BC rx1day for 2005-2099
# df.ecdf.total.aggregated     #16 BC Monthly Averages ecdf for for 1991-2005
# df.ecdf.total                #17 BC Daily Averages ecdf for for 1991-2005
#-------------------------------------------------------------------------------------------------------------------

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
require(precintcon)
require(nsRFA)
require(imputeTS)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: External sources are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# External scripts are sourced
source("BLOCK_biasCorrection1D.R")
#source("BLOCK_georeference.R")
#source("BLOCK_RCM_Extract.R")
#source("BLOCK_CRAN.R")
#source("BLOCK_CLC.R")
#source("BLOCK_OBS.R")

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Pre-processed binary datasets are loaded
#-------------------------------------------------------------------------------------------------------------

# RData objects are loaded from BRICK01
load("BRICK01_df_HadGEM2_RegCM47_RCP85_AIJS.RData")

# RData objects are loaded from BRICK02
load("BRICK02_df_obs_AIJS.RData")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX + Weather Station Metada + Input Data
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

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Preprocessed binary datasets are modified
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

# BEWARE!! irrelevant variables are deleted from data.frames
df_obs_1970_1990 <- df_obs_1970_1990[,c(1,2)]
df_obs_1970_2005 <- df_obs_1970_2005[,c(1,2)]
df_obs_1970_2022 <- df_obs_1970_2022[,c(1,2)]
df_obs_1991_2005 <- df_obs_1991_2005[,c(1,2)]
df.export.rcm01 <- df.export.rcm01[,c(1,2)]
df.export.rcm02 <- df.export.rcm02[,c(1,2)]

# A subset observational dataset is created for the 1971-2000 baseline
df_obs_1971_2000 <- df_obs_1970_2005[df_obs_1970_2005$date >= "1971-01-01" & df_obs_1970_2005$date <= "2000-12-31", ]

# {lubridate} functions are applied to create new variables monthly
df_obs_1970_1990$MONTH <- lubridate::month(df_obs_1970_1990$date)
df_obs_1970_2005$MONTH <- lubridate::month(df_obs_1970_2005$date)
df_obs_1991_2005$MONTH <- lubridate::month(df_obs_1991_2005$date)
df_obs_1971_2000$MONTH <- lubridate::month(df_obs_1971_2000$date)
df.export.rcm01$MONTH <- lubridate::month(df.export.rcm01$date)
df.export.rcm02$MONTH <- lubridate::month(df.export.rcm02$date)

# {lubridate} functions are applied to create new variables yearly
df_obs_1970_1990$YEAR <- lubridate::year(df_obs_1970_1990$date)
df_obs_1970_2005$YEAR <- lubridate::year(df_obs_1970_2005$date)
df_obs_1991_2005$YEAR <- lubridate::year(df_obs_1991_2005$date)
df_obs_1971_2000$YEAR <- lubridate::year(df_obs_1971_2000$date)
df.export.rcm01$YEAR <- lubridate::year(df.export.rcm01$date)
df.export.rcm02$YEAR <- lubridate::year(df.export.rcm02$date)

# A historical numerical monthly simulation data.frame is created
df.base.future <- df.export.rcm02[df.export.rcm02$date >= "2006-01-01" & df.export.rcm02$date <= "2099-12-31", ]
df.hindcast.val <- df.export.rcm01[df.export.rcm01$date >= "1991-01-01" & df.export.rcm01$date <= "2005-12-31", ]
df.hindcast.cal <- df.export.rcm01[df.export.rcm01$date >= "1970-01-01" & df.export.rcm01$date <= "1990-12-31", ]
df.hindcast.tot <- df.export.rcm01[df.export.rcm01$date >= "1970-01-01" & df.export.rcm01$date <= "2005-12-31", ]
df.hindcast.baseline <- df.export.rcm01[df.export.rcm01$date >= "1971-01-01" & df.export.rcm01$date <= "2000-12-31", ]

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Execution of Bias-Correction
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A maximum historical threshold factor is defined
thres.factor.max <- (max(df_obs_1970_2022$prec))*as.numeric(attri[13, 2])

# A list container is erased
list.cont_hist_rqm_qmap <- NULL
list.cont_hist_gpm_dR <- NULL
list.cont_hist_eqm_dR <- NULL
list.cont_hist_dqm_dR <- NULL
list.cont_hist_qdm_dR <- NULL
list.cont_sim_rqm_qmap <- NULL
list.cont_sim_gpm_dR <- NULL
list.cont_sim_eqm_dR <- NULL
list.cont_sim_dqm_dR <- NULL
list.cont_sim_qdm_dR <- NULL
list.cont_arch_rqm_qmap <- NULL
list.cont_arch_gpm_dR <- NULL
list.cont_arch_eqm_dR <- NULL
list.cont_arch_dqm_dR <- NULL
list.cont_arch_qdm_dR <- NULL
list.int.vect.obs.fut <- NULL
list.int.vect.hin.fut <- NULL
list.int.vect.sim.fut <- NULL
list.int.vect.obs.hist <-NULL
list.int.vect.hin.hist <- NULL
list.int.vect.sim.hist <- NULL

# A list container is defined
list.cont_hist_rqm_qmap <- list()
list.cont_hist_gpm_dR <- list()
list.cont_hist_eqm_dR <- list()
list.cont_hist_dqm_dR <- list()
list.cont_hist_qdm_dR <- list()
list.cont_sim_rqm_qmap <- list()
list.cont_sim_gpm_dR <- list()
list.cont_sim_eqm_dR <- list()
list.cont_sim_dqm_dR <- list()
list.cont_sim_qdm_dR <- list()
list.cont_arch_rqm_qmap <- list()
list.cont_arch_gpm_dR <- list()
list.cont_arch_eqm_dR <- list()
list.cont_arch_dqm_dR <- list()
list.cont_arch_qdm_dR <- list()
list.int.vect.obs.fut <- list()
list.int.vect.hin.fut <- list()
list.int.vect.sim.fut <- list()
list.int.vect.obs.hist <-list()
list.int.vect.hin.hist <- list()
list.int.vect.sim.hist <- list()

# ==============================
# Outermost loop is initialized
# ==============================
for (i in 1:12) { # length = 12: months
  
  # Data containers are erased
  bc_sim_eqm_dR <- NULL
  bc_hist_eqm_dR <- NULL
  bc_arch_eqm_dR <- NULL
  bc_sim_gpm_dR <- NULL
  bc_hist_gpm_dR <- NULL
  bc_arch_gpm_dR <- NULL
  bc_sim_dqm_dR <- NULL
  bc_hist_dqm_dR <- NULL
  bc_arch_dqm_dR <- NULL
  bc_sim_qdm_dR <- NULL
  bc_hist_qdm_dR <- NULL
  bc_arch_qdm_dR <- NULL
  qm.fit.sim <- NULL
  bc_sim_rqm_qmap <- NULL
  qm.fit.hist <- NULL
  bc_hist_rqm_qmap <- NULL
  qm.fit.arch <- NULL
  bc_arch_rqm_qmap <- NULL
  
  # A future simulated monthly subset is executed
  int.vect.obs.fut <- subset(df_obs_1970_2005, MONTH==i)   # 1970-2005
  int.vect.hin.fut <- subset(df.hindcast.tot, MONTH==i)    # 1970-2005
  int.vect.sim.fut <- subset(df.base.future, MONTH==i)     # 2006-2099
  
  # A historical monthly subset is executed
  int.vect.obs.hist <- subset(df_obs_1970_1990, MONTH==i)  # 1970-1990
  int.vect.hin.hist <- subset(df.hindcast.cal, MONTH==i)   # 1970-1990
  int.vect.sim.hist <- subset(df.hindcast.val, MONTH==i)   # 1991-2005
  
  # -------------------------------------------------------------------------------------------------------------------
  # biasCorrect {hyfo} is used for future bias correction
  bc_sim_eqm_dR <- biasCorrect(obs = int.vect.obs.fut,      # 1970-2005
                               hindcast = int.vect.hin.fut, # 1970-2005
                               frc =  int.vect.sim.fut,     # 2006-2099
                               preci = TRUE,
                               prThreshold = umbral.rain,
                               method = "eqm",
                               extrapolate = "constant")
  
  # Precipitation field is transformed to vector
  bc_sim_eqm_dR <- as.numeric(bc_sim_eqm_dR$prec)
  
  # NAs OR over-threshold values are replaced
  bc_sim_eqm_dR <- ifelse(bc_sim_eqm_dR > thres.factor.max, thres.factor.max, bc_sim_eqm_dR)
  bc_sim_eqm_dR <- ifelse(is.na(bc_sim_eqm_dR), 0, bc_sim_eqm_dR)
  
  # biasCorrect {hyfo} is used for for historical bias correction
  bc_hist_eqm_dR <- biasCorrect(obs = int.vect.obs.hist,      # 1970-1990
                                hindcast = int.vect.hin.hist, # 1970-1990
                                frc =  int.vect.sim.hist,     # 1991-2005
                                preci = TRUE,
                                prThreshold = umbral.rain,
                                method = "eqm",
                                extrapolate = "constant")
  
  # Precipitation field is transformed to vector
  bc_hist_eqm_dR <- as.numeric(bc_hist_eqm_dR$prec)
  
  # NAs OR over-threshold values are replaced
  bc_hist_eqm_dR <- ifelse(bc_hist_eqm_dR > thres.factor.max, thres.factor.max, bc_hist_eqm_dR)
  bc_hist_eqm_dR <- ifelse(is.na(bc_hist_eqm_dR), 0, bc_hist_eqm_dR)
  
  # biasCorrect {hyfo} is used for for historical bias correction
  bc_arch_eqm_dR <- biasCorrect(obs = int.vect.obs.fut,      # 1970-2005
                                hindcast = int.vect.hin.fut, # 1970-2005
                                frc =  int.vect.hin.fut,     # 1970-2005
                                preci = TRUE,
                                prThreshold = umbral.rain,
                                method = "eqm",
                                extrapolate = "constant")
  
  # Precipitation field is transformed to vector
  bc_arch_eqm_dR <- as.numeric(bc_arch_eqm_dR$prec)
  
  # NAs OR over-threshold values are replaced
  bc_arch_eqm_dR <- ifelse(bc_arch_eqm_dR > thres.factor.max, thres.factor.max, bc_arch_eqm_dR)
  bc_arch_eqm_dR <- ifelse(is.na(bc_arch_eqm_dR), 0, bc_arch_eqm_dR)
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # biasCorrection1D {downscaleR} is used for future bias correction
  bc_sim_gpm_dR <- gpqm(o = int.vect.obs.fut[, 2],         # 1970-2005
                        p = int.vect.hin.fut[, 2],         # 1970-2005
                        s = int.vect.sim.fut[, 2],         # 2006-2099
                        precip = TRUE,
                        pr.threshold = umbral.rain,
                        theta = p.theta)
  
  # NAs OR over-threshold values are replaced
  bc_sim_gpm_dR <- ifelse(bc_sim_gpm_dR > thres.factor.max, thres.factor.max, bc_sim_gpm_dR)
  bc_sim_gpm_dR <- ifelse(is.na(bc_sim_gpm_dR), 0, bc_sim_gpm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_hist_gpm_dR <- gpqm(o = int.vect.obs.hist[, 2],       # 1970-1990
                         p = int.vect.hin.hist[, 2],       # 1970-1990
                         s = int.vect.sim.hist[, 2],       # 1991-2005
                         precip = TRUE,
                         pr.threshold = umbral.rain,
                         theta = p.theta)
  
  # NAs OR over-threshold values are replaced
  bc_hist_gpm_dR <- ifelse(bc_hist_gpm_dR > thres.factor.max, thres.factor.max, bc_hist_gpm_dR)
  bc_hist_gpm_dR <- ifelse(is.na(bc_hist_gpm_dR), 0, bc_hist_gpm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_arch_gpm_dR <- gpqm(o = int.vect.obs.fut[, 2],        # 1970-2005
                         p = int.vect.hin.fut[, 2],        # 1970-2005
                         s = int.vect.hin.fut[, 2],        # 1970-2005
                         precip = TRUE,
                         pr.threshold = umbral.rain,
                         theta = p.theta)
  
  # NAs OR over-threshold values are replaced
  bc_arch_gpm_dR <- ifelse(bc_arch_gpm_dR > thres.factor.max, thres.factor.max, bc_arch_gpm_dR)
  bc_arch_gpm_dR <- ifelse(is.na(bc_arch_gpm_dR), 0, bc_arch_gpm_dR)
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # biasCorrection1D {downscaleR} is used for future bias correction
  bc_sim_dqm_dR <- dqm(o = int.vect.obs.fut[, 2],          # 1970-2005
                       p = int.vect.hin.fut[, 2],          # 1970-2005
                       s = int.vect.sim.fut[, 2],          # 2006-2099
                       precip = TRUE,
                       pr.threshold = umbral.rain,
                       n.quantiles = NULL,
                       detrend = TRUE)
  
  # NAs OR over-threshold values are replaced
  bc_sim_dqm_dR <- ifelse(bc_sim_dqm_dR > thres.factor.max, thres.factor.max, bc_sim_dqm_dR)
  bc_sim_dqm_dR <- ifelse(is.na(bc_sim_dqm_dR), 0, bc_sim_dqm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_hist_dqm_dR <- dqm(o = int.vect.obs.hist[, 2],        # 1970-1990
                        p = int.vect.hin.hist[, 2],        # 1970-1990
                        s = int.vect.sim.hist[, 2],        # 1991-2005
                        precip = TRUE,
                        pr.threshold = umbral.rain,
                        n.quantiles = NULL,
                        detrend = TRUE)
  
  # NAs OR over-threshold values are replaced
  bc_hist_dqm_dR <- ifelse(bc_hist_dqm_dR > thres.factor.max, thres.factor.max, bc_hist_dqm_dR)
  bc_hist_dqm_dR <- ifelse(is.na(bc_hist_dqm_dR), 0, bc_hist_dqm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_arch_dqm_dR <- dqm(o = int.vect.obs.fut[, 2],         # 1970-2005
                        p = int.vect.hin.fut[, 2],         # 1970-2005
                        s = int.vect.hin.fut[, 2],         # 1970-2005
                        precip = TRUE,
                        pr.threshold = umbral.rain,
                        n.quantiles = NULL,
                        detrend = TRUE)
  
  # NAs OR over-threshold values are replaced
  bc_arch_dqm_dR <- ifelse(bc_arch_dqm_dR > thres.factor.max, thres.factor.max, bc_arch_dqm_dR)
  bc_arch_dqm_dR <- ifelse(is.na(bc_arch_dqm_dR), 0, bc_arch_dqm_dR)
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # biasCorrection1D {downscaleR} is used for future bias correction
  bc_sim_qdm_dR <- qdm(o = int.vect.obs.fut[, 2],          # 1970-2005
                       p = int.vect.hin.fut[, 2],          # 1970-2005
                       s = int.vect.sim.fut[, 2],          # 2006-2099
                       precip = TRUE,
                       pr.threshold = umbral.rain,
                       n.quantiles = NULL,
                       jitter.factor = 0.01)
  
  # NAs OR over-threshold values are replaced
  bc_sim_qdm_dR <- ifelse(bc_sim_qdm_dR > thres.factor.max, thres.factor.max, bc_sim_qdm_dR)
  bc_sim_qdm_dR <- ifelse(is.na(bc_sim_qdm_dR), 0, bc_sim_qdm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_hist_qdm_dR <- qdm(o = int.vect.obs.hist[, 2],        # 1970-1990
                        p = int.vect.hin.hist[, 2],        # 1970-1990
                        s = int.vect.sim.hist[, 2],        # 1991-2005
                        precip = TRUE,
                        pr.threshold = umbral.rain,
                        n.quantiles = NULL,
                        jitter.factor = 0.01)
  
  # NAs OR over-threshold values are replaced
  bc_hist_qdm_dR <- ifelse(bc_hist_qdm_dR > thres.factor.max, thres.factor.max, bc_hist_qdm_dR)
  bc_hist_qdm_dR <- ifelse(is.na(bc_hist_qdm_dR), 0, bc_hist_qdm_dR)
  
  # biasCorrection1D {downscaleR} is used for historical bias correction
  bc_arch_qdm_dR <- qdm(o = int.vect.obs.fut[, 2],         # 1970-2005
                        p = int.vect.hin.fut[, 2],         # 1970-2005
                        s = int.vect.hin.fut[, 2],         # 1970-2005
                        precip = TRUE,
                        pr.threshold = umbral.rain,
                        n.quantiles = NULL,
                        jitter.factor = 0.01)
  
  # NAs OR over-threshold values are replaced
  bc_arch_qdm_dR <- ifelse(bc_arch_qdm_dR > thres.factor.max, thres.factor.max, bc_arch_qdm_dR)
  bc_arch_qdm_dR <- ifelse(is.na(bc_arch_qdm_dR), 0, bc_arch_qdm_dR)
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # fitQmap {qmap} is used for future bias correction
  qm.fit.sim <- fitQmapRQUANT(int.vect.obs.fut[, 2],       # 1970-2005 fitQmaprqmLIN/fitQmapRQUANT
                              int.vect.hin.fut[, 2],       # 1970-2005
                              wet.day = umbral.rain,
                              qstep = 0.01)
  
  bc_sim_rqm_qmap <- doQmapRQUANT(int.vect.sim.fut[, 2],   # 2006-2099 doQmaprqmLIN/doQmapRQUANT
                                  qm.fit.sim,
                                  type="linear2")
  
  # NAs OR over-threshold values are replaced
  bc_sim_rqm_qmap <- ifelse(bc_sim_rqm_qmap > thres.factor.max, thres.factor.max, bc_sim_rqm_qmap)
  bc_sim_rqm_qmap <- ifelse(is.na(bc_sim_rqm_qmap), 0, bc_sim_rqm_qmap)
  
  #fitQmap {qmap} is used for historical bias correction
  qm.fit.hist <- fitQmapRQUANT(int.vect.obs.hist[, 2],     # 1970-1990 fitQmaprqmLIN/fitQmapRQUANT
                               int.vect.hin.hist[, 2],     # 1970-1990
                               wet.day = umbral.rain,
                               qstep = 0.01)
  
  bc_hist_rqm_qmap <- doQmapRQUANT(int.vect.sim.hist[, 2], # 1991-2005 doQmaprqmLIN/doQmapRQUANT
                                   qm.fit.hist,
                                   type="linear2")
  
  # NAs OR over-threshold values are replaced
  bc_hist_rqm_qmap <- ifelse(bc_hist_rqm_qmap > thres.factor.max, thres.factor.max, bc_hist_rqm_qmap)
  bc_hist_rqm_qmap <- ifelse(is.na(bc_hist_rqm_qmap), 0, bc_hist_rqm_qmap)
  
  #fitQmap {qmap} is used for historical bias correction
  qm.fit.arch <- fitQmapRQUANT(int.vect.obs.fut[, 2],      # 1970-2005 fitQmaprqmLIN/fitQmapRQUANT
                               int.vect.hin.fut[, 2],      # 1970-2005
                               wet.day = umbral.rain,
                               qstep = 0.01)
  
  bc_arch_rqm_qmap <- doQmapRQUANT(int.vect.hin.fut[, 2],  # 1970-2005 doQmaprqmLIN/doQmapRQUANT
                                   qm.fit.arch,
                                   type="linear2")
  
  # NAs OR over-threshold values are replaced
  bc_arch_rqm_qmap <- ifelse(bc_arch_rqm_qmap > thres.factor.max, thres.factor.max, bc_arch_rqm_qmap)
  bc_arch_rqm_qmap <- ifelse(is.na(bc_arch_rqm_qmap), 0, bc_arch_rqm_qmap)
  # -------------------------------------------------------------------------------------------------------------------
  
  # A list of "prec" data.frames is created
  list.cont_hist_rqm_qmap[[i]] <- bc_hist_rqm_qmap
  list.cont_hist_gpm_dR[[i]] <- bc_hist_gpm_dR
  list.cont_hist_eqm_dR[[i]] <- bc_hist_eqm_dR
  list.cont_hist_dqm_dR[[i]] <- bc_hist_dqm_dR
  list.cont_hist_qdm_dR[[i]] <- bc_hist_qdm_dR
  list.cont_sim_rqm_qmap[[i]] <-bc_sim_rqm_qmap
  list.cont_sim_gpm_dR[[i]] <- bc_sim_gpm_dR
  list.cont_sim_eqm_dR[[i]] <- bc_sim_eqm_dR
  list.cont_sim_dqm_dR[[i]] <- bc_sim_dqm_dR
  list.cont_sim_qdm_dR[[i]] <- bc_sim_qdm_dR
  list.cont_arch_rqm_qmap[[i]] <- bc_arch_rqm_qmap
  list.cont_arch_gpm_dR[[i]] <- bc_arch_gpm_dR
  list.cont_arch_eqm_dR[[i]] <- bc_arch_eqm_dR
  list.cont_arch_dqm_dR[[i]] <- bc_arch_dqm_dR
  list.cont_arch_qdm_dR[[i]] <- bc_arch_qdm_dR
  
  # A list of "date" data.frames is created
  list.int.vect.obs.fut[[i]] <- int.vect.obs.fut$date
  list.int.vect.hin.fut[[i]] <- int.vect.hin.fut$date
  list.int.vect.sim.fut[[i]] <- int.vect.sim.fut$date
  list.int.vect.obs.hist[[i]] <-int.vect.obs.hist$date
  list.int.vect.hin.hist[[i]] <- int.vect.hin.hist$date
  list.int.vect.sim.hist[[i]] <- int.vect.sim.hist$date
  
  # ==============================
} # Outermost loop is closed
  # ==============================

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Comparison of future daily projections for each BC method
#-------------------------------------------------------------------------------------------------------------

# An external list of "date" data.frames is created
date.list.int.vect.obs.fut <- ldply(list.int.vect.obs.fut, data.frame)    # 13149 obs
date.list.int.vect.hin.fut <- ldply(list.int.vect.hin.fut, data.frame)    # 13149 obs
date.list.int.vect.sim.fut <- ldply(list.int.vect.sim.fut, data.frame)    # 34333 obs
date.list.int.vect.obs.hist <- ldply(list.int.vect.obs.hist, data.frame)  # 7670 obs
date.list.int.vect.hin.hist <- ldply(list.int.vect.hin.hist, data.frame)  # 7670 obs
date.list.int.vect.sim.hist <- ldply(list.int.vect.sim.hist, data.frame)  # 5479 obs

# An external list of validation data.frames is created
df.hist.temp.rqm <- ldply(list.cont_hist_rqm_qmap, data.frame)            # 5479 obs
df.hist.temp.gpm <- ldply(list.cont_hist_gpm_dR, data.frame)              # 5479 obs
df.hist.temp.eqm <- ldply(list.cont_hist_eqm_dR, data.frame)              # 5479 obs
df.hist.temp.dqm <- ldply(list.cont_hist_dqm_dR, data.frame)              # 5479 obs
df.hist.temp.qdm <- ldply(list.cont_hist_qdm_dR, data.frame)              # 5479 obs
df.sim.temp.rqm <- ldply(list.cont_sim_rqm_qmap, data.frame)              # 34333 obs
df.sim.temp.gpm <- ldply(list.cont_sim_gpm_dR, data.frame)                # 34333 obs
df.sim.temp.eqm <- ldply(list.cont_sim_eqm_dR, data.frame)                # 34333 obs
df.sim.temp.dqm <- ldply(list.cont_sim_dqm_dR, data.frame)                # 34333 obs
df.sim.temp.qdm <- ldply(list.cont_sim_qdm_dR, data.frame)                # 34333 obs
df.arch.temp.rqm <- ldply(list.cont_arch_rqm_qmap, data.frame)            # 7670 obs
df.arch.temp.gpm <- ldply(list.cont_arch_gpm_dR, data.frame)              # 7670 obs
df.arch.temp.eqm <- ldply(list.cont_arch_eqm_dR, data.frame)              # 7670 obs
df.arch.temp.dqm <- ldply(list.cont_arch_dqm_dR, data.frame)              # 7670 obs
df.arch.temp.qdm <- ldply(list.cont_arch_qdm_dR, data.frame)              # 7670 obs

# Variables "prec" and "date" are binded
df.hist.reintegrated.rqm <- cbind(date.list.int.vect.sim.hist, df.hist.temp.rqm)  # 5479 obs
df.hist.reintegrated.gpm <- cbind(date.list.int.vect.sim.hist, df.hist.temp.gpm)  # 5479 obs
df.hist.reintegrated.eqm <- cbind(date.list.int.vect.sim.hist, df.hist.temp.eqm)  # 5479 obs
df.hist.reintegrated.dqm <- cbind(date.list.int.vect.sim.hist, df.hist.temp.dqm)  # 5479 obs
df.hist.reintegrated.qdm <- cbind(date.list.int.vect.sim.hist, df.hist.temp.qdm)  # 5479 obs
df.sim.reintegrated.rqm <- cbind(date.list.int.vect.sim.fut, df.sim.temp.rqm)     # 34333 obs
df.sim.reintegrated.gpm <- cbind(date.list.int.vect.sim.fut, df.sim.temp.gpm)     # 34333 obs
df.sim.reintegrated.eqm <- cbind(date.list.int.vect.sim.fut, df.sim.temp.eqm)     # 34333 obs
df.sim.reintegrated.dqm <- cbind(date.list.int.vect.sim.fut, df.sim.temp.dqm)     # 34333 obs
df.sim.reintegrated.qdm <- cbind(date.list.int.vect.sim.fut, df.sim.temp.qdm)     # 34333 obs
df.arch.reintegrated.rqm <- cbind(date.list.int.vect.obs.fut, df.arch.temp.rqm)   # 7670 obs
df.arch.reintegrated.gpm <- cbind(date.list.int.vect.obs.fut, df.arch.temp.gpm)   # 7670 obs
df.arch.reintegrated.eqm <- cbind(date.list.int.vect.obs.fut, df.arch.temp.eqm)   # 7670 obs
df.arch.reintegrated.dqm <- cbind(date.list.int.vect.obs.fut, df.arch.temp.dqm)   # 7670 obs
df.arch.reintegrated.qdm <- cbind(date.list.int.vect.obs.fut, df.arch.temp.qdm)   # 7670 obs

# Variable names are defined
names(df.hist.reintegrated.rqm) <- c("dates", "prec")
names(df.hist.reintegrated.gpm) <- c("dates", "prec")
names(df.hist.reintegrated.eqm) <- c("dates", "prec")
names(df.hist.reintegrated.dqm) <- c("dates", "prec")
names(df.hist.reintegrated.qdm) <- c("dates", "prec")
names(df.sim.reintegrated.rqm) <- c("dates", "prec")
names(df.sim.reintegrated.gpm) <- c("dates", "prec")
names(df.sim.reintegrated.eqm) <- c("dates", "prec")
names(df.sim.reintegrated.dqm) <- c("dates", "prec")
names(df.sim.reintegrated.qdm) <- c("dates", "prec")
names(df.arch.reintegrated.rqm) <- c("dates", "prec")
names(df.arch.reintegrated.gpm) <- c("dates", "prec")
names(df.arch.reintegrated.eqm) <- c("dates", "prec")
names(df.arch.reintegrated.dqm) <- c("dates", "prec")
names(df.arch.reintegrated.qdm) <- c("dates", "prec")

# data.frames are sorted by "date" variable
bc_hist_rqm_qmap <- df.hist.reintegrated.rqm[order(df.hist.reintegrated.rqm$dates),]
bc_hist_gpm_dR <- df.hist.reintegrated.gpm[order(df.hist.reintegrated.gpm$dates),]
bc_hist_eqm_dR <- df.hist.reintegrated.eqm[order(df.hist.reintegrated.eqm$dates),]
bc_hist_dqm_dR <- df.hist.reintegrated.dqm[order(df.hist.reintegrated.dqm$dates),]
bc_hist_qdm_dR <- df.hist.reintegrated.qdm[order(df.hist.reintegrated.qdm$dates),]

# data.frames are sorted by "date" variable
bc_sim_rqm_qmap <- df.sim.reintegrated.rqm[order(df.sim.reintegrated.rqm$dates),]
bc_sim_gpm_dR <- df.sim.reintegrated.gpm[order(df.sim.reintegrated.gpm$dates),]
bc_sim_eqm_dR <- df.sim.reintegrated.eqm[order(df.sim.reintegrated.eqm$dates),]
bc_sim_dqm_dR <- df.sim.reintegrated.dqm[order(df.sim.reintegrated.dqm$dates),]
bc_sim_qdm_dR <- df.sim.reintegrated.qdm[order(df.sim.reintegrated.qdm$dates),]

# data.frames are sorted by "date" variable
bc_arch_rqm_qmap <- df.arch.reintegrated.rqm[order(df.arch.reintegrated.rqm$dates),]
bc_arch_gpm_dR <- df.arch.reintegrated.gpm[order(df.arch.reintegrated.gpm$dates),]
bc_arch_eqm_dR <- df.arch.reintegrated.eqm[order(df.arch.reintegrated.eqm$dates),]
bc_arch_dqm_dR <- df.arch.reintegrated.dqm[order(df.arch.reintegrated.dqm$dates),]
bc_arch_qdm_dR <- df.arch.reintegrated.qdm[order(df.arch.reintegrated.qdm$dates),]

# {lubridate} functions are applied to create new columns containing:
bc_sim_eqm_dR$YEAR <- lubridate::year(bc_sim_eqm_dR$date)                 # 2006-2099
bc_sim_rqm_qmap$YEAR <- lubridate::year(bc_sim_rqm_qmap$date)             # 2006-2099
bc_sim_gpm_dR$YEAR <- lubridate::year(bc_sim_gpm_dR$date)                 # 2006-2099
bc_sim_dqm_dR$YEAR <- lubridate::year(bc_sim_dqm_dR$date)                 # 2006-2099
bc_sim_qdm_dR$YEAR <- lubridate::year(bc_sim_qdm_dR$date)                 # 2006-2099

# {lubridate} functions are applied to create new columns containing:
bc_hist_eqm_dR$YEAR <- lubridate::year(bc_hist_eqm_dR$date)               # 1991-2005
bc_hist_rqm_qmap$YEAR <- lubridate::year(bc_hist_rqm_qmap$date)           # 1991-2005
bc_hist_gpm_dR$YEAR <- lubridate::year(bc_hist_gpm_dR$date)               # 1991-2005
bc_hist_dqm_dR$YEAR <- lubridate::year(bc_hist_dqm_dR$date)               # 1991-2005
bc_hist_qdm_dR$YEAR <- lubridate::year(bc_hist_qdm_dR$date)               # 1991-2005

# {lubridate} functions are applied to create new columns containing:
bc_arch_eqm_dR$YEAR <- lubridate::year(bc_arch_eqm_dR$date)               # 1970-1990
bc_arch_rqm_qmap$YEAR <- lubridate::year(bc_arch_rqm_qmap$date)           # 1970-1990
bc_arch_gpm_dR$YEAR <- lubridate::year(bc_arch_gpm_dR$date)               # 1970-1990
bc_arch_dqm_dR$YEAR <- lubridate::year(bc_arch_dqm_dR$date)               # 1970-1990
bc_arch_qdm_dR$YEAR <- lubridate::year(bc_arch_qdm_dR$date)               # 1970-1990

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Derivation of daily AMP24 per BC method 
#-------------------------------------------------------------------------------------------------------------

# Daily AMP24 observed data are aggregated by year
amp_24hr_sim_raw <- as.data.frame(tapply(df.base.future$prec, df.base.future$YEAR, max))       # 2006-2099
amp_24hr_sim_eqm <- as.data.frame(tapply(bc_sim_eqm_dR$prec, bc_sim_eqm_dR$YEAR, max))         # 2006-2099
amp_24hr_sim_rqm <- as.data.frame(tapply(bc_sim_rqm_qmap$prec, bc_sim_rqm_qmap$YEAR, max))     # 2006-2099
amp_24hr_sim_gpm <- as.data.frame(tapply(bc_sim_gpm_dR$prec, bc_sim_gpm_dR$YEAR, max))         # 2006-2099
amp_24hr_sim_dqm <- as.data.frame(tapply(bc_sim_dqm_dR$prec, bc_sim_dqm_dR$YEAR, max))         # 2006-2099
amp_24hr_sim_qdm <- as.data.frame(tapply(bc_sim_qdm_dR$prec, bc_sim_qdm_dR$YEAR, max))         # 2006-2099

# Daily AMP24 observed data are aggregated by year
amp_24hr_hist_raw <- as.data.frame(tapply(df.hindcast.val$prec, df.hindcast.val$YEAR, max))    # 1990-2005
amp_24hr_hist_eqm <- as.data.frame(tapply(bc_hist_eqm_dR$prec, bc_hist_eqm_dR$YEAR, max))      # 1990-2005
amp_24hr_hist_rqm <- as.data.frame(tapply(bc_hist_rqm_qmap$prec, bc_hist_rqm_qmap$YEAR, max))  # 1990-2005
amp_24hr_hist_gpm <- as.data.frame(tapply(bc_hist_gpm_dR$prec, bc_hist_gpm_dR$YEAR, max))      # 1990-2005
amp_24hr_hist_dqm <- as.data.frame(tapply(bc_hist_dqm_dR$prec, bc_hist_dqm_dR$YEAR, max))      # 1990-2005
amp_24hr_hist_qdm <- as.data.frame(tapply(bc_hist_qdm_dR$prec, bc_hist_qdm_dR$YEAR, max))      # 1990-2005

# Daily AMP24 observed data are aggregated by year
amp_24hr_arch_raw <- as.data.frame(tapply(df.hindcast.tot$prec, df.hindcast.tot$YEAR, max))    # 1970-2005
amp_24hr_arch_eqm <- as.data.frame(tapply(bc_arch_eqm_dR$prec, bc_arch_eqm_dR$YEAR, max))      # 1970-2005
amp_24hr_arch_rqm <- as.data.frame(tapply(bc_arch_rqm_qmap$prec, bc_arch_rqm_qmap$YEAR, max))  # 1970-2005
amp_24hr_arch_gpm <- as.data.frame(tapply(bc_arch_gpm_dR$prec, bc_arch_gpm_dR$YEAR, max))      # 1970-2005
amp_24hr_arch_dqm <- as.data.frame(tapply(bc_arch_dqm_dR$prec, bc_arch_dqm_dR$YEAR, max))      # 1970-2005
amp_24hr_arch_qdm <- as.data.frame(tapply(bc_arch_qdm_dR$prec, bc_arch_qdm_dR$YEAR, max))      # 1970-2005

# Variable names are added to data.frame 
amp_24hr_sim_comp_dates <- rownames(amp_24hr_sim_eqm)   # for future projections
amp_24hr_hist_comp_dates <- rownames(amp_24hr_hist_eqm) # for historical projections
amp_24hr_arch_comp_dates <- rownames(amp_24hr_arch_eqm) # for archival projections

# A comparison data.frame is created
amp_24hr_sim_comp <- data.frame(amp_24hr_sim_comp_dates,
                                amp_24hr_sim_eqm,
                                amp_24hr_sim_rqm,
                                amp_24hr_sim_gpm,
                                amp_24hr_sim_dqm,
                                amp_24hr_sim_qdm,
                                amp_24hr_sim_raw)

# A comparison data.frame is created
amp_24hr_hist_comp <- data.frame(amp_24hr_hist_comp_dates,
                                 amp_24hr_hist_eqm,
                                 amp_24hr_hist_rqm,
                                 amp_24hr_hist_gpm,
                                 amp_24hr_hist_dqm,
                                 amp_24hr_hist_qdm,
                                 amp_24hr_hist_raw)

# A comparison data.frame is created
amp_24hr_arch_comp <- data.frame(amp_24hr_arch_comp_dates,
                                 amp_24hr_arch_eqm,
                                 amp_24hr_arch_rqm,
                                 amp_24hr_arch_gpm,
                                 amp_24hr_arch_dqm,
                                 amp_24hr_arch_qdm,
                                 amp_24hr_arch_raw)

# data.frame variables are renamed
names(amp_24hr_sim_comp) <- c("dates", "eqm", "rqm", "gpm", "dqm", "qdm", "raw")
names(amp_24hr_hist_comp) <- c("dates", "eqm", "rqm", "gpm", "dqm", "qdm", "raw")
names(amp_24hr_arch_comp) <- c("dates", "eqm", "rqm", "gpm", "dqm", "qdm", "raw")

# ID variables are incorporated
amp_24hr_sim_comp$rcm_id <- rcm.id
amp_24hr_sim_comp$rcp_id <- c("sim")
amp_24hr_sim_comp$sta_id <- sta.id 

# ID variables are incorporated
amp_24hr_hist_comp$rcm_id <- rcm.id
amp_24hr_hist_comp$rcp_id <- c("hist")
amp_24hr_hist_comp$sta_id <- sta.id 

# ID variables are incorporated
amp_24hr_arch_comp$rcm_id <- rcm.id
amp_24hr_arch_comp$rcp_id <- c("arch")
amp_24hr_arch_comp$sta_id <- sta.id 

# Duplicated values are identified per BC method
df_amp_24hr_sim_duplicated <- data.frame(as.vector(sum(duplicated(amp_24hr_sim_comp$eqm))),
                                         as.vector(sum(duplicated(amp_24hr_sim_comp$rqm))),
                                         as.vector(sum(duplicated(amp_24hr_sim_comp$gpm))),
                                         as.vector(sum(duplicated(amp_24hr_sim_comp$dqm))),
                                         as.vector(sum(duplicated(amp_24hr_sim_comp$qdm))))

# Duplicated values are identified per BC method
df_amp_24hr_hist_duplicated <- data.frame(as.vector(sum(duplicated(amp_24hr_hist_comp$eqm))),
                                          as.vector(sum(duplicated(amp_24hr_hist_comp$rqm))),
                                          as.vector(sum(duplicated(amp_24hr_hist_comp$gpm))),
                                          as.vector(sum(duplicated(amp_24hr_hist_comp$dqm))),
                                          as.vector(sum(duplicated(amp_24hr_hist_comp$qdm))))

# Duplicated values are identified per BC method
df_amp_24hr_arch_duplicated <- data.frame(as.vector(sum(duplicated(amp_24hr_arch_comp$eqm))),
                                          as.vector(sum(duplicated(amp_24hr_arch_comp$rqm))),
                                          as.vector(sum(duplicated(amp_24hr_arch_comp$gpm))),
                                          as.vector(sum(duplicated(amp_24hr_arch_comp$dqm))),
                                          as.vector(sum(duplicated(amp_24hr_arch_comp$qdm))))

# data.frame variables are renamed
names(df_amp_24hr_sim_duplicated) <- c("eqm", "rqm", "gpm", "dqm", "qdm")
names(df_amp_24hr_hist_duplicated) <- c("eqm", "rqm", "gpm", "dqm", "qdm")
names(df_amp_24hr_arch_duplicated) <- c("eqm", "rqm", "gpm", "dqm", "qdm")

# ID variables are incorporated
df_amp_24hr_sim_duplicated$rcm_id <- rcm.id
df_amp_24hr_sim_duplicated$rcp_id <- c("sim")
df_amp_24hr_sim_duplicated$sta_id <- sta.id

# ID variables are incorporated
df_amp_24hr_hist_duplicated$rcm_id <- rcm.id
df_amp_24hr_hist_duplicated$rcp_id <- c("hist")
df_amp_24hr_hist_duplicated$sta_id <- sta.id 

# ID variables are incorporated
df_amp_24hr_arch_duplicated$rcm_id <- rcm.id
df_amp_24hr_arch_duplicated$rcp_id <- c("arch")
df_amp_24hr_arch_duplicated$sta_id <- sta.id 

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK:  Bias Correction (BC) performance evaluation
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Performance evaluation metrics for historical daily datasets
# ----------------------------------------------------------------------------------------------------------------------------------

# Descriptive statistics are requested for observed, historical and forecasted datasets
df.base.desc <- c(as.data.frame(round(stat.desc(df_obs_1970_2005[2], norm = FALSE, p=0.95), 3)), # 1970-2005
                  as.data.frame(round(stat.desc(df_obs_1991_2005[2], norm = FALSE, p=0.95), 3)), # 1991-2005
                  as.data.frame(round(stat.desc(df.hindcast.val[2], norm = FALSE, p=0.95), 3)),  # 1991-2005
                  as.data.frame(round(stat.desc(bc_hist_rqm_qmap[2], norm = FALSE, p=0.95), 3)), # 1991-2005
                  as.data.frame(round(stat.desc(bc_hist_gpm_dR[2], norm = FALSE, p=0.95), 3)),   # 1991-2005
                  as.data.frame(round(stat.desc(bc_hist_eqm_dR[2], norm = FALSE, p=0.95), 3)),   # 1991-2005
                  as.data.frame(round(stat.desc(bc_hist_dqm_dR[2], norm = FALSE, p=0.95), 3)),   # 1991-2005
                  as.data.frame(round(stat.desc(bc_hist_qdm_dR[2], norm = FALSE, p=0.95), 3)),   # 1991-2005
                  as.data.frame(round(stat.desc(bc_sim_rqm_qmap[2], norm = FALSE, p=0.95), 3)),  # 1991-2005
                  as.data.frame(round(stat.desc(bc_sim_gpm_dR[2], norm = FALSE, p=0.95), 3)),    # 1991-2005
                  as.data.frame(round(stat.desc(bc_sim_eqm_dR[2], norm = FALSE, p=0.95), 3)),    # 1991-2005
                  as.data.frame(round(stat.desc(bc_sim_dqm_dR[2], norm = FALSE, p=0.95), 3)),    # 1991-2005
                  as.data.frame(round(stat.desc(bc_sim_qdm_dR[2], norm = FALSE, p=0.95), 3)))    # 1991-2005

# An external list of data.frames is created
df.base.desc <- (sapply(df.base.desc, unlist))

# A data.frame is created
df.base.desc <- data.frame(df.base.desc)

# data.frame variables are renamed
names(df.base.desc) <- c("obs_1970_2005",        # 13149 obs 1970-2005
                         "obs_1991_2005",        # 5479 obs 1991-2005
                         "hindcast_1991_2005",   # 5479 obs 1991-2005
                         "bc_hist_rqm",          # 5479 obs 1991-2005
                         "bc_hist_gpm",          # 5479 obs 1991-2005
                         "bc_hist_eqm",          # 5479 obs 1991-2005
                         "bc_hist_dqm",          # 5479 obs 1991-2005
                         "bc_hist_qdm",          # 5479 obs 1991-2005
                         "bc_sim_rqm",           # 34333 obs 2006-2099
                         "bc_sim_gpm",           # 34333 obs 2006-2099
                         "bc_sim_eqm",           # 34333 obs 2006-2099
                         "bc_sim_dqm",           # 34333 obs 2006-2099
                         "bc_sim_qdm")           # 34333 obs 2006-2099

# data.frame row names are renamed
rownames(df.base.desc) <- c("nbr.val", "NBR.ZERO", "nbr.na", "min",
                            "MAX", "range", "SUM", "MEDIAN", "MEAN",
                            "SE.mean", "CI.mean.0.95", "var",
                            "STD.DEV", "coef.var")

# Special attention to parameters: 
# NBR.ZERO: number of days with null precipitation compared to hindcast_1991_2005
# MAX: max precipitation for obs hindcast_1991_2005
# SUM: sum of precipitation for obs hindcast_1991_2005
# MEDIAN : median precipitation for obs hindcast_1991_2005
# MEAN : mean precipitation for obs hindcast_1991_2005
# STD.DEV: Standard deviation for obs hindcast_1991_2005

# ID variables are incorporated
df.base.desc$rcm_id <- rcm.id
df.base.desc$rcp_id <- c("hist")
df.base.desc$sta_id <- sta.id

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Performance evaluation metrics for historical daily AMP24 datasets
# ----------------------------------------------------------------------------------------------------------------------------------

# {lubridate} functions are applied to create new columns containing year:
df_obs_1991_2005$YEAR <- lubridate::year(df_obs_1991_2005$date)     # 1991-2005
bc_hist_eqm_dR$YEAR <- lubridate::year(bc_hist_eqm_dR$date)         # 1991-2005
bc_hist_rqm_qmap$YEAR <- lubridate::year(bc_hist_rqm_qmap$date)     # 1991-2005
bc_hist_gpm_dR$YEAR <- lubridate::year(bc_hist_gpm_dR$date)         # 1991-2005
bc_hist_dqm_dR$YEAR <- lubridate::year(bc_hist_dqm_dR$date)         # 1991-2005
bc_hist_qdm_dR$YEAR <- lubridate::year(bc_hist_qdm_dR$date)         # 1991-2005
df.hindcast.val$YEAR <- lubridate::year(df.hindcast.val$date)       # 1991-2005

# {lubridate} functions are applied to create new columns containing month:
df_obs_1991_2005$MONTH <- lubridate::month(df_obs_1991_2005$date)   # 1991-2005
bc_hist_eqm_dR$MONTH <- lubridate::month(bc_hist_eqm_dR$date)       # 1991-2005
bc_hist_rqm_qmap$MONTH <- lubridate::month(bc_hist_rqm_qmap$date)   # 1991-2005
bc_hist_gpm_dR$MONTH <- lubridate::month(bc_hist_gpm_dR$date)       # 1991-2005
bc_hist_dqm_dR$MONTH <- lubridate::month(bc_hist_dqm_dR$date)       # 1991-2005
bc_hist_qdm_dR$MONTH <- lubridate::month(bc_hist_qdm_dR$date)       # 1991-2005
df.hindcast.val$MONTH <- lubridate::month(df.hindcast.val$date)     # 1991-2005

# Observed data are aggregated by year to determine daily AMP24s monthly
amp_24hrAMP_hist_obs <- df_obs_1991_2005 %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_eqm <- bc_hist_eqm_dR %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_rqm <- bc_hist_rqm_qmap %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_gpm <- bc_hist_gpm_dR %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_dqm <- bc_hist_dqm_dR %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_qdm <- bc_hist_qdm_dR %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))
amp_24hrAMP_hist_raw <- df.hindcast.val %>%
  group_by(YEAR, MONTH) %>%
  slice(which.max(prec))

# Variable names are added to data.frame
amp_24hrAMP_hist_comp_dates <- amp_24hrAMP_hist_eqm$dates# rownames(amp_24hrAMP_hist_eqm)

# A summary monthly data.frame is created
amp_24hrAMP_hist_comp <- data.frame(amp_24hrAMP_hist_comp_dates,  # 180 obs 1991-2005
                                    amp_24hrAMP_hist_obs$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_eqm$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_rqm$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_gpm$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_dqm$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_qdm$prec,    # 180 obs 1991-2005
                                    amp_24hrAMP_hist_raw$prec)    # 180 obs 1991-2005

# data.frame variables are renamed
names(amp_24hrAMP_hist_comp) <- c("dates", "obs", "eqm", "rqm", "gpm", "dqm", "qdm","raw")

# ID variables are incorporated
amp_24hrAMP_hist_comp$rcm_id <- rcm.id
amp_24hrAMP_hist_comp$rcp_id <- c("hist")
amp_24hrAMP_hist_comp$sta_id <- sta.id 

# Duplicated monthly values are identified per BC method
df_amp_24hrAMP_hist_duplicated <- data.frame(as.vector(sum(duplicated(amp_24hrAMP_hist_comp$eqm))),  # 180 obs 1991-2005
                                             as.vector(sum(duplicated(amp_24hrAMP_hist_comp$rqm))),  # 180 obs 1991-2005
                                             as.vector(sum(duplicated(amp_24hrAMP_hist_comp$gpm))),  # 180 obs 1991-2005
                                             as.vector(sum(duplicated(amp_24hrAMP_hist_comp$dqm))),  # 180 obs 1991-2005
                                             as.vector(sum(duplicated(amp_24hrAMP_hist_comp$qdm))))  # 180 obs 1991-2005

# data.frame variables are renamed
names(df_amp_24hrAMP_hist_duplicated) <- c("eqm", "rqm", "gpm", "dqm", "qdm")

# ID variables are incorporated
df_amp_24hrAMP_hist_duplicated$rcm_id <- rcm.id
df_amp_24hrAMP_hist_duplicated$rcp_id <- c("hist")
df_amp_24hrAMP_hist_duplicated$sta_id <- sta.id 

# A nRMSE historical monthly vector is created
v.rmse.24hrAMP <- c(hydroGOF::nrmse(amp_24hrAMP_hist_comp$rqm, amp_24hrAMP_hist_comp$obs, norm="maxmin"), # 180 obs
                    hydroGOF::nrmse(amp_24hrAMP_hist_comp$gpm, amp_24hrAMP_hist_comp$obs, norm="maxmin"), # 180 obs
                    hydroGOF::nrmse(amp_24hrAMP_hist_comp$eqm, amp_24hrAMP_hist_comp$obs, norm="maxmin"), # 180 obs
                    hydroGOF::nrmse(amp_24hrAMP_hist_comp$dqm, amp_24hrAMP_hist_comp$obs, norm="maxmin"), # 180 obs
                    hydroGOF::nrmse(amp_24hrAMP_hist_comp$qdm, amp_24hrAMP_hist_comp$obs, norm="maxmin"), # 180 obs
                    hydroGOF::nrmse(amp_24hrAMP_hist_comp$raw, amp_24hrAMP_hist_comp$obs, norm="maxmin")) # 180 obs

# A MBE Mean Bias Error historical monthly vector is created
v.bme.24hrAMP <- c(tdr::tdStats(amp_24hrAMP_hist_comp$rqm, amp_24hrAMP_hist_comp$obs, functions = "mbe"),
                   tdr::tdStats(amp_24hrAMP_hist_comp$gpm, amp_24hrAMP_hist_comp$obs, functions = "mbe"),
                   tdr::tdStats(amp_24hrAMP_hist_comp$eqm, amp_24hrAMP_hist_comp$obs, functions = "mbe"),
                   tdr::tdStats(amp_24hrAMP_hist_comp$dqm, amp_24hrAMP_hist_comp$obs, functions = "mbe"),
                   tdr::tdStats(amp_24hrAMP_hist_comp$qdm, amp_24hrAMP_hist_comp$obs, functions = "mbe"),
                   tdr::tdStats(amp_24hrAMP_hist_comp$raw, amp_24hrAMP_hist_comp$obs, functions = "mbe"))

# A MDA Modified index of agreement vector is created for historical monthly
v.mda.24hrAMP <- c(hydroGOF::md(amp_24hrAMP_hist_comp$rqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::md(amp_24hrAMP_hist_comp$gpm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::md(amp_24hrAMP_hist_comp$eqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::md(amp_24hrAMP_hist_comp$dqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::md(amp_24hrAMP_hist_comp$qdm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::md(amp_24hrAMP_hist_comp$raw, amp_24hrAMP_hist_comp$obs))

# A Percent Bias vector is created for historical monthly
v.pbias.24hrAMP <- c(hydroGOF::pbias(amp_24hrAMP_hist_comp$rqm, amp_24hrAMP_hist_comp$obs),
                     hydroGOF::pbias(amp_24hrAMP_hist_comp$gpm, amp_24hrAMP_hist_comp$obs),
                     hydroGOF::pbias(amp_24hrAMP_hist_comp$eqm, amp_24hrAMP_hist_comp$obs),
                     hydroGOF::pbias(amp_24hrAMP_hist_comp$dqm, amp_24hrAMP_hist_comp$obs),
                     hydroGOF::pbias(amp_24hrAMP_hist_comp$qdm, amp_24hrAMP_hist_comp$obs),
                     hydroGOF::pbias(amp_24hrAMP_hist_comp$raw, amp_24hrAMP_hist_comp$obs))

# A Mean Absolute Error vector is created for historical monthly
v.mae.24hrAMP <- c(hydroGOF::mae(amp_24hrAMP_hist_comp$rqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::mae(amp_24hrAMP_hist_comp$gpm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::mae(amp_24hrAMP_hist_comp$eqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::mae(amp_24hrAMP_hist_comp$dqm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::mae(amp_24hrAMP_hist_comp$qdm, amp_24hrAMP_hist_comp$obs),
                   hydroGOF::mae(amp_24hrAMP_hist_comp$raw, amp_24hrAMP_hist_comp$obs))

# A historical monthly performance data.frame is created
df.perf.desc.24hrAMP <- data.frame(v.rmse.24hrAMP, v.bme.24hrAMP, v.mda.24hrAMP, v.pbias.24hrAMP, v.mae.24hrAMP) # 180 obs

# data.frame variables are renamed
names(df.perf.desc.24hrAMP) <- c("nRMSE","MBE", "MDA", "PBIAS", "MAE")

# data.frame row names are renamed
rownames(df.perf.desc.24hrAMP) <- c("hist_rqm_AMP", "hist_gpm_AMP", "hist_eqm_AMP",
                                    "hist_dqm_AMP", "hist_qdm_AMP", "hist_raw_AMP")

# ID variables are incorporated
df.perf.desc.24hrAMP$rcm_id <- rcm.id
df.perf.desc.24hrAMP$rcp_id <- c("hist")
df.perf.desc.24hrAMP$sta_id <- sta.id 

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Performance evaluation metrics for CORDEX member during baseline 1971-2000
# ----------------------------------------------------------------------------------------------------------------------------------

# A nRMSE historical vector is created
v.rmse.baseline <- hydroGOF::nrmse(df.hindcast.baseline$prec, df_obs_1971_2000$prec, norm="maxmin") # 10958 obs

# A MBE Mean Bias Error historical vector is created
v.mbe.baseline <- tdr::tdStats(df.hindcast.baseline$prec, df_obs_1971_2000$prec, functions = "mbe") # 10958 obs

# A MDA Modified index of agreement vector is created
v.mda.baseline <- hydroGOF::md(df.hindcast.baseline$prec, df_obs_1971_2000$prec) # 10958 obs

# A Percent Bias vector is created
v.pbias.baseline <- hydroGOF::pbias(df.hindcast.baseline$prec, df_obs_1971_2000$prec) # 10958 obs

# A Mean Absolute Error vector is created
v.mae.baseline <- hydroGOF::mae(df.hindcast.baseline$prec, df_obs_1971_2000$prec) # 10958 obs

# A historical performance data.frame is created
df.perf.desc.daily.baseline <- data.frame(v.rmse.baseline, v.mbe.baseline, v.mda.baseline, v.pbias.baseline, v.mae.baseline) # 10958 obs

# data.frame variables are renamed
names(df.perf.desc.daily.baseline) <- c("nRMSE", "MBE", "MDA", "PBIAS", "MAE")

# data.frame row names are renamed
rownames(df.perf.desc.daily.baseline) <- rcm.id

# ID variables are incorporated
df.perf.desc.daily.baseline$rcm_id <- rcm.id
df.perf.desc.daily.baseline$rcp_id <- c("hist")
df.perf.desc.daily.baseline$sta_id <- sta.id 

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Performance evaluation metrics for historical daily datasets for each BC method
# ----------------------------------------------------------------------------------------------------------------------------------

# A nRMSE historical vector is created
v.rmse <- c(hydroGOF::nrmse(bc_hist_rqm_qmap$prec, df_obs_1991_2005$prec, norm="maxmin"), # 5479 obs
            hydroGOF::nrmse(bc_hist_gpm_dR$prec, df_obs_1991_2005$prec, norm="maxmin"),   # 5479 obs
            hydroGOF::nrmse(bc_hist_eqm_dR$prec, df_obs_1991_2005$prec, norm="maxmin"),   # 5479 obs
            hydroGOF::nrmse(bc_hist_dqm_dR$prec, df_obs_1991_2005$prec, norm="maxmin"),   # 5479 obs
            hydroGOF::nrmse(bc_hist_qdm_dR$prec, df_obs_1991_2005$prec, norm="maxmin"),   # 5479 obs
            hydroGOF::nrmse(df.hindcast.val$prec, df_obs_1991_2005$prec, norm="maxmin"))  # 5479 obs

# A MBE Mean Bias Error historical vector is created
v.mbe <- c(tdr::tdStats(bc_hist_rqm_qmap$prec, df_obs_1991_2005$prec, functions = "mbe"),
           tdr::tdStats(bc_hist_gpm_dR$prec, df_obs_1991_2005$prec, functions = "mbe"),
           tdr::tdStats(bc_hist_eqm_dR$prec, df_obs_1991_2005$prec, functions = "mbe"),
           tdr::tdStats(bc_hist_dqm_dR$prec, df_obs_1991_2005$prec, functions = "mbe"),
           tdr::tdStats(bc_hist_qdm_dR$prec, df_obs_1991_2005$prec, functions = "mbe"),
           tdr::tdStats(df.hindcast.val$prec, df_obs_1991_2005$prec, functions = "mbe"))

# A MDA Modified index of agreement vector is created
v.mda <- c(hydroGOF::md(bc_hist_rqm_qmap$prec, df_obs_1991_2005$prec),
           hydroGOF::md(bc_hist_gpm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::md(bc_hist_eqm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::md(bc_hist_dqm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::md(bc_hist_qdm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::md(df.hindcast.val$prec, df_obs_1991_2005$prec))

# A Percent Bias vector is created
v.pbias <- c(hydroGOF::pbias(bc_hist_rqm_qmap$prec, df_obs_1991_2005$prec),
             hydroGOF::pbias(bc_hist_gpm_dR$prec, df_obs_1991_2005$prec),
             hydroGOF::pbias(bc_hist_eqm_dR$prec, df_obs_1991_2005$prec),
             hydroGOF::pbias(bc_hist_dqm_dR$prec, df_obs_1991_2005$prec),
             hydroGOF::pbias(bc_hist_qdm_dR$prec, df_obs_1991_2005$prec),
             hydroGOF::pbias(df.hindcast.val$prec, df_obs_1991_2005$prec))

# A Mean Absolute Error vector is created
v.mae <- c(hydroGOF::mae(bc_hist_rqm_qmap$prec, df_obs_1991_2005$prec),
           hydroGOF::mae(bc_hist_gpm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::mae(bc_hist_eqm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::mae(bc_hist_dqm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::mae(bc_hist_qdm_dR$prec, df_obs_1991_2005$prec),
           hydroGOF::mae(df.hindcast.val$prec, df_obs_1991_2005$prec))

# A historical performance data.frame is created
df.perf.desc.daily <- data.frame(v.rmse, v.mbe, v.mda, v.pbias, v.mae) # 5479 obs

# data.frame variables are renamed
names(df.perf.desc.daily) <- c("nRMSE", "MBE", "MDA", "PBIAS", "MAE")

# data.frame row names are renamed
rownames(df.perf.desc.daily) <- c("hist_rqm_daily", "hist_gpm_daily", "hist_eqm_daily",
                                  "hist_dqm_daily", "hist_qdm_daily", "hist_raw_daily")

# ID variables are incorporated
df.perf.desc.daily$rcm_id <- rcm.id
df.perf.desc.daily$rcp_id <- c("hist")
df.perf.desc.daily$sta_id <- sta.id 

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: climdex.prcptot {climdex.pcic} ETCCDI during baseline 1971-2000
# ----------------------------------------------------------------------------------------------------------------------------------

# Dates are parsed into PCICt class
prec_dates_val_baseline <- as.PCICt(as.character(df_obs_1971_2000$date), cal="gregorian") # 10958 obs

# climdexInput object from vectors of data are created
ci_val_raw_baseline_baseline <- climdexInput.raw(prec = df.hindcast.baseline$prec,     # 10958 obs
                                                 prec.dates = prec_dates_val_baseline, # 10958 obs
                                                 base.range=c(1971, 2000))

# climdexInput object from vectors of data are created
ci_val_obs_baseline <- climdexInput.raw(prec = df_obs_1971_2000$prec,         # 10958 obs
                                        prec.dates = prec_dates_val_baseline, # 10958 obs
                                        base.range=c(1971, 2000))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# climdex.prcptot {climdex.pcic} Total Daily Prec. prcptot[mm/year]
# ETCCDI Climate Change Indices. Definitions of the 27 core indices
# http://etccdi.pacificclimate.org/list_27_indices.shtml
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total prec.wet days [mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prcptot_val_raw_baseline <- ceiling(mean(climdex.prcptot(ci_val_raw_baseline_baseline), na.rm = TRUE))
prcptot_val_obs_baseline <- ceiling(mean(climdex.prcptot(ci_val_obs_baseline), na.rm = TRUE))
v.prcptot_baseline <- c(prcptot_val_raw_baseline)

prcptot_val_raw_baseline_full <- climdex.prcptot(ci_val_raw_baseline_baseline)
prcptot_val_obs_baseline_full <- climdex.prcptot(ci_val_obs_baseline)
v.prcptot_baseline_full <- c(prcptot_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive wet days cwd RR ≥ 1mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd_val_raw_baseline <- ceiling(mean(climdex.cwd(ci_val_raw_baseline_baseline), na.rm = TRUE))
cwd_val_obs_baseline <- ceiling(mean(climdex.cwd(ci_val_obs_baseline), na.rm = TRUE))
v.cwd_baseline <- c(cwd_val_raw_baseline, cwd_val_obs_baseline)

cwd_val_raw_baseline_full <- climdex.cwd(ci_val_raw_baseline_baseline)
cwd_val_obs_baseline_full <- climdex.cwd(ci_val_obs_baseline)
v.cwd_baseline_full <- c(cwd_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive dry days cdd [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdd_val_raw_baseline <- ceiling(mean(climdex.cdd(ci_val_raw_baseline_baseline), na.rm = TRUE))
cdd_val_obs_baseline <- ceiling(mean(climdex.cdd(ci_val_obs_baseline), na.rm = TRUE))
v.cdd_baseline <- c(cdd_val_raw_baseline, cdd_val_obs_baseline)

cdd_val_raw_baseline_full <- climdex.cdd(ci_val_raw_baseline_baseline)
cdd_val_obs_baseline_full <- climdex.cdd(ci_val_obs_baseline)
v.cdd_baseline_full <- c(cdd_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 10mm r10mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r10mm_val_raw_baseline <- ceiling(mean(climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 10), na.rm = TRUE))
r10mm_val_obs_baseline <- ceiling(mean(climdex.rnnmm(ci_val_obs_baseline, threshold = 10), na.rm = TRUE))
v.r10mm_baseline <- c(r10mm_val_raw_baseline, r10mm_val_obs_baseline)

r10mm_val_raw_baseline_full <- climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 10)
r10mm_val_obs_baseline_full <- climdex.rnnmm(ci_val_obs_baseline, threshold = 10)
v.r10mm_baseline_full <- c(r10mm_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 20mm r20mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r20mm_val_raw_baseline <- ceiling(mean(climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 20), na.rm = TRUE))
r20mm_val_obs_baseline <- ceiling(mean(climdex.rnnmm(ci_val_obs_baseline, threshold = 20), na.rm = TRUE))
v.r20mm_baseline <- c(r20mm_val_raw_baseline, r20mm_val_obs_baseline)

r20mm_val_raw_baseline_full <- climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 20)
r20mm_val_obs_baseline_full <- climdex.rnnmm(ci_val_obs_baseline, threshold = 20)
v.r20mm_baseline_full <- c(r20mm_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 30mm r30mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r30mm_val_raw_baseline <- ceiling(mean(climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 30), na.rm = TRUE))
r30mm_val_obs_baseline <- ceiling(mean(climdex.rnnmm(ci_val_obs_baseline, threshold = 30), na.rm = TRUE))
v.r30mm_baseline <- c(r30mm_val_raw_baseline, r30mm_val_obs_baseline)

r30mm_val_raw_baseline_full <- climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 30)
r30mm_val_obs_baseline_full <- climdex.rnnmm(ci_val_obs_baseline, threshold = 30)
v.r30mm_baseline_full <- c(r30mm_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 50mm r50mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r50mm_val_raw_baseline <- ceiling(mean(climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 50), na.rm = TRUE))
r50mm_val_obs_baseline <- ceiling(mean(climdex.rnnmm(ci_val_obs_baseline, threshold = 50), na.rm = TRUE))
v.r50mm_baseline <- c(r50mm_val_raw_baseline, r50mm_val_obs_baseline)

r50mm_val_raw_baseline_full <- climdex.rnnmm(ci_val_raw_baseline_baseline, threshold = 50)
r50mm_val_obs_baseline_full <- climdex.rnnmm(ci_val_obs_baseline, threshold = 50)
v.r50mm_baseline_full <- c(r50mm_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 95p r95p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r95p_val_raw_baseline <- ceiling(mean(climdex.r95ptot(ci_val_raw_baseline_baseline), na.rm = TRUE))
r95p_val_obs_baseline <- ceiling(mean(climdex.r95ptot(ci_val_obs_baseline), na.rm = TRUE))
v.r95p_baseline <- c(r95p_val_raw_baseline, r95p_val_obs_baseline)

r95p_val_raw_baseline_full <- climdex.r95ptot(ci_val_raw_baseline_baseline)
r95p_val_obs_baseline_full <- climdex.r95ptot(ci_val_obs_baseline)
v.r95p_baseline_full <- c(r95p_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 99p r99p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r99p_val_raw_baseline <- ceiling(mean(climdex.r99ptot(ci_val_raw_baseline_baseline), na.rm = TRUE))
r99p_val_obs_baseline <- ceiling(mean(climdex.r99ptot(ci_val_obs_baseline), na.rm = TRUE))
v.r99p_baseline <- c(r99p_val_raw_baseline, r99p_val_obs_baseline)

r99p_val_raw_baseline_full <- climdex.r99ptot(ci_val_raw_baseline_baseline)
r99p_val_obs_baseline_full <- climdex.r99ptot(ci_val_obs_baseline)
v.r99p_baseline_full <- c(r99p_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum 1-day precipitation. rx1day[mm/1day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx1day_val_raw_baseline <- ceiling(max(climdex.rx1day(ci_val_raw_baseline_baseline, freq = c("monthly")), na.rm = TRUE))
rx1day_val_obs_baseline <- ceiling(max(climdex.rx1day(ci_val_obs_baseline, freq = c("monthly")), na.rm = TRUE))
v.rx1day_baseline <- c(rx1day_val_raw_baseline, rx1day_val_obs_baseline)

rx1day_val_raw_baseline_full <- climdex.rx1day(ci_val_raw_baseline_baseline, freq = c("monthly"))
rx1day_val_obs_baseline_full <- climdex.rx1day(ci_val_obs_baseline, freq = c("monthly"))
v.rx1day_baseline_full <- c(rx1day_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum consecutive 5-day precipitation. Rx5day [mm/5day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx5day_val_raw_baseline <- ceiling(max(climdex.rx5day(ci_val_raw_baseline_baseline, freq = c("monthly")), na.rm = TRUE))
rx5day_val_obs_baseline <- ceiling(max(climdex.rx5day(ci_val_obs_baseline, freq = c("monthly")), na.rm = TRUE))
v.rx5day_baseline <- c(rx5day_val_raw_baseline, rx5day_val_obs_baseline)

rx5day_val_raw_baseline_full <- climdex.rx5day(ci_val_raw_baseline_baseline, freq = c("monthly"))
rx5day_val_obs_baseline_full <- climdex.rx5day(ci_val_obs_baseline, freq = c("monthly"))
v.rx5day_baseline_full <- c(rx5day_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Simple precipitation intensity index sdii [mm/day] of wet days ONLY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdii_val_raw_baseline <- ceiling(mean(climdex.sdii(ci_val_raw_baseline_baseline), na.rm = TRUE))
sdii_val_obs_baseline <- ceiling(mean(climdex.sdii(ci_val_obs_baseline), na.rm = TRUE))
v.sdii_baseline <- c(sdii_val_raw_baseline, sdii_val_obs_baseline)

sdii_val_raw_baseline_full <- climdex.sdii(ci_val_raw_baseline_baseline)
sdii_val_obs_baseline_full <- climdex.sdii(ci_val_obs_baseline)
v.sdii_baseline_full <- c(sdii_val_raw_baseline_full)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# A climdex.baseline performance data.frame is created
df.perf.climdex.baseline <- data.frame(v.prcptot_baseline, v.cwd_baseline, v.cdd_baseline,
                                       v.r10mm_baseline, v.r20mm_baseline, v.r30mm_baseline,
                                       v.r50mm_baseline, v.r95p_baseline, v.r99p_baseline,
                                       v.rx1day_baseline, v.rx5day_baseline, v.sdii_baseline)

# data.frame variables are renamed
names(df.perf.climdex.baseline) <- c("prcptot", "cwd", "cdd",
                                     "r10mm", "r20mm", "r30mm",
                                     "r50mm", "r95p", "r99p",
                                     "rx1day", "rx5day", "sdii")

# data.frame row names are renamed
rownames(df.perf.climdex.baseline) <- c("raw", "obs")

# ID variables are incorporated
df.perf.climdex.baseline$rcm_id <- rcm.id
df.perf.climdex.baseline$rcp_id <- c("hist")
df.perf.climdex.baseline$sta_id <- sta.id 

# A climdex.baseline performance data.frame is created
df.perf.climdex.baseline_full <- c(v.prcptot_baseline_full, v.cwd_baseline_full, v.cdd_baseline_full,
                                   v.r10mm_baseline_full, v.r20mm_baseline_full, v.r30mm_baseline_full,
                                   v.r50mm_baseline_full, v.r95p_baseline_full, v.r99p_baseline_full,
                                   v.rx1day_baseline_full, v.rx5day_baseline_full, v.sdii_baseline_full)

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: climdex.prcptot {climdex.pcic} Total Daily Precipitation implementation of the ETCCDI climate change indices
# ----------------------------------------------------------------------------------------------------------------------------------

# Dates are parsed into PCICt class
prec_dates_val <- as.PCICt(as.character(df_obs_1991_2005$date), cal="gregorian")

# climdexInput object from vectors of data are created
ci_val_raw <- climdexInput.raw(prec = df.hindcast.val$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_obs <- climdexInput.raw(prec = df_obs_1991_2005$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_rqm <- climdexInput.raw(prec = bc_hist_rqm_qmap$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_gpm <- climdexInput.raw(prec = bc_hist_gpm_dR$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_eqm <- climdexInput.raw(prec = bc_hist_eqm_dR$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_dqm <- climdexInput.raw(prec = bc_hist_dqm_dR$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_val_qdm <- climdexInput.raw(prec = bc_hist_qdm_dR$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# climdex.prcptot {climdex.pcic} Total Daily Prec. prcptot[mm/year]
# ETCCDI Climate Change Indices. Definitions of the 27 core indices
# http://etccdi.pacificclimate.org/list_27_indices.shtml
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total prec.wet days [mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prcptot_val_raw <- ceiling(mean(climdex.prcptot(ci_val_raw), na.rm = TRUE))
prcptot_val_obs <- ceiling(mean(climdex.prcptot(ci_val_obs), na.rm = TRUE))
prcptot_val_rqm <- ceiling(mean(climdex.prcptot(ci_val_rqm), na.rm = TRUE))
prcptot_val_gpm <- ceiling(mean(climdex.prcptot(ci_val_gpm), na.rm = TRUE))
prcptot_val_eqm <- ceiling(mean(climdex.prcptot(ci_val_eqm), na.rm = TRUE))
prcptot_val_dqm <- ceiling(mean(climdex.prcptot(ci_val_dqm), na.rm = TRUE))
prcptot_val_qdm <- ceiling(mean(climdex.prcptot(ci_val_qdm), na.rm = TRUE))
v.prcptot <- c(prcptot_val_raw, prcptot_val_obs, prcptot_val_rqm, prcptot_val_gpm, 
               prcptot_val_eqm, prcptot_val_dqm, prcptot_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive wet days cwd RR ≥ 1mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd_val_raw <- ceiling(mean(climdex.cwd(ci_val_raw), na.rm = TRUE))
cwd_val_obs <- ceiling(mean(climdex.cwd(ci_val_obs), na.rm = TRUE))
cwd_val_rqm <- ceiling(mean(climdex.cwd(ci_val_rqm), na.rm = TRUE))
cwd_val_gpm <- ceiling(mean(climdex.cwd(ci_val_gpm), na.rm = TRUE))
cwd_val_eqm <- ceiling(mean(climdex.cwd(ci_val_eqm), na.rm = TRUE))
cwd_val_dqm <- ceiling(mean(climdex.cwd(ci_val_dqm), na.rm = TRUE))
cwd_val_qdm <- ceiling(mean(climdex.cwd(ci_val_qdm), na.rm = TRUE))
v.cwd <- c(cwd_val_raw, cwd_val_obs, cwd_val_rqm, cwd_val_gpm, 
           cwd_val_eqm, cwd_val_dqm, cwd_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive dry days cdd [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdd_val_raw <- ceiling(mean(climdex.cdd(ci_val_raw), na.rm = TRUE))
cdd_val_obs <- ceiling(mean(climdex.cdd(ci_val_obs), na.rm = TRUE))
cdd_val_rqm <- ceiling(mean(climdex.cdd(ci_val_rqm), na.rm = TRUE))
cdd_val_gpm <- ceiling(mean(climdex.cdd(ci_val_gpm), na.rm = TRUE))
cdd_val_eqm <- ceiling(mean(climdex.cdd(ci_val_eqm), na.rm = TRUE))
cdd_val_dqm <- ceiling(mean(climdex.cdd(ci_val_dqm), na.rm = TRUE))
cdd_val_qdm <- ceiling(mean(climdex.cdd(ci_val_qdm), na.rm = TRUE))
v.cdd <- c(cdd_val_raw, cdd_val_obs, cdd_val_rqm, cdd_val_gpm, 
           cdd_val_eqm, cdd_val_dqm, cdd_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 10mm r10mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r10mm_val_raw <- ceiling(mean(climdex.rnnmm(ci_val_raw, threshold = 10), na.rm = TRUE))
r10mm_val_obs <- ceiling(mean(climdex.rnnmm(ci_val_obs, threshold = 10), na.rm = TRUE))
r10mm_val_rqm <- ceiling(mean(climdex.rnnmm(ci_val_rqm, threshold = 10), na.rm = TRUE))
r10mm_val_gpm <- ceiling(mean(climdex.rnnmm(ci_val_gpm, threshold = 10), na.rm = TRUE))
r10mm_val_eqm <- ceiling(mean(climdex.rnnmm(ci_val_eqm, threshold = 10), na.rm = TRUE))
r10mm_val_dqm <- ceiling(mean(climdex.rnnmm(ci_val_dqm, threshold = 10), na.rm = TRUE))
r10mm_val_qdm <- ceiling(mean(climdex.rnnmm(ci_val_qdm, threshold = 10), na.rm = TRUE))
v.r10mm <- c(r10mm_val_raw, r10mm_val_obs, r10mm_val_rqm, r10mm_val_gpm, 
             r10mm_val_eqm, r10mm_val_dqm, r10mm_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 20mm r20mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r20mm_val_raw <- ceiling(mean(climdex.rnnmm(ci_val_raw, threshold = 20), na.rm = TRUE))
r20mm_val_obs <- ceiling(mean(climdex.rnnmm(ci_val_obs, threshold = 20), na.rm = TRUE))
r20mm_val_rqm <- ceiling(mean(climdex.rnnmm(ci_val_rqm, threshold = 20), na.rm = TRUE))
r20mm_val_gpm <- ceiling(mean(climdex.rnnmm(ci_val_gpm, threshold = 20), na.rm = TRUE))
r20mm_val_eqm <- ceiling(mean(climdex.rnnmm(ci_val_eqm, threshold = 20), na.rm = TRUE))
r20mm_val_dqm <- ceiling(mean(climdex.rnnmm(ci_val_dqm, threshold = 20), na.rm = TRUE))
r20mm_val_qdm <- ceiling(mean(climdex.rnnmm(ci_val_qdm, threshold = 20), na.rm = TRUE))
v.r20mm <- c(r20mm_val_raw, r20mm_val_obs, r20mm_val_rqm, r20mm_val_gpm, 
             r20mm_val_eqm, r20mm_val_dqm, r20mm_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 30mm r30mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r30mm_val_raw <- ceiling(mean(climdex.rnnmm(ci_val_raw, threshold = 30), na.rm = TRUE))
r30mm_val_obs <- ceiling(mean(climdex.rnnmm(ci_val_obs, threshold = 30), na.rm = TRUE))
r30mm_val_rqm <- ceiling(mean(climdex.rnnmm(ci_val_rqm, threshold = 30), na.rm = TRUE))
r30mm_val_gpm <- ceiling(mean(climdex.rnnmm(ci_val_gpm, threshold = 30), na.rm = TRUE))
r30mm_val_eqm <- ceiling(mean(climdex.rnnmm(ci_val_eqm, threshold = 30), na.rm = TRUE))
r30mm_val_dqm <- ceiling(mean(climdex.rnnmm(ci_val_dqm, threshold = 30), na.rm = TRUE))
r30mm_val_qdm <- ceiling(mean(climdex.rnnmm(ci_val_qdm, threshold = 30), na.rm = TRUE))
v.r30mm <- c(r30mm_val_raw, r30mm_val_obs, r30mm_val_rqm, r30mm_val_gpm, 
             r30mm_val_eqm, r30mm_val_dqm, r30mm_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 50mm r50mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r50mm_val_raw <- ceiling(mean(climdex.rnnmm(ci_val_raw, threshold = 50), na.rm = TRUE))
r50mm_val_obs <- ceiling(mean(climdex.rnnmm(ci_val_obs, threshold = 50), na.rm = TRUE))
r50mm_val_rqm <- ceiling(mean(climdex.rnnmm(ci_val_rqm, threshold = 50), na.rm = TRUE))
r50mm_val_gpm <- ceiling(mean(climdex.rnnmm(ci_val_gpm, threshold = 50), na.rm = TRUE))
r50mm_val_eqm <- ceiling(mean(climdex.rnnmm(ci_val_eqm, threshold = 50), na.rm = TRUE))
r50mm_val_dqm <- ceiling(mean(climdex.rnnmm(ci_val_dqm, threshold = 50), na.rm = TRUE))
r50mm_val_qdm <- ceiling(mean(climdex.rnnmm(ci_val_qdm, threshold = 50), na.rm = TRUE))
v.r50mm <- c(r50mm_val_raw, r50mm_val_obs, r50mm_val_rqm, r50mm_val_gpm, 
             r50mm_val_eqm, r50mm_val_dqm, r50mm_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 95p r95p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r95p_val_raw <- ceiling(mean(climdex.r95ptot(ci_val_raw), na.rm = TRUE))
r95p_val_obs <- ceiling(mean(climdex.r95ptot(ci_val_obs), na.rm = TRUE))
r95p_val_rqm <- ceiling(mean(climdex.r95ptot(ci_val_rqm), na.rm = TRUE))
r95p_val_gpm <- ceiling(mean(climdex.r95ptot(ci_val_gpm), na.rm = TRUE))
r95p_val_eqm <- ceiling(mean(climdex.r95ptot(ci_val_eqm), na.rm = TRUE))
r95p_val_dqm <- ceiling(mean(climdex.r95ptot(ci_val_dqm), na.rm = TRUE))
r95p_val_qdm <- ceiling(mean(climdex.r95ptot(ci_val_qdm), na.rm = TRUE))
v.r95p <- c(r95p_val_raw, r95p_val_obs, r95p_val_rqm, r95p_val_gpm, 
            r95p_val_eqm, r95p_val_dqm, r95p_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 99p r99p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r99p_val_raw <- ceiling(mean(climdex.r99ptot(ci_val_raw), na.rm = TRUE))
r99p_val_obs <- ceiling(mean(climdex.r99ptot(ci_val_obs), na.rm = TRUE))
r99p_val_rqm <- ceiling(mean(climdex.r99ptot(ci_val_rqm), na.rm = TRUE))
r99p_val_gpm <- ceiling(mean(climdex.r99ptot(ci_val_gpm), na.rm = TRUE))
r99p_val_eqm <- ceiling(mean(climdex.r99ptot(ci_val_eqm), na.rm = TRUE))
r99p_val_dqm <- ceiling(mean(climdex.r99ptot(ci_val_dqm), na.rm = TRUE))
r99p_val_qdm <- ceiling(mean(climdex.r99ptot(ci_val_qdm), na.rm = TRUE))
v.r99p <- c(r99p_val_raw, r99p_val_obs, r99p_val_rqm, r99p_val_gpm, 
            r99p_val_eqm, r99p_val_dqm, r99p_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum 1-day precipitation. rx1day[mm/1day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx1day_val_raw <- ceiling(max(climdex.rx1day(ci_val_raw, freq = c("monthly")), na.rm = TRUE))
rx1day_val_obs <- ceiling(max(climdex.rx1day(ci_val_obs, freq = c("monthly")), na.rm = TRUE))
rx1day_val_rqm <- ceiling(max(climdex.rx1day(ci_val_rqm, freq = c("monthly")), na.rm = TRUE))
rx1day_val_gpm <- ceiling(max(climdex.rx1day(ci_val_gpm, freq = c("monthly")), na.rm = TRUE))
rx1day_val_eqm <- ceiling(max(climdex.rx1day(ci_val_eqm, freq = c("monthly")), na.rm = TRUE))
rx1day_val_dqm <- ceiling(max(climdex.rx1day(ci_val_dqm, freq = c("monthly")), na.rm = TRUE))
rx1day_val_qdm <- ceiling(max(climdex.rx1day(ci_val_qdm, freq = c("monthly")), na.rm = TRUE))
v.rx1day <- c(rx1day_val_raw, rx1day_val_obs, rx1day_val_rqm, rx1day_val_gpm, 
              rx1day_val_eqm, rx1day_val_dqm, rx1day_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum consecutive 5-day precipitation. Rx5day [mm/5day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx5day_val_raw <- ceiling(max(climdex.rx5day(ci_val_raw, freq = c("monthly")), na.rm = TRUE))
rx5day_val_obs <- ceiling(max(climdex.rx5day(ci_val_obs, freq = c("monthly")), na.rm = TRUE))
rx5day_val_rqm <- ceiling(max(climdex.rx5day(ci_val_rqm, freq = c("monthly")), na.rm = TRUE))
rx5day_val_gpm <- ceiling(max(climdex.rx5day(ci_val_gpm, freq = c("monthly")), na.rm = TRUE))
rx5day_val_eqm <- ceiling(max(climdex.rx5day(ci_val_eqm, freq = c("monthly")), na.rm = TRUE))
rx5day_val_dqm <- ceiling(max(climdex.rx5day(ci_val_dqm, freq = c("monthly")), na.rm = TRUE))
rx5day_val_qdm <- ceiling(max(climdex.rx5day(ci_val_qdm, freq = c("monthly")), na.rm = TRUE))
v.rx5day <- c(rx5day_val_raw, rx5day_val_obs, rx5day_val_rqm, rx5day_val_gpm, 
              rx5day_val_eqm, rx5day_val_dqm, rx5day_val_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Simple precipitation intensity index sdii [mm/day] of wet days ONLY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdii_val_raw <- ceiling(mean(climdex.sdii(ci_val_raw), na.rm = TRUE))
sdii_val_obs <- ceiling(mean(climdex.sdii(ci_val_obs), na.rm = TRUE))
sdii_val_rqm <- ceiling(mean(climdex.sdii(ci_val_rqm), na.rm = TRUE))
sdii_val_gpm <- ceiling(mean(climdex.sdii(ci_val_gpm), na.rm = TRUE))
sdii_val_eqm <- ceiling(mean(climdex.sdii(ci_val_eqm), na.rm = TRUE))
sdii_val_dqm <- ceiling(mean(climdex.sdii(ci_val_dqm), na.rm = TRUE))
sdii_val_qdm <- ceiling(mean(climdex.sdii(ci_val_qdm), na.rm = TRUE))
v.sdii <- c(sdii_val_raw, sdii_val_obs, sdii_val_rqm, sdii_val_gpm, 
            sdii_val_eqm, sdii_val_dqm, sdii_val_qdm)


# A climdex performance data.frame is created
df.perf.climdex <- data.frame(v.prcptot, v.cwd, v.cdd,
                              v.r10mm, v.r20mm, v.r30mm,
                              v.r50mm, v.r95p, v.r99p,
                              v.rx1day, v.rx5day, v.sdii)

# data.frame variables are renamed
names(df.perf.climdex) <- c("prcptot", "cwd", "cdd",
                            "r10mm", "r20mm", "r30mm",
                            "r50mm", "r95p", "r99p",
                            "rx1day", "rx5day", "sdii")

# data.frame row names are renamed
rownames(df.perf.climdex) <- c("raw", "obs", "rqm", "gpm", "eqm", "dqm", "qdm")

# ID variables are incorporated
df.perf.climdex$rcm_id <- rcm.id
df.perf.climdex$rcp_id <- c("hist")
df.perf.climdex$sta_id <- sta.id 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum 1-day precipitation. rx1day[mm/1day]. Annual Max.ONLY!!!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx1day_valYEAR_raw <- ceiling((climdex.rx1day(ci_val_raw, freq = c("annual"))))
rx1day_valYEAR_obs <- ceiling((climdex.rx1day(ci_val_obs, freq = c("annual"))))
rx1day_valYEAR_rqm <- ceiling((climdex.rx1day(ci_val_rqm, freq = c("annual"))))
rx1day_valYEAR_gpm <- ceiling((climdex.rx1day(ci_val_gpm, freq = c("annual"))))
rx1day_valYEAR_eqm <- ceiling((climdex.rx1day(ci_val_eqm, freq = c("annual"))))
rx1day_valYEAR_dqm <- ceiling((climdex.rx1day(ci_val_dqm, freq = c("annual"))))
rx1day_valYEAR_qdm <- ceiling((climdex.rx1day(ci_val_qdm, freq = c("annual"))))

# A data.frame is created
v.rx1day_YEAR <- data.frame(rx1day_valYEAR_raw, rx1day_valYEAR_obs, rx1day_valYEAR_rqm, rx1day_valYEAR_gpm, 
                            rx1day_valYEAR_eqm, rx1day_valYEAR_dqm, rx1day_valYEAR_qdm)

# data.frame row names are renamed
names(v.rx1day_YEAR) <- c("raw", "obs", "rqm", "gpm", "eqm", "dqm", "qdm")

# melt {reshape} function is requested to convert data from wide
# to long format 
df.rx1day_YEAR.long <- reshape::melt(v.rx1day_YEAR)

# data.frame row names are renamed
names(df.rx1day_YEAR.long) <- c("BC", "AMP24")

# ID variables are incorporated
df.rx1day_YEAR.long$rcm_id <- rcm.id
df.rx1day_YEAR.long$rcp_id <- c("hist")
df.rx1day_YEAR.long$sta_id <- sta.id 

# A future simulated monthly subset is executed
df.rx1day_YEAR.long.subset.hist <- subset(df.rx1day_YEAR.long, BC != "raw")

# Maximum monthly AMPs24 are extracted for future simulations
max.sim.fut.rqm <- tapply(bc_sim_rqm_qmap$prec, bc_sim_rqm_qmap$YEAR, max)
max.sim.fut.gpm <- tapply(bc_sim_gpm_dR$prec, bc_sim_gpm_dR$YEAR, max)
max.sim.fut.eqm <- tapply(bc_sim_eqm_dR$prec, bc_sim_eqm_dR$YEAR, max)
max.sim.fut.dqm <- tapply(bc_sim_dqm_dR$prec, bc_sim_dqm_dR$YEAR, max)
max.sim.fut.qdm <- tapply(bc_sim_qdm_dR$prec, bc_sim_qdm_dR$YEAR, max)

# A data.frame is created
df.rx1day_YEAR_fut_sim <- data.frame(max.sim.fut.rqm,
                                     max.sim.fut.gpm,
                                     max.sim.fut.eqm,
                                     max.sim.fut.dqm,
                                     max.sim.fut.qdm)

# data.frame row names are renamed
names(df.rx1day_YEAR_fut_sim) <- c("rqm", "gpm", "eqm", "dqm", "qdm")

# melt {reshape} function is requested to convert data from wide
# to long format 
df.rx1day_YEAR_fut_sim.long <- reshape::melt(df.rx1day_YEAR_fut_sim)

# data.frame row names are renamed
names(df.rx1day_YEAR_fut_sim.long) <- c("BC", "AMP24")

# ID variables are incorporated
df.rx1day_YEAR_fut_sim.long$rcm_id <- rcm.id
df.rx1day_YEAR_fut_sim.long$rcp_id <- c("fut")
df.rx1day_YEAR_fut_sim.long$sta_id <- sta.id 

# An rbind data.frame is created
df.rx1day_YEAR_hist_fut_comp <- rbind(df.rx1day_YEAR.long.subset.hist,
                                      df.rx1day_YEAR_fut_sim.long)

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Future climdex.prcptot {climdex.pcic} Total Daily Precipitation implementation of the ETCCDI climate change indices
# ----------------------------------------------------------------------------------------------------------------------------------

# Dates are parsed into PCICt class
prec_dates_fut <- as.PCICt(as.character(df.base.future$date), cal="gregorian")

# climdexInput object from vectors of data are created
ci_fut_raw <- climdexInput.raw(prec = df.base.future$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# climdexInput object from vectors of data are created
ci_fut_obs <- climdexInput.raw(prec = df_obs_1991_2005$prec,
                               prec.dates = prec_dates_val,
                               base.range=c(1991, 2005))

# climdexInput object from vectors of data are created
ci_fut_rqm <- climdexInput.raw(prec = bc_sim_rqm_qmap$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# climdexInput object from vectors of data are created
ci_fut_gpm <- climdexInput.raw(prec = bc_sim_gpm_dR$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# climdexInput object from vectors of data are created
ci_fut_eqm <- climdexInput.raw(prec = bc_sim_eqm_dR$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# climdexInput object from vectors of data are created
ci_fut_dqm <- climdexInput.raw(prec = bc_sim_dqm_dR$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# climdexInput object from vectors of data are created
ci_fut_qdm <- climdexInput.raw(prec = bc_sim_qdm_dR$prec,
                               prec.dates = prec_dates_fut,
                               base.range=c(2006, 2099))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# climdex.prcptot {climdex.pcic} Total Daily Prec. prcptot[mm/year]
# ETCCDI Climate Change Indices. Definitions of the 27 core indices
# http://etccdi.pacificclimate.org/list_27_indices.shtml
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total prec.wet days [mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prcptot_fut_raw <- ceiling(mean(climdex.prcptot(ci_fut_raw), na.rm = TRUE))
prcptot_fut_obs <- ceiling(mean(climdex.prcptot(ci_fut_obs), na.rm = TRUE))
prcptot_fut_rqm <- ceiling(mean(climdex.prcptot(ci_fut_rqm), na.rm = TRUE))
prcptot_fut_gpm <- ceiling(mean(climdex.prcptot(ci_fut_gpm), na.rm = TRUE))
prcptot_fut_eqm <- ceiling(mean(climdex.prcptot(ci_fut_eqm), na.rm = TRUE))
prcptot_fut_dqm <- ceiling(mean(climdex.prcptot(ci_fut_dqm), na.rm = TRUE))
prcptot_fut_qdm <- ceiling(mean(climdex.prcptot(ci_fut_qdm), na.rm = TRUE))
v.prcptot <- c(prcptot_fut_raw, prcptot_fut_obs, prcptot_fut_rqm, prcptot_fut_gpm, 
               prcptot_fut_eqm, prcptot_fut_dqm, prcptot_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive wet days cwd RR ≥ 1mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd_fut_raw <- ceiling(mean(climdex.cwd(ci_fut_raw), na.rm = TRUE))
cwd_fut_obs <- ceiling(mean(climdex.cwd(ci_fut_obs), na.rm = TRUE))
cwd_fut_rqm <- ceiling(mean(climdex.cwd(ci_fut_rqm), na.rm = TRUE))
cwd_fut_gpm <- ceiling(mean(climdex.cwd(ci_fut_gpm), na.rm = TRUE))
cwd_fut_eqm <- ceiling(mean(climdex.cwd(ci_fut_eqm), na.rm = TRUE))
cwd_fut_dqm <- ceiling(mean(climdex.cwd(ci_fut_dqm), na.rm = TRUE))
cwd_fut_qdm <- ceiling(mean(climdex.cwd(ci_fut_qdm), na.rm = TRUE))
v.cwd <- c(cwd_fut_raw, cwd_fut_obs, cwd_fut_rqm, cwd_fut_gpm, 
           cwd_fut_eqm, cwd_fut_dqm, cwd_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Consecutive dry days cdd [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdd_fut_raw <- ceiling(mean(climdex.cdd(ci_fut_raw), na.rm = TRUE))
cdd_fut_obs <- ceiling(mean(climdex.cdd(ci_fut_obs), na.rm = TRUE))
cdd_fut_rqm <- ceiling(mean(climdex.cdd(ci_fut_rqm), na.rm = TRUE))
cdd_fut_gpm <- ceiling(mean(climdex.cdd(ci_fut_gpm), na.rm = TRUE))
cdd_fut_eqm <- ceiling(mean(climdex.cdd(ci_fut_eqm), na.rm = TRUE))
cdd_fut_dqm <- ceiling(mean(climdex.cdd(ci_fut_dqm), na.rm = TRUE))
cdd_fut_qdm <- ceiling(mean(climdex.cdd(ci_fut_qdm), na.rm = TRUE))
v.cdd <- c(cdd_fut_raw, cdd_fut_obs, cdd_fut_rqm, cdd_fut_gpm, 
           cdd_fut_eqm, cdd_fut_dqm, cdd_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 10mm r10mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r10mm_fut_raw <- ceiling(mean(climdex.rnnmm(ci_fut_raw, threshold = 10), na.rm = TRUE))
r10mm_fut_obs <- ceiling(mean(climdex.rnnmm(ci_fut_obs, threshold = 10), na.rm = TRUE))
r10mm_fut_rqm <- ceiling(mean(climdex.rnnmm(ci_fut_rqm, threshold = 10), na.rm = TRUE))
r10mm_fut_gpm <- ceiling(mean(climdex.rnnmm(ci_fut_gpm, threshold = 10), na.rm = TRUE))
r10mm_fut_eqm <- ceiling(mean(climdex.rnnmm(ci_fut_eqm, threshold = 10), na.rm = TRUE))
r10mm_fut_dqm <- ceiling(mean(climdex.rnnmm(ci_fut_dqm, threshold = 10), na.rm = TRUE))
r10mm_fut_qdm <- ceiling(mean(climdex.rnnmm(ci_fut_qdm, threshold = 10), na.rm = TRUE))
v.r10mm <- c(r10mm_fut_raw, r10mm_fut_obs, r10mm_fut_rqm, r10mm_fut_gpm, 
             r10mm_fut_eqm, r10mm_fut_dqm, r10mm_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 20mm r20mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r20mm_fut_raw <- ceiling(mean(climdex.rnnmm(ci_fut_raw, threshold = 20), na.rm = TRUE))
r20mm_fut_obs <- ceiling(mean(climdex.rnnmm(ci_fut_obs, threshold = 20), na.rm = TRUE))
r20mm_fut_rqm <- ceiling(mean(climdex.rnnmm(ci_fut_rqm, threshold = 20), na.rm = TRUE))
r20mm_fut_gpm <- ceiling(mean(climdex.rnnmm(ci_fut_gpm, threshold = 20), na.rm = TRUE))
r20mm_fut_eqm <- ceiling(mean(climdex.rnnmm(ci_fut_eqm, threshold = 20), na.rm = TRUE))
r20mm_fut_dqm <- ceiling(mean(climdex.rnnmm(ci_fut_dqm, threshold = 20), na.rm = TRUE))
r20mm_fut_qdm <- ceiling(mean(climdex.rnnmm(ci_fut_qdm, threshold = 20), na.rm = TRUE))
v.r20mm <- c(r20mm_fut_raw, r20mm_fut_obs, r20mm_fut_rqm, r20mm_fut_gpm, 
             r20mm_fut_eqm, r20mm_fut_dqm, r20mm_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 30mm r30mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r30mm_fut_raw <- ceiling(mean(climdex.rnnmm(ci_fut_raw, threshold = 30), na.rm = TRUE))
r30mm_fut_obs <- ceiling(mean(climdex.rnnmm(ci_fut_obs, threshold = 30), na.rm = TRUE))
r30mm_fut_rqm <- ceiling(mean(climdex.rnnmm(ci_fut_rqm, threshold = 30), na.rm = TRUE))
r30mm_fut_gpm <- ceiling(mean(climdex.rnnmm(ci_fut_gpm, threshold = 30), na.rm = TRUE))
r30mm_fut_eqm <- ceiling(mean(climdex.rnnmm(ci_fut_eqm, threshold = 30), na.rm = TRUE))
r30mm_fut_dqm <- ceiling(mean(climdex.rnnmm(ci_fut_dqm, threshold = 30), na.rm = TRUE))
r30mm_fut_qdm <- ceiling(mean(climdex.rnnmm(ci_fut_qdm, threshold = 30), na.rm = TRUE))
v.r30mm <- c(r30mm_fut_raw, r30mm_fut_obs, r30mm_fut_rqm, r30mm_fut_gpm, 
             r30mm_fut_eqm, r30mm_fut_dqm, r30mm_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual count of days when PRCP≥ 50mm r50mm [days/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r50mm_fut_raw <- ceiling(mean(climdex.rnnmm(ci_fut_raw, threshold = 50), na.rm = TRUE))
r50mm_fut_obs <- ceiling(mean(climdex.rnnmm(ci_fut_obs, threshold = 50), na.rm = TRUE))
r50mm_fut_rqm <- ceiling(mean(climdex.rnnmm(ci_fut_rqm, threshold = 50), na.rm = TRUE))
r50mm_fut_gpm <- ceiling(mean(climdex.rnnmm(ci_fut_gpm, threshold = 50), na.rm = TRUE))
r50mm_fut_eqm <- ceiling(mean(climdex.rnnmm(ci_fut_eqm, threshold = 50), na.rm = TRUE))
r50mm_fut_dqm <- ceiling(mean(climdex.rnnmm(ci_fut_dqm, threshold = 50), na.rm = TRUE))
r50mm_fut_qdm <- ceiling(mean(climdex.rnnmm(ci_fut_qdm, threshold = 50), na.rm = TRUE))
v.r50mm <- c(r50mm_fut_raw, r50mm_fut_obs, r50mm_fut_rqm, r50mm_fut_gpm, 
             r50mm_fut_eqm, r50mm_fut_dqm, r50mm_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 95p r95p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r95p_fut_raw <- ceiling(mean(climdex.r95ptot(ci_fut_raw), na.rm = TRUE))
r95p_fut_obs <- ceiling(mean(climdex.r95ptot(ci_fut_obs), na.rm = TRUE))
r95p_fut_rqm <- ceiling(mean(climdex.r95ptot(ci_fut_rqm), na.rm = TRUE))
r95p_fut_gpm <- ceiling(mean(climdex.r95ptot(ci_fut_gpm), na.rm = TRUE))
r95p_fut_eqm <- ceiling(mean(climdex.r95ptot(ci_fut_eqm), na.rm = TRUE))
r95p_fut_dqm <- ceiling(mean(climdex.r95ptot(ci_fut_dqm), na.rm = TRUE))
r95p_fut_qdm <- ceiling(mean(climdex.r95ptot(ci_fut_qdm), na.rm = TRUE))
v.r95p <- c(r95p_fut_raw, r95p_fut_obs, r95p_fut_rqm, r95p_fut_gpm, 
            r95p_fut_eqm, r95p_fut_dqm, r95p_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Annual total PRCP when RR > 99p r99p[mm/year]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r99p_fut_raw <- ceiling(mean(climdex.r99ptot(ci_fut_raw), na.rm = TRUE))
r99p_fut_obs <- ceiling(mean(climdex.r99ptot(ci_fut_obs), na.rm = TRUE))
r99p_fut_rqm <- ceiling(mean(climdex.r99ptot(ci_fut_rqm), na.rm = TRUE))
r99p_fut_gpm <- ceiling(mean(climdex.r99ptot(ci_fut_gpm), na.rm = TRUE))
r99p_fut_eqm <- ceiling(mean(climdex.r99ptot(ci_fut_eqm), na.rm = TRUE))
r99p_fut_dqm <- ceiling(mean(climdex.r99ptot(ci_fut_dqm), na.rm = TRUE))
r99p_fut_qdm <- ceiling(mean(climdex.r99ptot(ci_fut_qdm), na.rm = TRUE))
v.r99p <- c(r99p_fut_raw, r99p_fut_obs, r99p_fut_rqm, r99p_fut_gpm, 
            r99p_fut_eqm, r99p_fut_dqm, r99p_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum 1-day precipitation. rx1day[mm/1day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx1day_fut_raw <- ceiling(max(climdex.rx1day(ci_fut_raw, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_obs <- ceiling(max(climdex.rx1day(ci_fut_obs, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_rqm <- ceiling(max(climdex.rx1day(ci_fut_rqm, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_gpm <- ceiling(max(climdex.rx1day(ci_fut_gpm, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_eqm <- ceiling(max(climdex.rx1day(ci_fut_eqm, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_dqm <- ceiling(max(climdex.rx1day(ci_fut_dqm, freq = c("monthly")), na.rm = TRUE))
rx1day_fut_qdm <- ceiling(max(climdex.rx1day(ci_fut_qdm, freq = c("monthly")), na.rm = TRUE))
v.rx1day <- c(rx1day_fut_raw, rx1day_fut_obs, rx1day_fut_rqm, rx1day_fut_gpm, 
              rx1day_fut_eqm, rx1day_fut_dqm, rx1day_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Monthly maximum consecutive 5-day precipitation. Rx5day [mm/5day]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rx5day_fut_raw <- ceiling(max(climdex.rx5day(ci_fut_raw, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_obs <- ceiling(max(climdex.rx5day(ci_fut_obs, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_rqm <- ceiling(max(climdex.rx5day(ci_fut_rqm, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_gpm <- ceiling(max(climdex.rx5day(ci_fut_gpm, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_eqm <- ceiling(max(climdex.rx5day(ci_fut_eqm, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_dqm <- ceiling(max(climdex.rx5day(ci_fut_dqm, freq = c("monthly")), na.rm = TRUE))
rx5day_fut_qdm <- ceiling(max(climdex.rx5day(ci_fut_qdm, freq = c("monthly")), na.rm = TRUE))
v.rx5day <- c(rx5day_fut_raw, rx5day_fut_obs, rx5day_fut_rqm, rx5day_fut_gpm, 
              rx5day_fut_eqm, rx5day_fut_dqm, rx5day_fut_qdm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANNUAL
# Simple precipitation intensity index sdii [mm/day] of wet days ONLY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdii_fut_raw <- ceiling(mean(climdex.sdii(ci_fut_raw), na.rm = TRUE))
sdii_fut_obs <- ceiling(mean(climdex.sdii(ci_fut_obs), na.rm = TRUE))
sdii_fut_rqm <- ceiling(mean(climdex.sdii(ci_fut_rqm), na.rm = TRUE))
sdii_fut_gpm <- ceiling(mean(climdex.sdii(ci_fut_gpm), na.rm = TRUE))
sdii_fut_eqm <- ceiling(mean(climdex.sdii(ci_fut_eqm), na.rm = TRUE))
sdii_fut_dqm <- ceiling(mean(climdex.sdii(ci_fut_dqm), na.rm = TRUE))
sdii_fut_qdm <- ceiling(mean(climdex.sdii(ci_fut_qdm), na.rm = TRUE))
v.sdii <- c(sdii_fut_raw, sdii_fut_obs, sdii_fut_rqm, sdii_fut_gpm, 
            sdii_fut_eqm, sdii_fut_dqm, sdii_fut_qdm)


# A climdex performance data.frame is created
df.perf.climdex.fut <- data.frame(v.prcptot, v.cwd, v.cdd,
                                  v.r10mm, v.r20mm, v.r30mm,
                                  v.r50mm, v.r95p, v.r99p,
                                  v.rx1day, v.rx5day, v.sdii)

# data.frame variables are renamed
names(df.perf.climdex.fut) <- c("prcptot", "cwd", "cdd",
                                "r10mm", "r20mm", "r30mm",
                                "r50mm", "r95p", "r99p",
                                "rx1day", "rx5day", "sdii")

# data.frame row names are renamed
rownames(df.perf.climdex.fut) <- c("raw", "obs", "rqm", "gpm", "eqm", "dqm", "qdm")

# ID variables are incorporated
df.perf.climdex.fut$rcm_id <- rcm.id
df.perf.climdex.fut$rcp_id <- c("rcp85")
df.perf.climdex.fut$sta_id <- sta.id 

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK:  ECDF + graphical analysis
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A new scale color is defined
scale01 <- c("#a6761d","#ca0020","#4dac26", "#404040", "#f4a582", "#0571b0")  
scale01 <- c ("#377eb8","#e41a1c", "#4daf4a", "black", "#ff7f00", "blue", "orange")
scale01 <- c ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
scale01 <- c ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#666666")

# Bias corrected data.frames are selected
df.ecdf.obs <- df_obs_1991_2005  # 5479
df.ecdf.rqm <- bc_hist_rqm_qmap  # 5479
df.ecdf.gpm <- bc_hist_gpm_dR    # 5479
df.ecdf.eqm <- bc_hist_eqm_dR    # 5479
df.ecdf.dqm <- bc_hist_dqm_dR    # 5479
df.ecdf.qdm <- bc_hist_qdm_dR    # 5479
df.ecdf.raw <- df.hindcast.val   # 5479

# ID variables are incorporated
df.ecdf.obs$method <- c("obs") 
df.ecdf.rqm$method <- c("rqm")
df.ecdf.gpm$method <- c("gpm")
df.ecdf.eqm$method <- c("eqm")
df.ecdf.dqm$method <- c("dqm")
df.ecdf.qdm$method <- c("qdm")
df.ecdf.raw$method <- c("raw")

# {lubridate} functions are applied to create new variables
df.ecdf.obs$MONTH <- month(df.ecdf.obs$date)
df.ecdf.rqm$MONTH <- month(df.ecdf.rqm$date)
df.ecdf.gpm$MONTH <- month(df.ecdf.gpm$date)
df.ecdf.eqm$MONTH <- month(df.ecdf.eqm$date)
df.ecdf.dqm$MONTH <- month(df.ecdf.dqm$date)
df.ecdf.qdm$MONTH <- month(df.ecdf.qdm$date)
df.ecdf.raw$MONTH <- month(df.ecdf.raw$date)

# Variable names are reordered
names(df.ecdf.obs) <- c("dates", "prec", "MONTH", "YEAR","method")
names(df.ecdf.raw) <- c("dates", "prec", "MONTH", "YEAR","method")

# A compiled data.frame is created
df.ecdf.total <-rbind(df.ecdf.obs,
                      df.ecdf.rqm,
                      df.ecdf.gpm,
                      df.ecdf.eqm,
                      df.ecdf.dqm,
                      df.ecdf.qdm,
                      df.ecdf.raw)

# Historical observations are repeated for dispersion plots
df.ecdf.total$obs <- rep(df.ecdf.obs$prec, 7)

# Months component of a date-time as character are added
df.ecdf.total$MONTH_CH <- lubridate::month(df.ecdf.total$dates, label = TRUE) 

# Character month vector is created
v.select.month.low <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Months are transformed and ordered as factors
df.ecdf.total$MONTH_CH <- factor(df.ecdf.total$MONTH_CH, levels = v.select.month.low)

# Daily precipitation is aggregated by Month and BC method
df.ecdf.total.aggregated <- aggregate(prec ~ MONTH_CH + method + MONTH,
                                      df.ecdf.total,
                                      FUN = sum)

# Monthly precipitation is divided by length of validation period
df.ecdf.total.aggregated$prec <- (df.ecdf.total.aggregated$prec) / val.length

# ID variables are incorporated
df.ecdf.total.aggregated$rcm_id <- rcm.id
df.ecdf.total.aggregated$rcp_id <- c("hist")
df.ecdf.total.aggregated$sta_id <- sta.id

# ID variables are incorporated
df.ecdf.total$rcm_id <- rcm.id
df.ecdf.total$rcp_id <- c("hist")
df.ecdf.total$sta_id <- sta.id

# Daily precipitation is aggregated by Month and BC method
df.hindcast.val.aggregated <- aggregate(prec ~ MONTH,
                                        df.hindcast.val,
                                        FUN = sum)

# Monthly precipitation is divided by length of validation period
df.hindcast.val.aggregated$prec <- (df.hindcast.val.aggregated$prec) / val.length

# QQ plot dispersion observations vs. simulations is created
ggplot() +
  geom_point(aes(x = obs,y = prec,shape = method,colour = method),data=df.ecdf.total,size = 1.8) +
  facet_wrap(facets = ~MONTH_CH, scales = 'free') +
  geom_abline(data=df.ecdf.total,colour = '#666666',size = 0.55,linetype = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_color_manual(values = scale01) +
  scale_shape_manual(values = c(2,13,4,5,6,7,8,9,10,11,12,13)) +
  ggtitle(paste(sta.id, rcm.id, "QQ plot dispersion observations vs. simulations for historical validation period")) +
  xlab("Observed precipitation (mm/day)") +
  ylab("Simulated precipitation (mm/day)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 90.0), text=element_text(size=14, family="serif"))

# A time-series BC vs. Observation for validation period is created for monthly totals
ggplot() +
  geom_point(aes(x = MONTH_CH,y = prec,shape = method,colour = method),data=df.ecdf.total.aggregated,size = 2.2) +
  geom_line(aes(x = MONTH,y = prec,colour = method,linetype = method),data=df.ecdf.total.aggregated,size = 0.85) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_color_manual(values = scale01) +
  scale_shape_manual(values = c(2,13,4,19,6,7,8,9,10,11,12,13)) +
  ggtitle(paste(sta.id, rcm.id, "Bias correction methods vs. observations for historical validation period")) +
  xlab("Month of the year") +
  ylab("Precipitation (mm/month)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# A boxplot BC vs. Observation for validation period is created for monthly totals
ggplot() +
  geom_boxplot(aes(y = prec,x = method,colour = method),data=df.ecdf.total) +
  #facet_wrap(facets = ~ MONTH_CH) +
  scale_color_manual(values = scale01) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0), limits = c(0, 600)) +
  ggtitle(paste(sta.id, rcm.id, "Bias correction methods vs. observations for historical validation period")) +
  xlab("Bias correction method (BC)") +
  ylab("Precipitation (mm/day)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# ecdf Empirical cumulative distribution functions for historical validation period
ggplot() +
  stat_ecdf(aes(x = prec,colour = method, linetype = method), size = 0.85,data=df.ecdf.total) +
  facet_wrap(facets = ~ MONTH_CH, scales = 'free') +
  scale_color_manual(values = scale01) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0), limits = c(0, 200)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, rcm.id, "Empirical cumulative distribution functions for historical validation period")) +
  xlab("Precipitation (mm/day)") +
  ylab("Cumulative distribution (fraction)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# ecdf Empirical cumulative distribution functions for historical validation period for MONTH
ggplot() +
  stat_ecdf(aes(x = prec,colour = method, linetype = method), size = 0.85, data = subset(df.ecdf.total, MONTH == 8)) +
  scale_color_manual(values = scale01) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0), limits = c(0, 200)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, rcm.id, "Empirical cumulative distribution functions for historical validation period")) +
  xlab("Precipitation (mm/day)") +
  ylab("Cumulative distribution (fraction)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Obs and Raw historical vs. BC methods for AMP24 ONLY 
ggplot() +
  geom_boxplot(aes(y = AMP24,x = BC,colour = BC),data=df.rx1day_YEAR.long, size = 0.85) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_color_manual(values = scale01) +
  scale_shape_manual(values = c(2,13,4,19,6,7,8,9,10,11,12,13)) +
  ggtitle(paste(sta.id, rcm.id, "Obs and Raw historical vs. BC methods for AMP24 ONLY")) +
  xlab("Bias Correction Method (BC)") +
  ylab("AMP24 Precipitation (mm/month)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Obs and Raw historical vs. BC methods for AMP24 ONLY 
ggplot() +
  geom_boxplot(aes(y = AMP24,x = BC,colour = BC, fill = rcp_id), data = df.rx1day_YEAR_hist_fut_comp,
               size = 0.85, outlier.shape = 19, outlier.size = 4) +
  geom_hline(data=df.rx1day_YEAR.long,yintercept = max(rx1day_valYEAR_obs),size = 0.55,linetype = 2) +
  geom_hline(data=df.rx1day_YEAR.long,yintercept = mean(rx1day_valYEAR_obs),size = 0.55,linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  stat_summary(aes(x = BC,y = AMP24,fill = rcp_id,linetype = rcp_id,shape = rcp_id),
               data=df.rx1day_YEAR_hist_fut_comp,size = 0.75,
               fun=mean,fun.args = list(mult = 1), position=position_dodge(width=0.75)) +
  scale_color_manual(values = scale01) +
  scale_fill_manual(values=c("#999999", "white", "#E69F00", "#56B4E9")) +
  scale_shape_manual(values = c(2,13,4,19,6,7,8,9,10,11,12,13)) +
  ggtitle(paste(sta.id, rcm.id, "Obs and Raw historical vs. BC methods for AMP24 ONLY")) +
  xlab("Bias Correction Method (BC)") +
  ylab("AMP24 Precipitation (mm/month)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Comparative daily 24AMPs per BC method are requested
View(amp_24hr_sim_comp)
View(amp_24hr_hist_comp)
View(amp_24hr_arch_comp)

# Comparative MONTHLY 24AMPs duplicates per BC are requested
View(amp_24hrAMP_hist_comp)
View(df_amp_24hr_sim_duplicated)
View(df_amp_24hr_hist_duplicated)
View(df_amp_24hr_arch_duplicated)

# Bias Correction numerical performance data.frames are requested
View(df.base.desc)
View(df.perf.desc.24hrAMP)
View(df.perf.desc.daily)
View(df.perf.desc.daily.baseline)
View(df.perf.climdex)
View(df.perf.climdex.baseline)
View(df.perf.climdex.fut)
View(df.rx1day_YEAR.long)
View(df.rx1day_YEAR_hist_fut_comp)

# Bias Correction graphical performance data.frames are requested
View(df.ecdf.total.aggregated)
View(df.ecdf.total)

# An output file name is concatenated
t.output01 <- paste("BRICK03_matrix_bc_24h_par_",rcm.id,"_",rcp.id,"_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK03_matrix_bc_24h_par_",rcm.id,"_",rcp.id,"_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(amp_24hr_sim_comp,            #1
     amp_24hr_hist_comp,           #2
     amp_24hr_arch_comp,           #3
     amp_24hrAMP_hist_comp,        #4
     df_amp_24hr_sim_duplicated,   #5
     df_amp_24hr_hist_duplicated,  #6
     df_amp_24hr_arch_duplicated,  #7
     df.base.desc,                 #8
     df.perf.desc.24hrAMP,         #9
     df.perf.desc.daily,           #10 
     df.perf.desc.daily.baseline,  #11
     df.perf.climdex,              #12
     df.perf.climdex.baseline,     #13
     df.perf.climdex.fut,          #14
     df.perf.climdex.baseline_full,#15
     df.rx1day_YEAR.long,          #16
     df.rx1day_YEAR_hist_fut_comp, #17
     df.ecdf.total.aggregated,     #18
     df.ecdf.total,                #19
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("amp_24hr_sim_comp"= amp_24hr_sim_comp,                         #1
                   "amp_24hr_hist_comp"= amp_24hr_hist_comp,                       #2
                   "amp_24hr_arch_comp"= amp_24hr_arch_comp,                       #3
                   "amp_24hrAMP_hist_comp"= amp_24hrAMP_hist_comp,                 #4
                   "df_amp_24hr_sim_duplicated"= df_amp_24hr_sim_duplicated,       #5
                   "df_amp_24hr_hist_duplicated"= df_amp_24hr_hist_duplicated,     #6
                   "df_amp_24hr_arch_duplicated"= df_amp_24hr_arch_duplicated,     #7
                   "df_base_desc"= df.base.desc,                                   #8
                   "df_perf.desc_24hrAMP"= df.perf.desc.24hrAMP,                   #9
                   "df_perf_desc_daily"= df.perf.desc.daily,                       #10
                   "df_perf_desc_daily_baseline" = df.perf.desc.daily.baseline,    #11
                   "df_perf_climdex"= df.perf.climdex,                             #12
                   "df_perf_climdex_baseline" = df.perf.climdex.baseline,          #13
                   "df_perf_climdex_fut" = df.perf.climdex.fut,                    #14
                   "df_perf_climdex_baseline_full" = df.perf.climdex.baseline_full,#15
                   "df_rx1day_YEAR_long" = df.rx1day_YEAR.long,                    #16
                   "df_rx1day_YEAR_hist_fut_comp" = df.rx1day_YEAR_hist_fut_comp,  #17
                   "df_ecdf_total_aggregated" = df.ecdf.total.aggregated,          #18
                   "df_ecdf_total" = df.ecdf.total)                                #19
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
