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

#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:

# BRICK04_AMP24_PDFs_par_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing dataset containing:
# bc_sim_chosen_arch       #1 Yearly AMP24-BC prec.[mm/day] for 1970-2005
# bc_sim_chosen_fut        #2 Yearly AMP24-BC prec.[mm/day] for 2020-2099
# df.merge.metric.arch     #3 Yearly AMP24-BC KS/AD/CvM/AIC/BIC for 1970-2005
# df.merge.metric.sim      #4 Yearly AMP24-BC KS/AD/CvM/AIC/BIC for 2020-2099
# df.merge.metric.arch.par #5 Yearly AMP24-BC Location/Scale/Shape par. for 1970-2005
# df.merge.metric.sim.par  #6 early AMP24-BC Location/Scale/Shape par. for 2020-2099
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
require(tdr)
require(dplyr)
require(ggplot2)
require(ggradar)
require(lubridate)
require(plyr)
require(pastecs)
require(IDFtool)
require(openxlsx)
require(evd)
require(nsRFA)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Pre-processed binary datasets are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# RData objects are loaded from BRICK03
load("BRICK03_matrix_bc_24h_par_HadGEM2_RegCM47_RCP85_AIJS.RData")

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

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK:  Daily Annual Maximum Precipitation (AMP) 
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

# A temporal subset is defined for period 2020-2099
bc_sim_chosen_fut <- amp_24hr_sim_comp[amp_24hr_sim_comp$date >= "2019-01-01" & amp_24hr_sim_comp$date <= "2099-12-31", ]
bc_sim_chosen_arch <- amp_24hr_arch_comp # 1970-2005

# data.frame variables are transformed
bc_sim_chosen_fut$dates <- as.numeric(bc_sim_chosen_fut$dates)   # 2020-2099
bc_sim_chosen_arch$dates <- as.numeric(bc_sim_chosen_arch$dates) # 1970-2005

# -------------------------------------------------------------------------------------------------
# A for-in loop is defined to generate AMP parameters
# -------------------------------------------------------------------------------------------------

# A list container is erased
list.merge.arch <- NULL
list.merge.sim <- NULL
list.merge.arch.par <- NULL
list.merge.sim.par <- NULL

# List containers are defined
list.merge.arch <- list()
list.merge.sim <- list()
list.merge.arch.par <- list()
list.merge.sim.par <- list()

# ============================
# Middle loop is initialized
# ============================
for (n in 2:7) { # length = 6: bc-methods
  
  # data.frame containers are erased
  matrix.base.sim <- NULL
  matrix.base.arch <- NULL
  
  # A data.frame is transformed to matrix object
  matrix.base.sim <- as.matrix(bc_sim_chosen_fut[, c(1, n)])
  matrix.base.arch <- as.matrix(bc_sim_chosen_arch[, c(1, n)])
  
  # Names are assigned to matrix columns
  colnames(matrix.base.sim) <- c("year", "1440")
  colnames(matrix.base.arch) <- c("year", "1440")
  
  # A duration vector is created [min]
  v.duration <- c(1440)
  
  # A periods vector is created [years]
  v.periods <- c(2, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100)
  
  # stat.desc {pastecs} function is requested
  df.stat.desc.sim <- as.data.frame(round(stat.desc(matrix.base.sim[ , ], norm = TRUE), 6))
  df.stat.desc.arch <- as.data.frame(round(stat.desc(matrix.base.arch[ , ], norm = TRUE), 6))
  
  # ----------------------------------------------------------------------------------------------------------------------------------
  # SUBBLOCK: Simulated Archival
  # ----------------------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.arch.lp3 <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "log.pearson3", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "log.pearson3", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.arch.lp3$Test.fit) == TRUE) {
    fit.arch.lp3 <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.arch.lp3[,] <- 9999
    names(fit.arch.lp3) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                             "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.arch.lp3 <- as.data.frame(IDF.Auto.arch.lp3$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.arch[ , 2],
                                    dist = c("LP3"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.arch.lp3$AIC <- temp.metric$AIC[1]
  fit.arch.lp3$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.arch.lp3 <- IDF.Auto.arch.lp3$Distribution$`1440 min`$Parameters$para
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.arch.ev1 <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gumbel", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gumbel", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.arch.ev1$Test.fit) == TRUE) {
    fit.arch.ev1 <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.arch.ev1[,] <- 9999
    names(fit.arch.ev1) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                             "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.arch.ev1 <- as.data.frame(IDF.Auto.arch.ev1$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.arch[ , 2],
                                    dist = c("GUMBEL"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.arch.ev1$AIC <- temp.metric$AIC[1]
  fit.arch.ev1$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.arch.ev1 <- IDF.Auto.arch.ev1$Distribution$`1440 min`$Parameters$para
  
  # A dummy shape variable is created for ev1 ONLY
  fit.par.arch.ev1[3] <- 0
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.arch.gev <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gev", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.arch, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gev", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.arch.gev$Test.fit) == TRUE) {
    fit.arch.gev <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.arch.gev[,] <- 9999
    names(fit.arch.gev) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                             "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.arch.gev <- as.data.frame(IDF.Auto.arch.gev$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.arch[ , 2],
                                    dist = c("GEV"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.arch.gev$AIC <- temp.metric$AIC[1]
  fit.arch.gev$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.arch.gev <- IDF.Auto.arch.gev$Distribution$`1440 min`$Parameters$para
  
  # data.frames objects are merged
  df.merge.arch <- rbind(fit.arch.lp3, fit.arch.ev1, fit.arch.gev)
  
  # data.frame row names are renamed
  rownames(df.merge.arch) <- c("lp3", "ev1", "gev")
  
  # data.frames objects are merged
  df.merge.arch.par <- rbind(fit.par.arch.lp3, fit.par.arch.ev1, fit.par.arch.gev)
  
  # data.frame row names are renamed
  rownames(df.merge.arch.par) <- c("lp3", "ev1", "gev")
  
  # ----------------------------------------------------------------------------------------------------------------------------------
  # SUBBLOCK: Simulated Future Projections
  # ----------------------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.sim.lp3 <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "log.pearson3", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "log.pearson3", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.sim.lp3$Test.fit) == TRUE) {
    fit.sim.lp3 <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.sim.lp3[,] <- 9999
    names(fit.sim.lp3) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                            "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.sim.lp3 <- as.data.frame(IDF.Auto.sim.lp3$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.sim[ , 2],
                                    dist = c("LP3"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.sim.lp3$AIC <- temp.metric$AIC[1]
  fit.sim.lp3$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.sim.lp3 <- IDF.Auto.sim.lp3$Distribution$`1440 min`$Parameters$para
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.sim.ev1 <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gumbel", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gumbel", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.sim.ev1$Test.fit) == TRUE) {
    fit.sim.ev1 <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.sim.ev1[,] <- 9999
    names(fit.sim.ev1) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                            "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.sim.ev1 <- as.data.frame(IDF.Auto.sim.ev1$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.sim[ , 2],
                                    dist = c("GUMBEL"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.sim.ev1$AIC <- temp.metric$AIC[1]
  fit.sim.ev1$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.sim.ev1 <- IDF.Auto.sim.ev1$Distribution$`1440 min`$Parameters$para
  
  # A dummy shape variable is created for ev1 ONLY
  fit.par.sim.ev1[3] <- 0
  
  # -------------------------------------------------------------------------------------------------------------------
  # IDFCurve {IDFtool} function is called and an new object is created
  IDF.Auto.sim.gev <- tryCatch( 
    { # if{}
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gev", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = TRUE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }, # otherwise{}
    error = function(e) {
      IDFCurve(Data = matrix.base.sim, Station=sta.id,
               Duration = v.duration, Periods = v.periods, 
               Type = "gev", M.fit = pdf.fit, 
               Plot = FALSE, Strategy = 1, logaxe = "", 
               CI = FALSE, CIpdf = FALSE, iter = 500,
               goodtest = FALSE, Resolution = 600, 
               SAVE = FALSE, name = TRUE)
    }) # End of try
  # -------------------------------------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------------------------------------------------
  # Goodness-of-fit tests are requested
  if (is.null(IDF.Auto.sim.gev$Test.fit) == TRUE) {
    fit.sim.gev <- data.frame(matrix(ncol = 10, nrow = 1))
    fit.sim.gev[,] <- 9999
    names(fit.sim.gev) <- c("KS", "p.value1", "AD", "p.value2", "Omega2",
                            "p.value3", "AIC", "BIC", "AICc", "KIC")
  } else {
    fit.sim.gev <- as.data.frame(IDF.Auto.sim.gev$Test.fit)
  }
  # -------------------------------------------------------------------------------------------------------------------
  
  # MSClaio2008 {nsRFA} is used to recalculate AIC and BIC
  temp.metric <- nsRFA::MSClaio2008(sample = matrix.base.sim[ , 2],
                                    dist = c("GEV"),
                                    crit = c("AIC", "BIC"))
  
  # AIC and BIC are replaced
  fit.sim.gev$AIC <- temp.metric$AIC[1]
  fit.sim.gev$BIC <- temp.metric$BIC[1]
  
  # Distribution parameters are requested
  fit.par.sim.gev <- IDF.Auto.sim.gev$Distribution$`1440 min`$Parameters$para
  
  # data.frames objects are merged
  df.merge.sim <- rbind(fit.sim.lp3, fit.sim.ev1, fit.sim.gev)
  
  # data.frame row names are renamed
  rownames(df.merge.sim) <- c("lp3", "ev1", "gev")
  
  # data.frames objects are merged
  df.merge.sim.par <- rbind(fit.par.sim.lp3, fit.par.sim.ev1, fit.par.sim.gev)
  
  # data.frame row names are renamed
  rownames(df.merge.sim.par) <- c("lp3", "ev1", "gev")
  
  # Rownames are added as data.frame variables
  df.merge.arch$method <- rownames(df.merge.arch)
  df.merge.sim$method <- rownames(df.merge.sim)
  
  # Rownames are added as data.frame variables
  df.merge.arch.par <- as.data.frame(df.merge.arch.par)
  df.merge.sim.par <- as.data.frame(df.merge.sim.par)
  df.merge.arch.par$method <- rownames(df.merge.arch.par)
  df.merge.sim.par$method <- rownames(df.merge.sim.par)
  
  # BC methods are added as data.frame variables
  df.merge.arch$bc <- names(bc_sim_chosen_arch[n])
  df.merge.sim$bc <- names(bc_sim_chosen_fut[n])
  df.merge.arch.par$bc <- names(bc_sim_chosen_arch[n])
  df.merge.sim.par$bc <- names(bc_sim_chosen_fut[n])
  
  # Data containers are filled
  list.merge.arch[[n]] <- df.merge.arch
  list.merge.arch.par[[n]] <- as.data.frame(df.merge.arch.par)
  list.merge.sim[[n]] <- df.merge.sim
  list.merge.sim.par[[n]] <- df.merge.sim.par
  
  # ==============================   
} # Middle loop is closed
# ==============================

# An external list of "date" data.frames is created
df.merge.metric.arch <- ldply(list.merge.arch, data.frame)
df.merge.metric.sim <- ldply(list.merge.sim, data.frame)
df.merge.metric.arch.par <- ldply(list.merge.arch.par, data.frame)
df.merge.metric.sim.par <- ldply(list.merge.sim.par, data.frame)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Objective functions critical values
#-------------------------------------------------------------------------------------------------------------

# data.frame dummy variables are created
df.merge.metric.arch$ruleKS <- NA
df.merge.metric.arch$ruleAD <- NA
df.merge.metric.arch$ruleOmega2 <- NA
df.merge.metric.sim$ruleKS <- NA
df.merge.metric.sim$ruleAD <- NA
df.merge.metric.sim$ruleOmega2 <- NA

# A for-in loop is defined to evaluate KS, AD and CVM TRUE or FALSE
# Kolmogorov-Smirnov (KS) [critical value = 0.20517] REF???
# Anderson-Darling (AD)   [critical value = 2.5018]  REF???
# Cramer-von Mises (CVM)  [critical value = 0.22101] REF???

# Internal counter is defined
counter <- length(df.merge.metric.arch$KS)

# ============================
# Inner loop is initialized
# ============================
for (w in 1:counter) {
  
  # Kolmogorov-Smirnov (KS) [critical value = 0.20517]
  if (df.merge.metric.arch$KS[w] < 0.20517) {
    df.merge.metric.arch$ruleKS[w] = "OK"
  }
  else {
    df.merge.metric.arch$ruleKS[w] = "FAILED"
  }
  
  # Anderson-Darling (AD) [critical value = 2.5018]
  if (df.merge.metric.arch$AD[w] < 2.5018) {
    df.merge.metric.arch$ruleAD[w] = "OK"
  }
  else {
    df.merge.metric.arch$ruleAD[w] = "FAILED"
  }
  
  # Cramer-von Mises (CVM) [critical value = 0.22101]
  if (df.merge.metric.arch$Omega2[w] < 0.22101) {
    df.merge.metric.arch$ruleOmega2[w] = "OK"
  }
  else {
    df.merge.metric.arch$ruleOmega2[w] = "FAILED"
  }  
  
  # Kolmogorov-Smirnov (KS) [critical value = 0.20517]
  if (df.merge.metric.sim$KS[w] < 0.20517) {
    df.merge.metric.sim$ruleKS[w] = "OK"
  }
  else {
    df.merge.metric.sim$ruleKS[w] = "FAILED"
  }
  
  # Anderson-Darling (AD) [critical value = 2.5018]
  if (df.merge.metric.sim$AD[w] < 2.5018) {
    df.merge.metric.sim$ruleAD[w] = "OK"
  }
  else {
    df.merge.metric.sim$ruleAD[w] = "FAILED"
  }
  
  # Cramer-von Mises (CVM) [critical value = 0.22101]
  if (df.merge.metric.sim$Omega2[w] < 0.22101) {
    df.merge.metric.sim$ruleOmega2[w] = "OK"
  }
  else {
    df.merge.metric.sim$ruleOmega2[w] = "FAILED"
  }  
  
  # ==============================
} # Inner loop is closed
# ==============================

# ID variables are incorporated
df.merge.metric.arch$rcm_id <- rcm.id
df.merge.metric.arch$sta_id <- sta.id 
df.merge.metric.arch$rcp_id <- c("hist")

# ID variables are incorporated
df.merge.metric.arch.par$rcm_id <- rcm.id
df.merge.metric.arch.par$sta_id <- sta.id 
df.merge.metric.arch.par$rcp_id <- c("hist")

# ID variables are incorporated
df.merge.metric.sim $rcm_id <- rcm.id
df.merge.metric.sim $sta_id <- sta.id 
df.merge.metric.sim $rcp_id <- c("sim")

# ID variables are incorporated
df.merge.metric.sim.par$rcm_id <- rcm.id
df.merge.metric.sim.par$sta_id <- sta.id 
df.merge.metric.sim.par$rcp_id <- c("sim")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Archival and future BC selected values are requested
View(bc_sim_chosen_arch)
View(bc_sim_chosen_fut)

# PDFs parameter data.frames are requested
View(df.merge.metric.arch)
View(df.merge.metric.sim)
View(df.merge.metric.arch.par)
View(df.merge.metric.sim.par)

# An output file name is concatenated
t.output01 <- paste("BRICK04_AMP24_PDFs_par_",rcm.id,"_",rcp.id,"_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK04_AMP24_PDFs_par_",rcm.id,"_",rcp.id,"_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(bc_sim_chosen_arch,       #1
     bc_sim_chosen_fut,        #2
     df.merge.metric.arch,     #3
     df.merge.metric.sim,      #4
     df.merge.metric.arch.par, #5
     df.merge.metric.sim.par,  #6
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("bc_sim_chosen_arch"= bc_sim_chosen_arch,             #1
                   "bc_sim_chosen_fut"= bc_sim_chosen_fut,               #2
                   "df_merge_metric_arch"= df.merge.metric.arch,         #3
                   "df_merge_metric_sim"= df.merge.metric.sim,           #4
                   "df_merge_metric_arch_par"= df.merge.metric.arch.par, #5
                   "df_merge_metric_sim_par"= df.merge.metric.sim.par)   #6
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
