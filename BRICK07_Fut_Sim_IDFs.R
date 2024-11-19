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
# /////////////////////////////////////////////////////////// /////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INFO: This script is intended for the generation of future IDFs curves on specific Lat/Long position
# of weather-stations for durations [5, 10, 15, 30, 60, 120, 180, 360, 720, 1440] minutes
# using Generalized Extreme Value (GEV) Distribution ONLY for future CORDEX member under RCP85
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
#
# BRICK04_AMP24_PDFs_par_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing dataset containing:
# bc_sim_chosen_arch       #1 Yearly AMP24-BC prec.[mm/day] for 1970-2005
# bc_sim_chosen_fut        #2 Yearly AMP24-BC prec.[mm/day] for 2020-2099
# df.merge.metric.arch     #3 Yearly AMP24-BC KS/AD/CvM/AIC/BIC for 1970-2005
# df.merge.metric.sim      #4 Yearly AMP24-BC KS/AD/CvM/AIC/BIC for 2020-2099
# df.merge.metric.arch.par #5 Yearly AMP24-BC Location/Scale/Shape par. for 1970-2005
# df.merge.metric.sim.par  #6 early AMP24-BC Location/Scale/Shape par. for 2020-2099
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

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
#
# BRICK07_matrix_bc_AMP_future_par_HadGEM2_RegCM47_RCP85_AIJS.RData: R-object containing dataset vectors 
# of future records with the following objects:
# df.int.obs.pdf.adjusted        #1: Hist. AMPs expressed as int.[mm/day] for 1970-2022 pdf-adjusted
# df.int.obs.model.adjusted      #2: Hist. AMPs expressed as int.[mm/day] for 1970-2022 model-adjusted
# df.outer.ratio.idf             #3: Fut.Rel.Ratio.Change for 1970-2022 per BC method
# df.outer.ratio.idf.pdf.adj     #4: Fut.Rel.Ratio.Change for 1970-2022 per BC method pdf.adj
# df.outer.ratio.idf.model.adj   #5: Fut.Rel.Ratio.Change for 1970-2022 per BC method model.adj
# df.outer.int.idf               #6: Fut.AMPs expressed as int.[mm/day] per BC method
# df.outer.int.idf.pdf.adj       #7: Fut.AMPs expressed as int.[mm/day] per BC method pdf.adj
# df.outer.int.idf.model.adj     #8: Fut.AMPs expressed as int.[mm/day] per BC method model.adj
# df.para.int.gev.coeff.select   #9: Fut.Par. A/B/C and perf. and BR2/RMSE for RQM ONLY model.adj
# df.fit.int.gev.select          #10: Fut.AMPs KS/AD/CvM/AIC/BIC for RQM ONLY model.adj

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
require(dplyr)
require(devtools)
require(investr)
require(IDFtool)
require(IDF)
require(DescTools)
require(ggplot2)
require(lubridate)
require(matrixStats)
require(pastecs)
require(tidyr)
require(fixIDF)
require(stringr)
require(nsRFA)
require(EnvStats)
require(openxlsx)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Preprocessed binary datasets are loaded
#-------------------------------------------------------------------------------------------------------------

# RData objects are loaded from BRICK04
load("BRICK04_AMP24_PDFs_par_HadGEM2_RegCM47_RCP85_AIJS.RData")

# RData objects are loaded from BRICK04
load("BRICK06_matrix_bc_AMP_par_HadGEM2_RegCM47_RCP85_AIJS.RData")

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

# PDF fitting method is selected
pdf.fit <- "lmoments"

# A duration vector is created [min]
v.duration <- c(5, 10, 15, 30, 60, 120, 360, 720, 1440)

# A periods vector is created [years]
v.periods <- c(2, 3, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Selection of BC Method for Final-Production
#-------------------------------------------------------------------------------------------------------------

# Selector (m), with length = 6: selection of BC-methods are considered
# Ouputs: "eqm[1]" "rqm[2]" "gpm[3]" "dqm[4]" "qdm[5]" "raw[6]"
# Info contained in vector = print(v.bc.arch)
# BEWARE !!! for the purpose of this work, RQM has been selected for Final-Production
bc.model.selector <- 2
# -------------------------------------------------------------------------------------------------

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Future GEV Function Distribution Execution
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# List containers are created
list.outer.ratio.idf <- NULL
list.outer.int.idf <- NULL
list.outer.dur.par <- NULL

# List containers are defined
list.outer.ratio.idf <- list()
list.outer.int.idf <- list()
list.outer.dur.par <- list()

# ==============================
# Outermost loop is initialized
# ==============================
for(p in 1:2) { # length = 2: IDFCurve either PDF-adjusted or Model-adjusted
  
  # Either PDF-adjusted or Model-adjusted selector is defined
  pfd.selector <- p
  
  # A future simulated monthly subset is executed for GEV only
  # BEWARE !!! LP3 and EV1 are neglected
  df.par.metric.amp.arch <- subset(df.merge.metric.arch.par, method == "gev") # gev ONLY
  df.par.metric.amp.sim <- subset(df.merge.metric.sim.par, method == "gev")   # gev ONLY
  
  # Bias-Corrected vectors are isolated
  v.bc.arch <- df.par.metric.amp.arch$bc
  v.bc.sim <- df.par.metric.amp.sim$bc
  
  # Bias-Corrected vectors are compared
  # Comparison should return TRUE
  v.bc.sim == v.bc.arch
  
  # List containers are created
  list.idf.bc.ratio <- NULL
  list.idf.bc.int <- NULL
  list.idf.dur.par <- NULL
  
  # List containers are defined
  list.idf.bc.ratio <- list()
  list.idf.bc.int <- list()
  list.idf.dur.par <- list()
  
  # ============================
  # Outer loop is initialized
  # ============================
  for(m in 1:length(v.bc.arch)) { # length = 6: selection of BC-methods considered
    
    # BC method is selected
    bc.selector <- v.bc.arch[m]
    
    # BC parameters are subset
    temp.amp.arch <- subset(df.par.metric.amp.arch, bc == bc.selector)
    temp.amp.sim <- subset(df.par.metric.amp.sim, bc == bc.selector)
    
    # BC GEV: LOCATION, SCALE and SHAPE parameters are isolated
    par.bc.gcm.baseline <- as.numeric(temp.amp.arch[1, c(1,2,3)]) # e.g. c(67.589330, 16.02118421,-0.0312708596)
    par.bc.gcm.future <-   as.numeric(temp.amp.sim[1, c(1,2,3)])  # e.g. c(73.457975, 18.34188741,-0.061771582)
    
    # BC future RCP Target days [mm/day] are isolated individually
    rcp.target.day.bc <-  as.vector(bc_sim_chosen_fut[, m+1])
    
    # List containers are created
    list.vol.cont.bc <- NULL
    list.int.cont.bc <- NULL
    
    # List containers are defined
    list.vol.cont.bc <- list()
    list.int.cont.bc <- list()
    
    # ============================
    # Middle loop is initialized
    # ============================
    for(j in 1:length(rcp.target.day.bc)) { # length = 80: AMPs per BC-method for 2020-2099

      #-------------------------------------------------------------------------------------------------------------------
      # SUB-BLOCK: Temporal BC for GCM-RCM ONLY
      #-------------------------------------------------------------------------------------------------------------------
      
      # GCM-RCM Model Future RCP [-]. It MUST yield: probability [-]
      vec.pgev.fut.bc <- pgevd(rcp.target.day.bc[j],
                               location = par.bc.gcm.future[1],
                               scale = par.bc.gcm.future[2],
                               shape = par.bc.gcm.future[3])
      
      # Temporal probability is requested
      #print(vec.pgev.fut.bc)
      
      # GCM-RCM Model Future is reversed with respect to Historical GCM-RCM [mm/day]. 
      # It MUST yield [mm/day]
      vol.qgev.fut.bc <- qgevd(vec.pgev.fut.bc,
                               location = par.bc.gcm.baseline[1],
                               scale = par.bc.gcm.baseline[2],
                               shape = par.bc.gcm.baseline[3]) 
      
      # Temporal intensity is requested
      #print(vol.qgev.fut.bc)
      
      # GCM-RCM Model Delta-Ratio (Future/Historical) unitless [-] for 24 hr. Unitless
      delta.ratio.fut.bc <- rcp.target.day.bc[j]/vol.qgev.fut.bc
      
      # Temporal Future/Historical is requested
      #print(delta.ratio.fut.bc)
      
      #-------------------------------------------------------------------------------------------------------------------
      # SUB-BLOCK: Spatial BC for GCM-RCM vs Observations
      #-------------------------------------------------------------------------------------------------------------------
      
      # GCM-RCM Model Future with respect to Historical GCM-RCM [24 hr]. 
      # It MUST yield: probability [-]
      vec.pgev.fut.bc.02 <- pgevd(vol.qgev.fut.bc,
                                  location = par.bc.gcm.baseline[1],
                                  scale = par.bc.gcm.baseline[2],
                                  shape = par.bc.gcm.baseline[3]) 
      
      # Temporal probability is requested
      #print(vec.pgev.fut.bc.02)
      
      # A duration vector is created [min]
      cont.sel.int <- c(5, 10, 15, 30, 60, 120, 360, 720, 1440)
      
      # List containers are created
      cont.vol.cont.bc <- NULL
      cont.int.cont.bc <- NULL
      
      # List containers are defined
      cont.vol.cont.bc <- vector()
      cont.int.cont.bc <- vector()
      
      # ============================
      # Inner loop is initialized
      # ============================
      for(i in 1:9) { # length = 9: GEV parameters, from 5-min to 1440-min per BC-method
        
        # vector data is erased
        par.obs.hist.min <- NULL
        
        # AMP geometric parameters are selected
        par.obs.hist.min <- as.numeric(c(df.para.vol.gev.iso[2, i],
                                         df.para.vol.gev.iso[3, i],
                                         df.para.vol.gev.iso[4, i]))
        
        # AMP geometric parameters are requested
        #print(par.obs.hist.min)
        
        # GCM-RCM Model Future reversed with respect to Historical Observations e.g. [mm/5-min]
        # in order to obtain a Model Quasi-Historical for e.g. 5 minutes ONLY. 
        # Therefore, it MUST yield [mm/5-min]
        vol.qgev.fut.bc.02 <- qgevd(vec.pgev.fut.bc.02,
                                    location = par.obs.hist.min[1],
                                    scale = par.obs.hist.min[2],
                                    shape = par.obs.hist.min[3])
        
        # Temporal volume is requested
        #print(vol.qgev.fut.bc.02)
        
        #-------------------------------------------------------------------------------------------------------------------
        # SUB-BLOCK: Temporal BC for Future GCM-RCM ONLY
        #-------------------------------------------------------------------------------------------------------------------
        
        # Future Model Quasi-Historical is projected by applying the Delta-Ratio
        # This is for future e.g. [mm/5-min] for RCP target[j] ONLY 
        # It MUST yield: e.g. [mm/5-min] for year target[j] 
        cont.vol.cont.bc[i] <- vol.qgev.fut.bc.02*delta.ratio.fut.bc 
        
        # Intensities are calculated [mm/hr]
        cont.int.cont.bc[i] <- 60*(vol.qgev.fut.bc.02*delta.ratio.fut.bc)/(cont.sel.int[i])
        
        # Temporal intensity is requested
        #print(cont.int.cont.bc[i])
        
        # ============================  
      } # End of inner-loop
        # ============================
      
      # List containers are filled
      list.vol.cont.bc[[j]] <- cont.vol.cont.bc
      list.int.cont.bc[[j]] <- cont.int.cont.bc
      
      #------------------------------------------------------------------------
      # Current loop is printed
      # Beware !!! this print message is computationally expensive
      #print(paste0("Target-day/ Year: ", (j+2019),
      #             " / BC-Method: ", bc.selector,
      #             " / PDF-selector: ", pfd.selector))
      #------------------------------------------------------------------------
      
      # ============================
    } # End of middle-loop
      # ============================
    
    # An external list of data.frames is created
    df.vol.cont.bc <- as.data.frame(do.call(rbind, list.vol.cont.bc))
    df.int.cont.bc <- as.data.frame(do.call(rbind, list.int.cont.bc))
    
    # Year variable is added to data.frames
    df.vol.cont.bc$date <- bc_sim_chosen_fut[, 1]
    df.int.cont.bc$date <- bc_sim_chosen_fut[, 1]
    
    # Names are assigned to matrix columns
    names(df.vol.cont.bc) <- c("min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720", "min_1440", "year")
    names(df.int.cont.bc) <- c("min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720", "min_1440", "year")
    
    # data.frame variables are reorganized
    df.vol.cont.bc <- df.vol.cont.bc %>% relocate("year", .before = "min_5")
    df.int.cont.bc <- df.int.cont.bc %>% relocate("year", .before = "min_5")
    
    # data.frames are transformed in to matrix objects
    matrix.base.int.sim.bc <- as.matrix(df.int.cont.bc)
    matrix.base.vol.sim.bc <- as.matrix(df.vol.cont.bc)
    
    # Names are assigned to matrix columns
    colnames(matrix.base.int.sim.bc) <- c("year", "5",	"10",	"15",	"30", "60", "120", "360", "720", "1440")
    colnames(matrix.base.vol.sim.bc) <- c("year", "5",	"10",	"15",	"30", "60", "120", "360", "720", "1440")
    
    #-------------------------------------------------------------------------------------------------------------------
    # SUB-BLOCK: Sub-Daily Annual Maximum Precipitation (AMP) 
    #-------------------------------------------------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------------------------------
    # GENERALIZED EXTREME VALUE (GEV) DISTRIBUTION [GEV]
    # -------------------------------------------------------------------------------------------------
    # Generic                LOCATION   SCALE     SHAPE
    # IDFtool                xi         alpha     kappa
    # -------------------------------------------------------------------------------------------------

    # IDFCurve {IDFtool} function is called and an new object is created
    IDF.Auto.int.gev.sim.bc <- IDFCurve(Data = matrix.base.int.sim.bc, Station = sta.id,
                                        Duration = v.duration, Periods = v.periods, 
                                        Type = "gev", M.fit = pdf.fit, 
                                        Plot = FALSE, Strategy = 1, logaxe = "", 
                                        CI = FALSE, CIpdf = FALSE, iter = 100,
                                        goodtest = FALSE, Resolution = 600, 
                                        SAVE = FALSE, name = TRUE)

    #----------------------------------------------------------------------------
    # SUB-BLOCK: PDF/Model intensities selection
    #----------------------------------------------------------------------------
    if (pfd.selector == 1) {
      df.intensity.hist <- IDF.Auto.int.gev$Intensities
      df.intensity.sim.bc <- IDF.Auto.int.gev.sim.bc$Intensities
      pdf.id <- c("pdf-adjusted")
    } else {
      df.intensity.hist <- IDF.Auto.int.gev$Models$HIDFUN$Predict
      df.intensity.sim.bc <- IDF.Auto.int.gev.sim.bc$Models$HIDFUN$Predict
      pdf.id <- c("model-adjusted")
    }
    #----------------------------------------------------------------------------
    
    # IDF delta intensities are created
    df.ratio.intensity.bc <-  round(df.intensity.sim.bc/df.intensity.hist, 3)
    
    # Matrix objects are transformed to data.frame objects
    df.intensity.sim.bc <- as.data.frame(df.intensity.sim.bc)
    df.ratio.intensity.bc <- as.data.frame(df.ratio.intensity.bc)
    
    # ID variables are incorporated
    df.intensity.sim.bc$bc <- v.bc.arch[m]
    df.ratio.intensity.bc$bc <- v.bc.arch[m]
    
    # ID variables are incorporated
    df.intensity.sim.bc$rcp_id <- rcp.id
    df.intensity.sim.bc$rcm_id <- rcm.id
    df.intensity.sim.bc$sta_id <- sta.id
    df.intensity.sim.bc$adj_id <- pdf.id
    
    # ID variables are incorporated
    df.ratio.intensity.bc$rcp_id <- rcp.id
    df.ratio.intensity.bc$rcm_id <- rcm.id
    df.ratio.intensity.bc$sta_id <- sta.id
    df.ratio.intensity.bc$adj_id <- pdf.id
    
    # data.frame is converted from wide to long format for ratio ONLY
    df.compost.bc.ratio <- data.frame(newcol = c(t(df.ratio.intensity.bc[, 1:length(v.periods)])), stringsAsFactors=FALSE)
    
    # Repetitions from 5-min to 1440-min per BC-method are created
    df.compost.bc.ratio$return <- rep(v.periods, 9)
    
    # Repetitions are incorporated
    df.compost.bc.ratio$duration <- c(rep((v.duration[1]), length(v.periods)),
                                      rep((v.duration[2]), length(v.periods)),
                                      rep((v.duration[3]), length(v.periods)),
                                      rep((v.duration[4]), length(v.periods)),
                                      rep((v.duration[5]), length(v.periods)),
                                      rep((v.duration[6]), length(v.periods)),
                                      rep((v.duration[7]), length(v.periods)),
                                      rep((v.duration[8]), length(v.periods)),
                                      rep((v.duration[9]), length(v.periods)))
    
    # data.frame variables are renamed
    colnames(df.compost.bc.ratio) <- c("ratio", "period", "duration")
    
    # ID variables are incorporated
    df.compost.bc.ratio$bc <- v.bc.arch[m]
    df.compost.bc.ratio$rcp_id <- rcp.id
    df.compost.bc.ratio$rcm_id <- rcm.id
    df.compost.bc.ratio$sta_id <- sta.id
    df.compost.bc.ratio$adj_id <- pdf.id
    
    # data.frame is converted form wide to long format for intensity ONLY
    df.compost.bc.int <- data.frame(newcol = c(t(df.intensity.sim.bc[, 1:length(v.periods)])), stringsAsFactors=FALSE)
    
    # Repetitions from 5-min to 1440-min per BC-method are created
    df.compost.bc.int$return <- rep(v.periods, 9)
    
    # Repetitions are incorporated
    df.compost.bc.int$duration <- c(rep((v.duration[1]), length(v.periods)),
                                    rep((v.duration[2]), length(v.periods)),
                                    rep((v.duration[3]), length(v.periods)),
                                    rep((v.duration[4]), length(v.periods)),
                                    rep((v.duration[5]), length(v.periods)),
                                    rep((v.duration[6]), length(v.periods)),
                                    rep((v.duration[7]), length(v.periods)),
                                    rep((v.duration[8]), length(v.periods)),
                                    rep((v.duration[9]), length(v.periods)))
    
    # data.frame variables are renamed
    colnames(df.compost.bc.int) <- c("ratio", "period", "duration")
    
    # ID variables are incorporated
    df.compost.bc.int$bc <- v.bc.arch[m]
    df.compost.bc.int$rcp_id <- rcp.id
    df.compost.bc.int$rcm_id <- rcm.id
    df.compost.bc.int$sta_id <- sta.id
    df.compost.bc.int$adj_id <- pdf.id
    
    # List containers are filled
    list.idf.bc.ratio[[m]] <- df.compost.bc.ratio
    list.idf.bc.int[[m]] <- df.compost.bc.int
    
    #----------------------------------------------------------------------------
    # SUB-BLOCK: PDF/Model parameter extraction xi, alpha, kappa per BC
    #----------------------------------------------------------------------------
    
    # bc.model.selector is evaluated
    if (p == 1) {
      
      # Duration vector for xi, alpha, kappa parameters are created [-]
      v.5min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`5 min`$Parameters$para)
      v.10min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`10 min`$Parameters$para)
      v.15min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`15 min`$Parameters$para)
      v.30min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`30 min`$Parameters$para)
      v.60min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`60 min`$Parameters$para)
      v.120min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`120 min`$Parameters$para)
      v.360min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`360 min`$Parameters$para)
      v.720min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`720 min`$Parameters$para)
      v.1440min.gev <- as.numeric (IDF.Auto.int.gev.sim.bc$Distribution$`1440 min`$Parameters$para)
      
      # A rbind matrix object is created
      df.tot.min.gev <- rbind(v.5min.gev,
                              v.10min.gev,
                              v.15min.gev,
                              v.30min.gev,
                              v.60min.gev,
                              v.120min.gev,
                              v.360min.gev,
                              v.720min.gev,
                              v.1440min.gev)
      
      # Matrix object is transformed to data.frame
      df.tot.min.gev <- as.data.frame(df.tot.min.gev)
      
      # data.frame variables are renamed
      names(df.tot.min.gev) <- c("xi", "alpha", "kappa")
      
      # ID variables are incorporated
      df.tot.min.gev$bc <- bc.selector
      df.tot.min.gev$pdf <- as.vector(p)
      df.tot.min.gev$rcp_id <- rcp.id
      df.tot.min.gev$rcm_id <- rcm.id
      df.tot.min.gev$sta_id <- sta.id
      df.tot.min.gev$dur <- c("min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720", "min_1440")
      
      # List containers are filled
      list.idf.dur.par[[m]] <- df.tot.min.gev
      
      # ==============================
    } # End of middle-loop
      # ==============================
    
    #----------------------------------------------------------------------------
    # SUB-BLOCK: PDF/Model parameter extraction for selected production BC
    #----------------------------------------------------------------------------
    
    # bc.model.selector is evaluated
    if (m == bc.model.selector) {
      
      # IDFCurve {IDFtool} function is called and an new object is created
      IDF.Auto.gev.select <- IDFCurve(Data = matrix.base.int.sim.bc, Station = sta.id,
                                      Duration = v.duration, Periods = v.periods, 
                                      Type = "gev", M.fit = pdf.fit, 
                                      Plot = FALSE, Strategy = 1, logaxe = "", 
                                      CI = FALSE, CIpdf = FALSE, iter = 100,
                                      goodtest = TRUE, Resolution = 600, 
                                      SAVE = FALSE, name = TRUE)
      
      # Goodness-of-fit tests are requested
      df.fit.int.gev.select <- as.data.frame(IDF.Auto.gev.select$Test.fit)
      
      # Rownames are added as variables
      df.fit.int.gev.select <- tibble::rownames_to_column(df.fit.int.gev.select, var = "dur")
      
      # A summary IDF coefficients data.frame is created for selected BC method
      # model: Inten ~ ((A)/((B + Dura)^C))
      df.para.int.gev.coeff.select <- as.data.frame(IDF.Auto.gev.select$Models$HIDFUN$Coefficients)
      
      # Rownames are added as variables
      # model: Inten ~ ((A)/((B + Dura)^C))
      df.para.int.gev.coeff.select <- tibble::rownames_to_column(df.para.int.gev.coeff.select, var = "period")
      
      # BR2 and RMSE metrics are isolated
      temp.metrics <- as.data.frame(IDF.Auto.gev.select$Models$HIDFUN$test.fit.reg)
      
      # BR2 and RMSE metrics are incorporated
      df.para.int.gev.coeff.select$BR2 <- temp.metrics$BR2
      df.para.int.gev.coeff.select$RMSE <- temp.metrics$RMSE

      # ID variables are incorporated
      df.para.int.gev.coeff.select$bc <- c("rqm")
      df.para.int.gev.coeff.select$rcp_id <- rcp.id
      df.para.int.gev.coeff.select$rcm_id <- rcm.id
      df.para.int.gev.coeff.select$sta_id <- sta.id
      df.para.int.gev.coeff.select$adj_id <- pdf.id
       
      # ==============================
    } # End of middle-loop
      # ==============================
    
    #------------------------------------------------------------------------
    # Current loop is printed
    # Beware !!! this print message is computationally inexpensive
    print(paste0("Target-day/ Year: ", (j+2019),
                 " / BC-Method: ", bc.selector,
                 " / PDF-selector: ", pfd.selector))
    #------------------------------------------------------------------------
  
    # ============================
  } # End of outer-loop
    # ============================
  
  # An external list of data.frames is created
  df.outer.ratio.idf <- as.data.frame(do.call(rbind, list.idf.bc.ratio))
  df.outer.int.idf <- as.data.frame(do.call(rbind, list.idf.bc.int))
  df.outer.dur.par <- as.data.frame(do.call(rbind, list.idf.dur.par))
  
  # List containers are filled
  list.outer.ratio.idf[[p]] <- df.outer.ratio.idf
  list.outer.int.idf[[p]] <- df.outer.int.idf
  list.outer.dur.par[[p]] <- df.outer.dur.par
  
  # ============================
} # End of outermost-loop
  # ============================





# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$
# &^%&%&%&%&^%$^%^&&^%&%&^%&$%&%^%!#%%&&%^%**$








# A combined external list of data.frames is created
df.outermost.ratio.idf <- as.data.frame(do.call(rbind, list.outer.ratio.idf))
df.outermost.int.idf <- as.data.frame(do.call(rbind, list.outer.int.idf))
df.outermost.dur.par <- as.data.frame(do.call(rbind, list.outer.dur.par))

# A pdf-adjusted only subset is created
df.outermost.ratio.idf.pdf.adjusted <- subset(df.outermost.ratio.idf, adj_id == "pdf-adjusted")
df.outermost.int.idf.pdf.adjusted <- subset(df.outermost.int.idf, adj_id == "pdf-adjusted")

# A model-adjusted only subset is created
df.outermost.ratio.idf.model.adjusted <- subset(df.outermost.ratio.idf, adj_id == "model-adjusted")
df.outermost.int.idf.model.adjusted <- subset(df.outermost.int.idf, adj_id == "model-adjusted")

# pdf-adjusted and model-adjusted historical observed intensities are isolated
df.intensity.hist.probe.pdf.adjusted <- as.data.frame(IDF.Auto.int.gev$Intensities)
df.intensity.hist.probe.model.adjusted <- as.data.frame(IDF.Auto.int.gev$Models$HIDFUN$Predict)

# ID variables are incorporated
df.intensity.hist.probe.pdf.adjusted$bc <- c("hist")
df.intensity.hist.probe.pdf.adjusted$rcp_id <- rcp.id
df.intensity.hist.probe.pdf.adjusted$rcm_id <- rcm.id
df.intensity.hist.probe.pdf.adjusted$sta_id <- sta.id
df.intensity.hist.probe.pdf.adjusted$adj_id <- c("pdf-adjusted")

df.intensity.hist.probe.model.adjusted$bc <- c("hist")
df.intensity.hist.probe.model.adjusted$rcp_id <- rcp.id
df.intensity.hist.probe.model.adjusted$rcm_id <- rcm.id
df.intensity.hist.probe.model.adjusted$sta_id <- sta.id
df.intensity.hist.probe.model.adjusted$adj_id <- c("model-adjusted")

# data.frame is converted form wide to long format for intensity ONLY
df.temp.pdf.adjusted <- data.frame(newcol = c(t(df.intensity.hist.probe.pdf.adjusted[, 1:length(v.periods)])), stringsAsFactors=FALSE)

# Repetitions from 5-min to 1440-min per BC-method are created
df.temp.pdf.adjusted$return <- rep(v.periods, 9)

# Repetitions are incorporated
df.temp.pdf.adjusted$duration <- c(rep((v.duration[1]), length(v.periods)),
                                   rep((v.duration[2]), length(v.periods)),
                                   rep((v.duration[3]), length(v.periods)),
                                   rep((v.duration[4]), length(v.periods)),
                                   rep((v.duration[5]), length(v.periods)),
                                   rep((v.duration[6]), length(v.periods)),
                                   rep((v.duration[7]), length(v.periods)),
                                   rep((v.duration[8]), length(v.periods)),
                                   rep((v.duration[9]), length(v.periods)))

# data.frame variables are renamed
colnames(df.temp.pdf.adjusted) <- c("ratio", "period", "duration")

# ID variables are incorporated
df.temp.pdf.adjusted$bc <- c("hist")
df.temp.pdf.adjusted$rcp_id <- rcp.id
df.temp.pdf.adjusted$rcm_id <- rcm.id
df.temp.pdf.adjusted$sta_id <- sta.id
df.temp.pdf.adjusted$adj_id <- c("pdf-adjusted")

# data.frame is converted form wide to long format for intensity ONLY
df.temp.model.adjusted <- data.frame(newcol = c(t(df.intensity.hist.probe.model.adjusted[, 1:length(v.periods)])), stringsAsFactors=FALSE)

# Repetitions from 5-min to 1440-min per BC-method are created
df.temp.model.adjusted$return <- rep(v.periods, 9)

# Repetitions are incorporated
df.temp.model.adjusted$duration <- c(rep((v.duration[1]), length(v.periods)),
                                     rep((v.duration[2]), length(v.periods)),
                                     rep((v.duration[3]), length(v.periods)),
                                     rep((v.duration[4]), length(v.periods)),
                                     rep((v.duration[5]), length(v.periods)),
                                     rep((v.duration[6]), length(v.periods)),
                                     rep((v.duration[7]), length(v.periods)),
                                     rep((v.duration[8]), length(v.periods)),
                                     rep((v.duration[9]), length(v.periods)))

# data.frame variables are renamed
colnames(df.temp.model.adjusted) <- c("ratio", "period", "duration")

# ID variables are incorporated
df.temp.model.adjusted$bc <- c("hist")
df.temp.model.adjusted$rcp_id <- rcp.id
df.temp.model.adjusted$rcm_id <- rcm.id
df.temp.model.adjusted$sta_id <- sta.id
df.temp.model.adjusted$adj_id <- c("pdf-adjusted")

# Relevant data.frames are rbinded
df.outermost.int.idf.pdf.adjusted <- rbind(df.outermost.int.idf.pdf.adjusted, df.temp.pdf.adjusted)
df.outermost.int.idf.model.adjusted <- rbind(df.outermost.int.idf.model.adjusted, df.temp.model.adjusted)

# data.frame variables kappa, bc and rcm_id are isolated 
df.radar.dur.par <- df.outermost.dur.par[, c(3,4,9)]

# Kappa is rounded to 3 decimals
df.radar.dur.par$kappa <- round(df.radar.dur.par$kappa, 3)

# data.frame is transformed from long to wide format
df.wide.radar.dur.par <- reshape2::dcast(data = df.radar.dur.par,
                                         formula = bc ~ dur,
                                         fun.aggregate = sum,
                                         value.var = "kappa")

# data.frame variables order is defined
col_order <- c("bc","min_5",	"min_10",	"min_15",	"min_30", "min_60", "min_120", "min_360", "min_720", "min_1440")

# data.frame variables are ordered
df.wide.radar.dur.par <- df.wide.radar.dur.par[, col_order]

# BC methods are coverted to uppercase
df.wide.radar.dur.par$bc <- toupper(df.wide.radar.dur.par$bc)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK:  ECDF + graphical analysis
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A new scale color is defined
scale01 <- c("#a6761d","#ca0020","#4dac26", "#404040", "#f4a582", "#0571b0")  
scale01 <- c ("#377eb8","#e41a1c", "#4daf4a", "black", "#ff7f00", "blue", "orange")
scale01 <- c ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02" , "#666666")
scale01 <- c ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
scale02 <- c ("#1B9E77", "#D95F02", "#7570B3", "#666666", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

# Ratios pdf-adjusted
ggplot() +
  geom_point(aes(x = duration,y = ratio,shape = bc,colour = bc),df.outermost.ratio.idf.pdf.adjusted, size = 2.00) +
  facet_wrap(facets = ~period) +
  geom_line(aes(x = duration,y = ratio,colour = bc,linetype = "solid"),df.outermost.ratio.idf.pdf.adjusted, linewidth = 0.50) +
  scale_color_manual(values = scale01) +
  scale_shape_manual(values = c(2,13,4,19,6,7,8,9,10,11,12,13)) +
  scale_x_continuous(trans='log10', breaks=v.duration) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  geom_hline(yintercept = 1.0,colour = 'black',size = 1.25,linetype = 2) +
  ggtitle(paste(sta.id, rcm.id, "Future IDF ratios vs. observations per return period", unique(df.outermost.int.idf.pdf.adjusted$adj_id))) +
  ylab("Projected/Observed Intensity Ratio (-)") +
  xlab("Duration (min)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Ratios model-adjusted
ggplot() +
  geom_point(aes(x = duration,y = ratio,shape = bc,colour = bc),df.outermost.ratio.idf.model.adjusted, size = 2.00) +
  facet_wrap(facets = ~period) +
  geom_line(aes(x = duration,y = ratio,colour = bc,linetype = "solid"),df.outermost.ratio.idf.model.adjusted, linewidth = 0.50) +
  scale_color_manual(values = scale01) +
  scale_shape_manual(values = c(2,13,4,19,6,7,8,9,10,11,12,13)) +
  scale_x_continuous(trans='log10', breaks=v.duration) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  geom_hline(yintercept = 1.0,colour = 'black',size = 1.25,linetype = 2) +
  ggtitle(paste(sta.id, rcm.id, "Future IDF ratios vs. observations per return period", unique(df.outermost.ratio.idf.model.adjusted$adj_id))) +
  ylab("Projected/Observed Intensity Ratio (-)") +
  xlab("Duration (min)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Intensities pdf-adjusted
ggplot() +
  geom_point(aes(x = duration,y = ratio,shape = bc,colour = bc),df.outermost.int.idf.pdf.adjusted, size = 2.00) +
  facet_wrap(facets = ~ period, scales = 'free') +
  geom_line(aes(x = duration,y = ratio,colour = bc,linetype = "solid"),df.outermost.int.idf.pdf.adjusted, size = 0.50) +
  scale_color_manual(values = scale02) +
  scale_shape_manual(values = c(2,13,4,6,7,8,9,10,11,12,13,19)) +
  scale_x_continuous(trans='log10', breaks=v.duration) +
  #scale_y_continuous(trans='log10', breaks=v.duration) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, rcm.id, "Future IDF vs. observations per return period", unique(df.outermost.int.idf.pdf.adjusted$adj_id))) +
  ylab("Intensity (mm/hr)") +
  xlab("Duration (min)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Intensities model-adjusted
ggplot() +
  geom_point(aes(x = duration,y = ratio,shape = bc,colour = bc),df.outermost.int.idf.model.adjusted, size = 2.00) +
  facet_wrap(facets = ~ period, scales = 'free') +
  geom_line(aes(x = duration,y = ratio,colour = bc,linetype = "solid"),df.outermost.int.idf.model.adjusted, size = 0.50) +
  scale_color_manual(values = scale02) +
  scale_shape_manual(values = c(2,13,4,6,7,8,9,10,11,12,13,19)) +
  scale_x_continuous(trans='log10', breaks=v.duration) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, rcm.id, "Future IDF vs. observations per return period", unique(df.outermost.int.idf.model.adjusted$adj_id))) +
  ylab("Intensity (mm/hr)") +
  xlab("Duration (min)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# Intensities model-adjusted
ggplot() +
  geom_point(aes(x = duration,y = ratio,shape = bc,colour = bc),df.outermost.int.idf.model.adjusted, size = 2.00) +
  facet_wrap(facets = ~ period, scales = 'free') +
  geom_line(aes(x = duration,y = ratio,colour = bc,linetype = "solid"),df.outermost.int.idf.model.adjusted, size = 0.50) +
  scale_color_manual(values = scale02) +
  scale_shape_manual(values = c(2,13,4,6,7,8,9,10,11,12,13,19)) +
  scale_x_continuous(trans='log10', breaks=v.duration, limits = c(4,35)) +
  scale_y_continuous(trans='log10', breaks=v.duration, limits = c(50,400)) +
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle(paste(sta.id, rcm.id, "Future IDF vs. observations per return period", unique(df.outermost.int.idf.model.adjusted$adj_id))) +
  ylab("Intensity (mm/hr)") +
  xlab("Duration (min)") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# IDF historical intensities are requested
View(df.intensity.hist.probe.pdf.adjusted)
View(df.intensity.hist.probe.model.adjusted)

# Combined ratio data.frames are requested
View(df.outermost.ratio.idf)
View(df.outermost.ratio.idf.pdf.adjusted)
View(df.outermost.ratio.idf.model.adjusted)

# Combined intensity data.frames are requested
View(df.outermost.int.idf)
View(df.outermost.int.idf.pdf.adjusted)
View(df.outermost.int.idf.model.adjusted)

# Final selected BC method geometric parameters is requested
# are requested for GEV ONLY
View(df.para.int.gev.coeff.select)

# Parameters xi, alpha, kappa pdf.adjusted is requested
View(df.outermost.dur.par)
View(df.wide.radar.dur.par)

# An output file name is concatenated
t.output01 <- paste("BRICK07_matrix_bc_AMP_future_par_",rcm.id,"_",rcp.id,"_",sta.id,".RData", sep = "")
t.output02 <- paste("BRICK07_matrix_bc_AMP_future_par_",rcm.id,"_",rcp.id,"_",sta.id,".xlsx", sep = "")

# An export/import data.frame is created to free RAM memory
save(df.intensity.hist.probe.pdf.adjusted,      #1
     df.intensity.hist.probe.model.adjusted,    #2
     df.outermost.ratio.idf,                    #3
     df.outermost.ratio.idf.pdf.adjusted,       #4
     df.outermost.ratio.idf.model.adjusted,     #5
     df.outermost.int.idf,                      #6
     df.outermost.int.idf.pdf.adjusted,         #7
     df.outermost.int.idf.model.adjusted,       #8
     df.para.int.gev.coeff.select,              #9
     df.outermost.dur.par,                      #10
     df.wide.radar.dur.par,                     #11
     file = t.output01)

# An export/import XLS object is created
xls.object <- list("df.int.obs.pdf.adjusted" = df.intensity.hist.probe.pdf.adjusted,       #1
                   "df.int.obs.model.adjusted" = df.intensity.hist.probe.model.adjusted,   #2
                   "df.outer.ratio.idf" = df.outermost.ratio.idf,                          #3
                   "df.outer.ratio.idf.pdf.adj" = df.outermost.ratio.idf.pdf.adjusted,     #4
                   "df.outer.ratio.idf.model.adj" = df.outermost.ratio.idf.model.adjusted, #5
                   "df.outer.int.idf" = df.outermost.int.idf,                              #6
                   "df.outer.int.idf.pdf.adj" = df.outermost.int.idf.pdf.adjusted,         #7
                   "df.outer.int.idf.model.adj" = df.outermost.int.idf.model.adjusted,     #8
                   "df.para.int.gev.coeff.select" = df.para.int.gev.coeff.select,          #9
                   "df.outermost.dur.par" = df.outermost.dur.par,                          #10
                   "df.wide.radar.dur.par"= df.wide.radar.dur.par)                         #11
write.xlsx(xls.object, file = t.output02, rowNames = TRUE)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
