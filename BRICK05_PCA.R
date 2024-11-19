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
# INFO: This script is intended for 
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

# Workspace is cleared
#gc()
#rm(list = ls())

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/HadGEM2_RegCM47_AIJS")

# Scientific notation is disabled
options(scipen=999)

# Start time is recorded
start.time <- Sys.time()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(devtools)
require(investr)
require(DescTools)
require(dplyr)
require(ggplot2)
require(matrixStats)
require(tidyr)
require(openxlsx)
require(factoextra)
require(FactoMineR)
require(tidyverse)
require(corrplot)
require(PerformanceAnalytics)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CORDEX + Weather Station Metadata + Input Data
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Weather station ID
#-------------------------------------------------------------------------------------------------------------

# Weather station ID is defined
sta.id <- "AIJS"

# RCM ID is defined
rcm.id <- "HadGEM2_RegCM47"

# RCM ID is defined
rcp.id <- "RCP85"

# Weather station coordinates are defined
sta.long <- -84.200 # Aeropuerto Internacional Juan Santamaria
sta.lat <- 9.997    # Aeropuerto Internacional Juan Santamaria

#-------------------------------------------------------------------------------------------------------------
# SUB-BLOCK: Preprocessed binary datasets are loaded
#-------------------------------------------------------------------------------------------------------------

# RData objects are loaded from BRICK03
load("BRICK03_matrix_bc_24h_par_HadGEM2_RegCM47_RCP85_AIJS.RData")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: PCA Analysis
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------
# References
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/
# https://www.r-bloggers.com/2021/05/principal-component-analysis-pca-in-r/
# http://www.sthda.com/english/wiki/wiki.php?id_contents=7851
# https://www.kaggle.com/code/agailloty/comprehensive-pca-with-r-using-factominer
# http://factominer.free.fr/factomethods/principal-components-analysis.html
# https://rpubs.com/gokul108/pca1
#-------------------------------------------------------------------------------------------------------------

# A AMP24 metrics performance subset is executed
subset01 <- df.perf.desc.24hrAMP[, c(1:5)]

# A climdex metrics performance subset is executed
subset02 <- df.perf.climdex[c(3:7,1), (1:12)]

# A cbind data.frame is created
subset03 <- cbind(subset02, subset01)

# CORDEX member raw simulations are excluded
subset03 <- subset03[c(1:5),]

# Irrelevant variables are excluded
#subset03 <- subset03[, -c(1:3,8,9,12)] # BEWARE !!! this is further tests ONLY

# Relative weights are assigned to variables
#w.vector <- c(1,1,1,1,1,1,1,1,1,1,1) # BEWARE !!! this is further tests ONLY

# A correlation matrix is created
cor.mat <- round(cor(subset03), 2)

# A corrplot is created
corrplot(cor.mat,
         type="upper", 
         order="original", 
         tl.col="black", tl.srt=45,
         method="pie")

# A correlation chart is created
#chart.Correlation(subset03, histogram=TRUE, pch=19) # BEWARE !!! this is further tests ONLY

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: FactoMineR PCA
# ----------------------------------------------------------------------------------------------------------------------------------

# The goal of principal component analysis is to transform the initial variables 
# into a new set of variables which explain the variation in the data.

# PCA {FactoMineR} Principal Component Analysis function is requested
result.factor.NS <- PCA(X = subset03,
                        ncp = 5,
                        graph = TRUE,
                        col.w = NULL) # or NULL or w.vector

# Analysis:
# This plots shows the structural relationship between the variables and the components.
# The projection of a variable vector onto the component axis allows to directly
# read the correlation between the variable and the component.
# The axis that represents PC 1 and PC 2 is the Pearson coefficient of correlation
# which goes from - 1 to 1.
# First PC Read the correlation from left to right. The idea behind this plot 
# is to show in which direction the variables correlate.
# Second PC Read the correlation from the bottom to the up.

# eigenvalue {factoextra} is called to Extract and visualize the eigenvalues/variances of dimensions
fviz_eig(result.factor.NS,
         addlabels=TRUE,
         hjust = -0.3,
         barfill = "gray",
         barcolor = "black",
         linecolor = "#FC4E07",
         ncp = 5) +
  theme_bw(base_size = 12.0)

# A PCA plot is requested
plot(result.factor.NS, axes=c(1, 2), choix="var", col.var="blue",new.plot=TRUE)

# fviz_pca {factoextra} function is requested to visualize PCA
fviz_pca_var(result.factor.NS, repel = TRUE, labelsize = 6,
             col.var = "black", col.circle = "black",) +
  theme_bw(base_size = 22.0) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 90.0),
        panel.grid.major = element_line(colour = '#cccccc'))

# fviz_pca {factoextra} function is requested to visualize PCA
fviz_pca_var(result.factor.NS, col.var="contrib", # "contrib" or "cos"
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # c("#FFCC00", "#CC9933", "#660033", "#330033"),
             ggtheme = theme_minimal())

# An object summary is requested
summary(result.factor.NS)

# Eigenvalue percentage of variance cumulative is requested
result.factor.NS$eig

# Analysis:
# An eigenvalue > 1 indicates that PCs account for more variance 
# than accounted by one of the original variables in standardized data.
# Hence, the components with eigenvalue > 1 are retained. 

# A list of matrices containing all the results for the active variables is requested
result.factor.NS$var

# A list of matrices containing all the results for the contribution is requested
result.factor.NS$var$contrib

# Analysis:
# This data.frame gives the contribution of each variable in the building of the dimensions.
# The sum of all the contributions must be 100. Variables that are not correlated with any 
# PC or correlated with the last dimensions are variables with low contribution and might be 
# removed to simplify the overall analysis. The larger the value of the contribution, the more 
# the variable contributes to the component.

# A list of matrices containing all the results for cos2 is requested
result.factor.NS$var$cos2

# Analysis:
# The quality of representation of the variables on factor map is called cos2
# (square cosine, squared coordinates). A high cos2 indicates a good representation 
# of the variable on the principal component. In this case the variable is positioned 
# close to the circumference of the correlation circle. A low cos2 indicates that 
# the variable is not perfectly represented by the component.


# a list of matrices containing all the results for the active individuals is requested
result.factor.NS$ind

# dimdesc {FactoMineR} is used to identify the most significantly associated variables
# with a given principal component. The output will be sorted by p-values.
desc.var <- dimdesc(result.factor.NS, axes = c(1,2), proba = 0.05)
desc.var$Dim.1
desc.var$Dim.2

# ----------------------------------------------------------------------------------------------------------------------------------
# SUBBLOCK: FactoMineR Clustering
# ----------------------------------------------------------------------------------------------------------------------------------

# HCPC{FactoMineR} Hierarchical Clustering on Principle Components function is requested
cluster01 <-  HCPC(res = result.factor.NS,
                   nb.clust = -1, # -1 Auto
                   metric = "euclidean",
                   method = "ward",
                   graph = FALSE,
                   cluster.CA = "rows" , # "rows" or "columns"
                   iter.max = 10)

# fviz_dend {factoextra} is used to visualize dendrogram generated by hierarchical clustering
fviz_dend(cluster01, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8)      # Augment the room for labels

# Object cluster01 is plotted in 3D
plot(cluster01, choice = "3D.map", ind.names = TRUE)

# fviz_cluster {factoextra} is used to visualize individuals on the principal component map
fviz_cluster(cluster01,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE,  # Show cluster centers
             palette = "jco",         # Color palette see. ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map")

# An object summary is requested
summary(cluster01)

# data.clust: The original data with a supplementary column called class containing the partition.
# desc.var: The variables describing clusters
# desc.ind: The more typical individuals of each cluster
# desc.axes: The axes describing clusters

# BEWARE !!! there is missing info here as NULLs are returned !!!
cluster01$desc.var      # WHY ????
cluster01$desc.axes     # WHY ????
cluster01$desc.ind$para # WHY ????

# fviz_cluster {factoextra} is used to visualize individuals on the principal component map
g.cluster01 <- fviz_cluster(cluster01)

# A ggplot2 object is create
g.cluster01 +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw(base_size = 18.0) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 90.0),panel.grid.major = element_line(colour = '#cccccc'))

# fviz_pca {factoextra} is called to visualize Principal Component Analysis using biplot
fviz_pca_biplot(result.factor.NS, repel = TRUE, labelsize = 6, col.var = "black")

# fviz_pca {factoextra} is called to visualize Principal Component Analysis using biplot
fviz_pca_biplot(result.factor.NS,  geom = "text")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
