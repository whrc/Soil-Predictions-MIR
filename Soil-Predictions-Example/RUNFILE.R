#Title: RUNFILE-Simple Soil-Predictions-Example
#Authors: Charlotte Rivard & Shree Dangal
#Date: 5/8/20
#Summary: #The following script predicts values for a single soil property (Organic Carbon)

#----------------------------------------------#
# Install Packages
#----------------------------------------------#
#install.packages(stringr)      # processing spectra
#install.packages(foreach)      # processing spectra
#install.packages(prospectr)    # processing spectral set
#install.packages(clhs)         # processing large sets
#install.packages(matrixStats)  # preprocessing baseline transformation
#install.packages(pls)          # pls models
#install.packages(resemble)     # mbl models


#----------------------------------------------#
# Data Preprocessing #
#----------------------------------------------#
source("Functions/preprocess_functions.R")

# Get Spectral Library Set
ALL_data <- getSpecLib(SAVENAME="ALL_data")

# Refine Spectral Library
OC_data <- refineSpecLib(SPECLIB=ALL_data, PROP="OC", CALVAL=TRUE, SAVENAME="OC_data")

# Define Reference and Prediction Sets
refSet <- OC_data[OC_data$calib==1,]
predSet <- OC_data[OC_data$calib==0,]

#----------------------------------------------#
# Partial Least Squares Regression #
#----------------------------------------------#
source("Functions/plsr_functions.R")
source("Functions/perform_functions.R")

# Make Model
plsr.OC <- makePLSModel(PROP="OC", REFNAME="refSet")

# Make Predictions
pls.predictions <- getModResults(PROP="OC", MODTYPE="PLS", MODNAME= "plsr.OC", PREDNAME= "predSet")

#----------------------------------------------#
# Memory Based Learner Model #
#----------------------------------------------#
source("Functions/mbl_functions.R")
source("Functions/perform_functions.R")

# Make Model
mbl.OC <- runMBL(PROP="OC", REFNAME="refSet", PREDNAME="predSet")

# Extract Predictions
mbl.predictions <- getModResults(PROP="OC", MODTYPE="MBL", MODNAME= "mbl.OC", PREDNAME= "predSet")


