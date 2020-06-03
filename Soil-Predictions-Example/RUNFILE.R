#Title: RUNFILE-Simple Soil-Predictions-Example
#Authors: Charlotte Rivard & Shree Dangal
#Date: 5/8/20
#Summary: #The following script predicts values for a single soil property (Organic Carbon)

#----------------------------------------------#
# Install Packages
#----------------------------------------------#
#install.packages(stringr) # preprocessing
#install.packages(foreach) # preprocessing
#install.packages(pls) # pls models
#install.packages(resemble) # mbl models


#----------------------------------------------#
# Reference Set Preprocessing #
#----------------------------------------------#
source("Functions/preprocess_functions.R")

# Get Spectral Library Set
all_data <- getSpecLib(SAVENAME="ALL_DATA")
View(all_data)

# Split into Reference Set and Prediction Set
split <- calValSplit(SPECLIB=all_data) 
refSet <- split$calib # Reference Set
predSet <- split$valid # Prediction Set

# Refine Spectral Library
refSet <- refineSpecLib(SPECLIB=refSet, PROP="OC", SAVENAME="refSet.OC")


#----------------------------------------------#
# Partial Least Squares Regression #
#----------------------------------------------#
source("Functions/runfile_functions.R")

# Make or Load Model
plsr.model <- makePLSModel(PROP="OC", REFSET=refSet)
#plsr.model <- makePLSModel(PROP="OC", REFPATH="Data_Processed/refSet.OC.RData")

# Make Predictions
pls.predictions <- getModResults(PROP="OC", MODTYPE="PLS", MODNAME= "plsr.model", PREDNAME= "predSet")
#pls.predictions <- getModResults(PROP="OC", MODTYPE="PLS", MODNAME="plsr.model", PREDPATH="./Data_Processed/valid.ALL.RData")


#----------------------------------------------#
# Memory Based Learner Model #
#----------------------------------------------#
source("Functions/runfile_functions.R")

# Make or Load Model
mbl.sqrt <- runMBL(PROP="OC", REFNAME="refSet", PREDNAME="predSet")
#mbl.sqrt <- runMBL(PROP="OC", REFPATH="Data_Processed/refSet.OC.RData", PREDPATH="Data_Processed/valid.ALL.RData")

# Extract Predictions
mbl.predictions <- getModResults(PROP="OC", MODTYPE="MBL", MODNAME= "mbl.sqrt", PREDNAME= "predSet")
#mbl.predictions <- getModResults(PROP="OC", MODTYPE="MBL", MODNAME="mbl.sqrt", PREDPATH="Data_Processed/valid.ALL.RData")



