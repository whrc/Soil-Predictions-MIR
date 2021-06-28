**Scripts**\
\
_setname_prep.R_\
Performs the calibration transfer on the spectra and saves as RData file in ‘spc’ folder\
Change the input csv file, the columns being selected as spectra (lines 11-12), and output name/location\
\
_setname_oc.R_\
Submit as a job through cloudops, creates all the mbl models with different parameter combinations, to output/oc folder\
Change input validation and calibration sets (line 32-38), property (oc) throughout the file, output location (line 107) and create output folder for soil property\
\
_setname-fratio.R_\
Calculates the fratio for all samples in the calibration and validation sets and outputs a list of outlier indices from the combined dataset.\
Ex: calibration set indices are 1-15000, validation set indices are from 15001-15240\
Change input calibration and validation spectra (5-6 and throughout), number of directories in line 8 as needed, output location- currently ‘fratio’ subfolder.\
\
_setname-extract.R_\
Creates comprehensive files containing all predictions for each mbl model by property. (ie. pred.oc.csv, pred.bd.csv)\
Creates a file containing the lower, mean and upper prediction estimates for each property across all models (all-predictions.csv)\

