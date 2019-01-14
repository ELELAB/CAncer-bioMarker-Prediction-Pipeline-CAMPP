# The-CAncer-bioMarker-Prediction-Pipeline-CAMPP-

'''Created by Thilde Bagger Terkelsen, 07-01-2018, thilde@cancer.dk, thildebate@gmail.com"'''

This GitHub repository was created for the The CAncer bioMarker Prediction Pipeline (CAMPP). 
A collection of R scripts made for computational prediction of serum cancer biomarkers. 
The pipeline was created by Thilde Bagger Terkelsen for internal use at the Danish Cancer Society Research Center. 
The pipeline consists of three R scrips which are currently run from the linux command-line using flags. The scrips are named:
                                  
                                  CAMPPInstall.R
                                  CAMPPFunctions.R
                                  CAMPP.R
                                  
The CAMPPInstall.R script is run the first time the pipeline is used to ensure that all needed R-packages are installed. 
This part may require the user to open R and manually pick a CRAN-mirror. The user guide explains this in detail.

The CAMPPFunctions.R script contains custum functions used for analysis and is sourced within the CAMPP.R. The CAMPPFunctions.R script should be located in the folder from which the pipeline is run.

The CAMPP.R script acts as the actual pipeline and is run with flags from the linux command-line. 
The flag -h is implemented for user help.
