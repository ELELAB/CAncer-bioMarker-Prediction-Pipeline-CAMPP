# The-CAncer-bioMarker-Prediction-Pipeline-CAMPP-
# Created by Thilde Bagger Terkelsen, 04-09-2017, thilde@cancer.dk, thildebate@gmail.com"

This GitHub repository was created for the The CAncer bioMarker Prediction Pipeline (CAMPP). 
A collection of R scripts made for computational prediction of serum cancer biomarkers. 
The pipeline was created by Thilde Bagger Terkelsen for internal use at the Danish Cancer Society Research Center. 
The pipeline consists of three R scrips which are currently run from the linux commandline using flags. The scrips are named:
                                  
                                  PipelineInstall.R
                                  PipelineFunctions.R
                                  Pipeline.R
                                  
The PipelineInstall.R script is run the first time the pipeline is used to ensure that all needed R-packages are install. 
This part may require the user to open R and manually pick a crane-mirror. The user guide explains this in detail.

The PipelineFunctions.R script contains custum functions used for analysis and is sourced directly from the Pipeline.R, 
- this script does not need to be run.

The Pipeline.R script acts as the actual pipeline and is run with user flags from the linux commandline. 
The flag -h is implemented for help.
