Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark.


This GitHub repository was created for the The CAncer bioMarker Prediction Pipeline (CAMPP). 
A collection of R scripts made for computational prediction of serum cancer biomarkers. 

The publication associated to CAMPP and to be cited in case of its usage is currently on biorxiv as a pre-print:

CAncer bioMarker Prediction Pipeline (CAMPP) - A standardised and user-friendly framework for the analysis of quantitative biological data. biorxiv, doi: https://doi.org/10.1101/608422

Thilde Terkelsen, Anders Krogh and Elena Papaleo

corresponding author: elenap@cancer.dk

contacts for software/scripts: thilde@cancer.dk, elenap@cancer.dk

The pipeline consists of two R scripts, which are currently run from the terminal command-line using flags. The scripts are named:
                                  
                                  CAMPPFunctions.R
                                  CAMPP.R
  
The CAMPPFunctions.R script contains custum functions used for analysis and is sourced within the CAMPP.R. The CAMPPFunctions.R script should be located in the folder from which the pipeline is run.

The CAMPP.R script acts as the actual pipeline and is run with flags from the command-line. 
The flag -h is implemented for user help.

The pipeline automatically checks for R-package dependencies and installs them if needed. However, as R-package updates may potentially break the code, we also provide the user with a "renv" library freeze (renv.lock). The user can specify running CAMPP with renv by setting the flag -e to "stable". See CAMPP Manual for specifics on this. 

Data example used in the CAMPPManual.pdf is located at https://github.com/ELELAB/N-glycan-TIF with the original publication, in the folder https://github.com/ELELAB/N-glycan-TIF/tree/master/Data/DataExamples.
