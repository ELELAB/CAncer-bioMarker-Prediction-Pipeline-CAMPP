#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                            ### USER ARGUMENTS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

install.packages("getopt", repos="https://cloud.r-project.org")
install.packages("BiocManager", repos="https://cloud.r-project.org")
library("getopt")
library("BiocManager")



# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "help", "h", 0, "logical",
  "rlib", "e", 2, "character",
  "variant", "v", 2, "character",
  "data", "d", 2, "character",
  "metadata", "m", 2, "character",
  "groups", "g", 2, "character",
  "datacheck", "j", 2, "logical",
  "kmeans", "k", 2, "character",
  "corr", "o", 2, "character",
  "survival", "u", 2, "character",
  "survplot", "q", 2, "numeric",
  "transform", "t", 2, "character",
  "standardize", "z", 2, "character",
  "batches", "b", 2, "character",
  "filename", "n", 2, "character",
  "sig", "f", 2, "character",
  "covar", "r", 2, "character",
  "stratify", "y", 2, "character",
  "plotmds", "s", 2, "logical",
  "colors", "c", 2, "character",
  "plotheatmap", "a", 2, "character",
  "lasso", "l", 2, "numeric",
  "WGCNA", "w", 2, "character",
  "cutoffWGCNA", "x", 2, "numeric",
  "PPint", "p", 2, "character",
  "GenemiRint", "i", 2, "character"), byrow=TRUE, ncol=4)


opt = getopt(spec)





# Help
if(!is.null(opt$help)) {
    cat("\nFlags:\n\n-e --rlib: Specify by if you want to re-install the newest versions of R-packages or if you want to run with the provided 'Renv' stable library. If the argument is omitted, R-packages will be installed (if nessesary) from an already specified CRAN mirror. The argument may be set to either a CRAN mirror of choice, see 'https://cran.r-project.org/mirrors.html'. If the argument it set to 'stable' (not case sensitive), the 'Renv' library freeze will be used instead.\n\n-d --data: file (xlsx or txt) with expression values, rows as features (genes, miRNAs, proteins, N-features) and columns as samples.\n\n-m --metadata: file (xlsx or txt) with metadata, minimum two columns one with ids (sample names matching those in the object above) and groups (diagnosis, tumor stage, ect.).\n(I) If the data comes from experimental batches and you want to correct for this, a column named 'batch' specifying which batch each sample belongs to (A,B,C,D, time1, time2, time3 ect) should also be included in the metadata. N.B specifying batches by numbers alone is not allowed.\n(II) If you are interested in performing survival analysis a column named 'survival' must be included specifying (in a binary way) which samples have survival information (1) and which do not (0). N.B. if you have paired cancer and normal samples the column 'survival' should only have the value 1/0 for tumour samples (NA or other character values should be used for normal samples.\n(IV) If you want to include covariates in your analysis these should be included in the metadata file as a column(s).\n\n-v --variant: Data 'variant'. Current options are 'array', 'seq', 'ms' or 'other'. This argument is mandatory and depeding on which option is chosen, data is transformed differently. If a second dataset is provided the -v option should be specified for each dataset, provided as a comma seperated list (no quotes, no paranthesis etc.).\n\n-g --groups: Argument -g should be specified as a comma separated list of length two (without quotes and parenthesis!). The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.\n\n-j Distributional Checks: Logical argument (TRUE/FALSE) which specifies whether Cullen-Frey graphs should be made for 10 randomly selected variables to check data distributions. This flag is per default set to TRUE.\n\n-k --kmeans: Argument for kmeans clustering. The -k flag must be specified as a character matching the name of a column in the metadata file, denoting the labeling of points on the MDS plot(s). If -k is empty (no column name specified) no labels will be added to the plot. \n\n-a --plotheatmap: Argument for heatmap specified as either; DE, DA, LASSO, EN or Consensus. \n\n-o --corr: Argument for correlation analysis. String specify which features should be correlated, options are: ALL, DE, DA, LASSO, EN or Consensus.\n\n-u --survival: Survival analysis may be performed on differentially expressed/abundant variables, variables from LASSO/EN regression or the consensus of these, e.g. flag -u must be specified as either; DE, DA, LASSO, EN or Consensus. In principal the full dataframe of variables may be used as well (if argument is set to ALL), HOWEVER this is not advisable unless the dataset is small with very few variables. Survival info must be included in the metadata excel file. The metadata file must contain at least four columns named; 'ids'(sample identifiers), 'age' (age in years at diagnosis, surgery or entry into trail), 'outcome.time' (time until end of follow-up in weeks, months or years, censuring, death) and 'outcome' (numeric 0 = censuring, 1=death). N.B. if you have (paired) normal samples the columns with survival information for these samples should contain NA values.\n\n-q --survplot: Arguments which specifies number of features to include per survival plot, e.g. many features requires splitting of the plot, default features per plot is 50.\n\n-z --standardize: Data centering. This option may be set to mean or median. If two datasets are provided the -z option should be specified for each dataset, provided as a comma seperated list (no quotes, no paranthesis etc.). If the flag -z is not specified and -v = array, then quantile normalization will be performed.\n\n-t --transform: should data be transformed? Current options are 'log2', 'log10' or 'logit'. If two datasets are provided the -t option should be specified for each dataset, provided as a comma seperated list (no quotes, no paranthesis etc.). If argument is left out, no transformation of data will occur.\n\n-b --databatch: Specifies if you want to correct for experimental sample (tissue/interstitial fluid) batches. Argument takes a string of length 1 (one dataset) or 2 (two datasets), where the string(s) match a column name(s) in the metadata file(s).\n\n-n --filename: Name of result files from analysis.\n\n-f --sig: Cut-offs for log fold change (logFC) and corrected p-value (fdr), defining significant hits (proteins, genes, miRNAs or N-features). If argument -f is set, it must be a comma separated list of length two (without quotes and parenthesis!), where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr). If omitted cutoffs will be set to -1 > logFC > 1 and corrected p-values < 0.05.\n\n-s --plotmds: TRUE or FALSE specifies if a preliminary MDSplot should be made for data overview.\n\n-r --covar: Covariates to include in analysis. If multiple of these, they should be specified with commas as separator, i.e. Covar1,Covar2,Covar3, (without quotes and parenthesis!). The first element in this list must be either TRUE or FALSE. If TRUE is specified then covariates will be included in both DE/DA analysis and Survival Analsysis. If FALSE is specified covariates will ONLY be used for Survival Analsysis. Names of covariates should match the desired columns in the metadata file.\n\n-y --stratify: This flag may be used if some of the categorical (NOT continious) covariates violate the cox proportional assumption. The pipline checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the pipeline with this flag followed by the names of the categorical covariates which failed and these will be stratified.\n\n-c --colors: Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) must be separted by commas, example: green,red,blue. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.\n\n-l --lasso: Argument specifying if LASSO or elestic net regression should be performed. This argument may be set to a value between 1 (LASSO) and > 0 (Elastic Net), but NOT to 0 exactly (Ridge Regression).\n\n-w --WGCNA: Argument specifying if Weighed Gene Co-expression Network Analysis should be performed. It takes a string, either DA, DE or ALL specifying if all variables should be included in WGCNA or only differentially expressed / abundant variables. \n\n-x --cutoffWGCNA: Argument specifying the cutoff value, in %, for top most interconnected genes (or other features) from each modules identified in the Weighed Gene Co-expression Network Analysis.\n\n-p --PPint: Argument specifying that protein-protein interaction networks should be generated using the results of the differential expression analysis. This argument must be a list of length two (without quotes and parenthesis!). The first element in this list must be a string specifying the type of gene identifier in the DE/DA results file provided, allowed identifiers are:\nensembl_peptide_id\nhgnc_symbol\nensembl_gene_id\nensembl_transcript_id\nuniprotswissprot\nThe second element in this list must be a string specifying the version of the stringDB to use. Currently only version supported is:\n11.0.\n\n-i--GenemiRint: Argument specifying that gene-miRNA interaction networks should be generated using the results of the differential expression analysis. This argument must be a list of length three (without quotes and parenthesis!). The first element in this list must be a string specifying the type of miRNA identifier in the expression data file, allowed identifiers are:\nmature_mirna_ids\nmature_mirna_accession.\nThe second element in this list must be a string specifying the miRNA-gene database to use, currently options are:\ntargetscan (validated miRNAs)\nmirtarbase (predicted miRNAs)\ntarscanbase (validated + predicted miRNAs)\n\n")
    
    stop("\n- Argument -h (help) selected, exiting script.")
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# START LOG FILE
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sink("./CAMPPlog.txt", append=TRUE, split=TRUE)
cat("\nCAMPP Running Messages:\n")


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                         ### INSTALL AND LOAD R-PACKAGES ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LIST OF R-PACKAGES
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


list.of.packages <- c("limma", "edgeR", "glmnet", "sva", "openxlsx", "ggplot2", "gridExtra", "heatmap.plus", "plyr", "data.table", "viridis", "squash", "survcomp", "survminer", "scales", "rms", "WGCNA", "fitdistrplus", "impute", "pcaMethods", "pROC", "VennDiagram", "mclust", "multiMiR", "biomaRt", "devtools", "arcdiagram", "doMC", "knor", "nnet")

# "stackoverflow"


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 1
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


source("CAMPPFunctions.R")


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NEWEST R-PACKAGES INSTALL AND LOAD OR PACKRAT LIBRARY
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(is.null(opt$rlib)) {
    arg.rlib <- "https://mirrors.dotsrc.org/cran/"
    is.rlib <- FALSE
} else if (length(grep("stable|Stable|STABLE",opt$rlib, ignore.case=TRUE,value=TRUE)) == 1) {
    is.rlib <- TRUE
} else {
    arg.rlib <- opt$rlib
    is.rlib <- FALSE
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (is.rlib == FALSE) {
    cat("\npackrat library not used. CAMPP will use available R-packages from user library and install these if needed.\n")
    
    missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    
    if (length(missing.packages) > 0) {
        cat("\nThe following package(s) is(are) not installed:\n")
        
        for (p in missing.packages) {
            print(p)
        }
        
        if ("arcdiagram" %in% missing.packages) {
            install_github("https://github.com/gastonstat/arcdiagram")
        }
        
        for (idx in 1:length(missing.packages)) {
            install.packages.auto(as.character(missing.packages[[idx]]), arg.rlib)
        }
    
    } else {
        cat("All packages are installed.\n")
    }
    
    
    file <- try(lapply(list.of.packages, library, character.only=T))
    
    if (class(file) == "try-error") {
        
        stop("Ups! Something went wrong, one or more R-package dependencies are not installed correctly. Check the script CAMPPmissingpackages.R. Alternatively you can download and use the Packrat library freeze, see manual for specifics. Packrat library in brief:\n\n (1.) Download the Packrat library from the CAMPP repository on github (if you have not already), and make sure it is located in the same folder as CAMPP.R.\n\n (2.) set the flat -e to TRUE, and run the pipeline.")
        
    } else {
        rm(file)
        cat("\nPACKAGES HAVE BEEN INSTALLED - READY TO RUN CAMPP!\n")
        cat("\n---------------------------------------------------------------------------------------------\n")
    }
    
} else if (is.rlib == TRUE) {
    
    if (!require("renv")) {
        install.packages("renv", repos="https://cloud.r-project.org")
    }
    library("renv")
    cat("\nStable 'Renv' library is being used for analysis!\n")
    renv::consent(provided = TRUE)
    renv::restore()
    file <- try(lapply(list.of.packages, library, character.only=T))
    cat("\nPACKAGES HAVE BEEN INSTALLED - READY TO RUN CAMPP!\n")
    cat("\n---------------------------------------------------------------------------------------------\n")
} else {
    cat("\nFlag -e must be either null (omitted) or TRUE or FALSE!.\n")
}









# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### ARGUMENTS SPECIFYING DATASETS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Vector to check user inputs
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# For DE/DA analysis, survival analysis and correlation analysis - user input must be in list:
must.be.type <- c("ALL","EN", "LASSO", "DA", "DE", "Consensus")

# For survival analysis, must be in metadata file:
must.contain <- c("survival", "outcome", "outcome.time")

# For miRNA-gene and protein-protein network analysis - user inputs must be in lists:
approvedGeneIDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
approvedmiRIDs <- c("mature_mirna_ids", "mature_mirna_accession")
genequery <- c("stringdatabase")
miRNAquery <- c("targetscan", "mirtarbase", "tarscanbase")





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Datasets
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (is.null(opt$data)) {
    stop("\n- Argument data (-d) is missing. Data provided must be an excel file with features (genes, proteins ect.) as rows and samples as columns.\n")
} else {
    arg.sets <- SplitList(opt$data)
    arg.data <- arg.sets[1]
    arg.data <- ReadMyFile(arg.data, TRUE)
    if (length(arg.sets) == 2) {
        arg.sdata <- arg.sets[2]
        arg.sdata <- ReadMyFile(arg.sdata, TRUE)
    }
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Missing Values (NAs)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (TRUE %in% unique(as.vector(is.na(arg.data)))) {
    arg.data <- ReplaceNAs(arg.data)
    
}

if (exists("arg.sdata")) {
    if (TRUE %in% unique(as.vector(is.na(arg.sdata)))) {
        arg.sdata <- ReplaceNAs(arg.sdata)
    }
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Meta-datasets
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (is.null(opt$metadata)){
    stop("\n- Argument metadata (-m) is missing. Metadata provided must be an excel file with minimum two columns named 'ids' (sample names matching those in the object above) and 'group' (diagnosis, tumor stage, ect.).\n")
} else {
    arg.metasets <- SplitList(opt$metadata)
    arg.metadata <- arg.metasets[1]
    arg.metadata <- ReadMyFile(arg.metadata, FALSE)
    if (length(arg.metasets) == 2) {
        arg.smetadata <- arg.metasets[2]
        arg.smetadata <- ReadMyFile(arg.smetadata, FALSE)
    }
    
    # IDs and Groups to check user input!
    if (is.null(opt$groups)){
        stop("Argument -g should be specified as a comma separated list of length two (without quotes and parenthesis!). The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.")
    } else {
        arg.groups <- SplitList(opt$groups)
        if (length(arg.groups) < 2) {
                stop("Argument -g should be specified as a comma separated list of length two (without quotes and parenthesis!). The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.")
        }
    }
    
    
        
    # IDs and Groups to contrast
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    # IDs
    acall <- parse(text = paste0("arg.metadata$", as.character(arg.groups[[1]])))
    arg.ids <- as.character(eval(acall))
    
    if (length(arg.ids) <= 1) {
        stop(paste0("No column in metadata file called ",arg.groups[[1]]))
    } else {
        arg.metadata$ids <- arg.ids
    }
    
    
    # Match Data and Metadata
    arg.metadata <- arg.metadata[arg.metadata$ids %in% colnames(arg.data),]
    
    # Groups
    acall <- parse(text = paste0("arg.metadata$", as.character(arg.groups[[2]])))
    arg.group <- as.factor(as.character(eval(acall)))
    
    if (length(arg.group) <= 1) {
        stop(paste0("No column in metadata file called ",arg.groups[[2]]))
    }
    
    
    
        
    if (length(arg.groups) == 4) {
            
        # IDs
        acall <- parse(text = paste0("arg.smetadata$", as.character(arg.groups[[3]])))
        arg.ids <- as.character(eval(acall))
            
        if (length(arg.ids) <= 1) {
            stop(paste0("No column in metadata file called ",arg.groups[[3]]))
        } else {
            arg.smetadata$ids <- arg.ids
        }
            
        # Match Data and Metadata
        arg.smetadata <- arg.smetadata[arg.smetadata$ids %in% colnames(arg.sdata),]
            
        # Groups
        acall <- parse(text = paste0("arg.smetadata$", as.character(arg.groups[[4]])))
        arg.sgroup <- as.factor(as.character(eval(acall)))
            
        if (length(arg.sgroup) <= 1) {
            stop(paste0("No column in metadata file called ",arg.groups[[4]]))
        }
    }
    
    
    # Batches
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if (is.null(opt$batches)){
        arg.databatch <- FALSE
        arg.sdatabatch <- FALSE
    } else {
        arg.batches <- SplitList(opt$batches)
        acall <- parse(text = paste0("arg.metadata$", as.character(arg.batches[[1]])))
        arg.batch <- as.factor(as.character(eval(acall)))
        arg.databatch <- TRUE
        if (length(arg.batch) <= 1) {
            stop(paste0("No column in metadata file called ",as.character(arg.batches[[1]])))
        }
        if (length(arg.batches) > 1 & exists("arg.smetadata")) {
            acall <- parse(text = paste0("arg.smetadata$", as.character(arg.batches[[2]])))
            arg.sbatch <- as.factor(as.character(eval(acall)))
            arg.sdatabatch <- TRUE
            if (length(arg.sbatch) <= 1) {
                stop(paste0("No column in metadata file called ",as.character(arg.batches[[2]])))
            }
        } else {
            arg.sdatabatch <- FALSE
        }
    }
}
    
    


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### ADDITIONAL ARGUMENTS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Variant (datatype)
if (is.null(opt$variant)){
    stop("\n- Argument data (-v) is missing. -v specifies type of input data (and sdata, if included).\n")
} else {
    arg.variant <- SplitList(opt$variant)
}



# Standardize data
if (is.null(opt$standardize)){
    arg.standardize <- c("none", "none")
} else {
    arg.standardize <- SplitList(opt$standardize)
}



# Transform
if (is.null(opt$transform)){
    arg.transform <- c("none", "none")
} else {
    arg.transform <- SplitList(opt$transform)
}


# Check data distribution
if (is.null(opt$datacheck)){
    arg.datacheck <- TRUE
} else {
    arg.datacheck <- opt$datacheck
}


# MDS plot
if (is.null(opt$plotmds)){
    arg.plotmds <- FALSE
} else {
    arg.plotmds <- opt$plotmds
}


# Kmeans
if (is.null(opt$kmeans)){
    arg.kmeans <- FALSE
} else {
    arg.kmeans <- TRUE
    if (opt$kmeans == TRUE) {
        labels.kmeans <- ""
    } else {
        file <- try(labels.kmeans <- as.character(eval(parse(text = paste0("arg.metadata$", as.character(opt$kmeans))))))
        if (class(file) == "try-error") {
            labels.kmeans <- ""
            rm(file)
        }
    }
}



# Significance
if (is.null(opt$sig)){
   cat("\n- No cut-off for significant hits has been chosen. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
   arg.logFC <- 1
   arg.FDR <- 0.05
   arg.slogFC <- 1
   arg.sFDR <- 0.05
} else {
    arg.sig <- as.numeric(SplitList(opt$sig))
    if (length(arg.sig) > 1) {
        arg.logFC <- arg.sig[1]
        arg.FDR <- arg.sig[2]
        if (length(arg.sig) > 3) {
            arg.slogFC <- arg.sig[3]
            arg.sFDR <- arg.sig[4]
        }
    } else {
        stop("If argument -f is set, it must be a comma separated list of length 2 OR 2*2 = 4 , if two datasets are used, (without quotes and parenthesis!) where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr) for each set. If -f is not specified defaults will be used. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
    }
}



# Colors
if (is.null(opt$colors)){
    arg.colors <- viridis(length(levels(arg.group)), begin = 0.2, end = 0.8)
} else {
    arg.colors <- SplitList(opt$colors)
}


# Filename
if (is.null(opt$filename)){
    arg.filename <- "Results"
} else {
    arg.filename <- opt$filename
}



# Heatmap
if (is.null(opt$plotheatmap)){
    arg.plotheatmap <- NULL
} else {
    arg.plotheatmap <- opt$plotheatmap
}



# Correlation
if (is.null(opt$corr)){
    arg.corrby <- NULL
} else if (opt$corr %in% must.be.type) {
    arg.corrby <- opt$corr
} else {
    stop(paste0("Argument -o (correlation analysis) is specified. This argument takes a string denoting which set of variables to use for correlation analysis, options are: ", must.be.type,"."))
}



# LASSO
if (is.null(opt$lasso)){
    arg.lasso <- NULL
} else {
    arg.lasso <- opt$lasso
}



# WGCNA

if (is.null(opt$WGCNA)){
    arg.WGCNA <- NULL
} else if (opt$WGCNA %in% c("DA", "DE", "ALL")) {
    arg.WGCNA <- opt$WGCNA
} else {
    stop("WGCNA may be performed with either results of differential abundance / expression analysis (DA / DE) or with all variables (ALL). N.B It is not advisable to run WGCNA with all variables if n > 5000. This will make WGCNA slow and plots will be difficult to iterpret. If 'ALL' is chosen and n > 5000, NO plots will be generated, but module variable interconnectivity scores will still be computed.")
}



# CutoffWGCNA
if (is.null(opt$cutoffWGCNA)){
    arg.cutoffWGCNA <- c(min(10, nrow(arg.data)/2), 25, 25)
} else {
    arg.cutoffWGCNA <- as.numeric(SplitList(opt$cutoffWGCNA))
}



# Survival Analysis
if (is.null(opt$survival)){
    arg.survival <- NULL
} else {
    arg.survival <- opt$survival
    if ((!arg.survival %in% must.be.type)) {
        stop("Options for survival analysis variable sets are; DA, LASSO, EN or Consensus. Please re-run pipeline with one of these!")
    }
}



# Covariates (DEA and survival)
if (is.null(opt$covar)){
    arg.covarD <- NULL
    arg.scovarD <- NULL
    arg.covarS <- NULL
} else {
    arg.covar <- SplitList(opt$covar)
    arg.covarS <-  arg.covar[-1]
    if (arg.covar[1] == TRUE) {
        arg.covarD <- arg.covar[-1]
        if (exists("arg.sdata")) {
            arg.scovarD <- arg.covar[-1]
        }
    } else if (arg.covar[1] == FALSE) {
        arg.covarD <- NULL
        arg.scovarD <- NULL
    } else {
        stop("First argument in '-r' must be TRUE or FALSE. If TRUE, covariates will be used for both DE analysis and survival analysis. If FALSE, covariates will be used only for survival analysis.")
    }
}





# Stratify
if (is.null(opt$stratify)){
    arg.stratify <- NULL
} else {
    arg.stratify <- SplitList(opt$stratify)
}



# Survival Plot
if (is.null(opt$survplot)){
    arg.survplot <- 50
} else {
    arg.survplot <- opt$survplot
}





# Network Analysis

if (is.null(opt$PPint)){
    arg.PPI <- NULL
} else {
    arg.PPI <- SplitList(opt$PPint)
    if (!arg.PPI[[1]] %in% approvedGeneIDs | !arg.PPI[[2]] %in% genequery | length(arg.PPI) != 2) {
        stop(paste0("Argument x must be a comma separated list (no quote or parenthesis) of length two. First element in list must specify type of gene ID matching the type of IDs in expression file, approved options are: ", approvedGeneIDs, ".Second element must specify which PPI database to use, currently options are: ", genequery, "."))
    }
}


if (is.null(opt$GenemiRint)){
    arg.GmiRI <- NULL
} else {
    arg.GmiRI <- SplitList(opt$GenemiRint)
    if (!arg.GmiRI[[1]] %in% approvedmiRIDs | !arg.GmiRI[[2]] %in% miRNAquery | length(arg.GmiRI) != 2) {
        stop(paste0("Argument x must be a comma separated list (no quote or parenthesis) of length two. First element in lit must specify type of miRNA ID matching the type of IDs in expression file, approved options are: ", approvedmiRIDs, ".Second element must specify which miRNA-gene database to use, currently options are: ", miRNAquery, "."))
    }
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create directory Results
dir.create(arg.filename)
setwd(paste0(arg.filename, "/"))





# Check if data contains zeros and negative values.
hasZeroD <- unique(as.vector(arg.data == 0))
hasNegD <- unique(as.vector(arg.data < 0))

if(arg.transform[1] %in% c("log2", "log10", "logit")) {
    if (TRUE %in% hasNegD) {
        stop("\n- Data contains negative values and cannot be log transformed. Re-run command WITHOUT flag -t, or alternatively if using two datasets, specify 'none' as the -t input for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
    } else {
        if (TRUE %in% hasZeroD) {
            arg.data.original <- arg.data
            arg.data <- ReplaceZero(arg.data, arg.group)
        }
    }
}




if (exists("arg.sdata")){
    hasZeroS <- unique(as.vector(arg.sdata == 0))
    hasNegS <- unique(as.vector(arg.sdata < 0))
}

if(exists("arg.sdata") & arg.transform[2] %in% c("log2", "log10", "logit")) {
    if (TRUE %in% hasNegS) {
        stop("\n- Second dataset contains negative values and cannot be log transformed. Re-run command WITHOUT flag -t, or alternatively if using two datasets, specify 'none' as the -t input for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
    } else {
        if (TRUE %in% hasZeroS) {
            arg.sdata.original <- arg.sdata
            arg.sdata <- ReplaceZero(arg.sdata, arg.group)
        }
    }
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                         ## Normalization, Filtering and Transformation ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataset

if (exists("arg.data.original")) {
    arg.data <- NormalizeData(arg.variant[1], arg.data, arg.group, arg.transform[1], arg.standardize[1], arg.data.original)
} else {
    arg.data <- NormalizeData(arg.variant[1], arg.data, arg.group, arg.transform[1], arg.standardize[1])
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Second Dataset

if(exists("arg.sdata")) {
    if (length(arg.variant) < 2) {
        stop("\nTwo datasets are input for correlation analysis, BUT argument -v only has length one. Length of -v must be two, see -h.\n")
    }
    if (length(arg.transform) < 2) {
        stop("\nTwo datasets are input for correlation analysis, BUT argument -t only has length one. Length of -t must be two, see -h.\n")
    }
}



if (exists("arg.sdata")) {
    if (exists("arg.sdata.original")) {
        arg.sdata <- NormalizeData(arg.variant[2], arg.sdata, arg.sgroup, arg.transform[2], arg.standardize[2], arg.sdata.original)
    } else {
        arg.sdata <- NormalizeData(arg.variant[2], arg.sdata, arg.sgroup, arg.transform[2], arg.standardize[2])
    }
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### BATCH CORRECTION ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



if (arg.databatch == TRUE){
    if (length(arg.batch) > 0) {
        arg.design <-  model.matrix(~arg.group)
        
        if (arg.variant[1] == "seq") {
            data.batch <- ComBat(as.matrix(arg.data$E), arg.batch, arg.design, par.prior=TRUE,prior.plots=FALSE)
        } else {
            data.batch <- ComBat(as.matrix(arg.data), arg.batch, arg.design, par.prior=TRUE,prior.plots=FALSE)
        }
        
    } else {
        data.batch <- arg.data
        cat("\n- No column names match specified batchs for dataset.\n")
    }
} else {
    cat("\n- No batch correction requested.\n")
}




if (arg.sdatabatch == TRUE){
    if (length(arg.sbatch) > 0) {
        arg.sdesign <- model.matrix(~arg.sgroup)
        
        if (arg.variant[2] == "seq") {
            sdata.batch <- ComBat(as.matrix(arg.sdata$E), arg.sbatch, arg.sdesign, par.prior=TRUE,prior.plots=FALSE)
        } else {
            sdata.batch <- ComBat(as.matrix(arg.sdata), arg.sbatch, arg.sdesign, par.prior=TRUE,prior.plots=FALSE)
        }
    } else {
        sdata.batch <- arg.sdata
        cat("\n- No column names match specified batchs for second dataset. Continuing without batch correction.\n")
    }
} else {
    cat("\n- No batch correction requested for second dataset.\n")
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                        ### Distributional Checks ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (arg.datacheck == TRUE) {
    if (arg.databatch == TRUE) {
        subset.data <- data.batch[sample(nrow(data.batch), 10),]
    } else {
        if (arg.variant[1] == "seq") {
            subset.data <- arg.data$E[sample(nrow(arg.data$E), 10),]
        } else {
            subset.data <- arg.data[sample(nrow(arg.data), 10),]
        }
    }
    
    list.of.dists <- FitDistributions(subset.data)
    
    dir.create("DataChecks")
    setwd("DataChecks/")
    PlotDistributions(subset.data, list.of.dists)
    setwd("..")
    
    rm(subset.data, list.of.dists)
}



if (arg.datacheck == TRUE & exists("arg.sdata")) {
    if (arg.sdatabatch == TRUE) {
        subset.data <- sdata.batch[sample(nrow(sdata.batch), 10),]
    } else {
        
        if (arg.variant[2] == "seq") {
            subset.data <- arg.sdata$E[sample(nrow(arg.sdata$E), 10),]
        } else {
            subset.data <- arg.sdata[sample(nrow(arg.sdata), 10),]
        }
    }
    
    list.of.dists <- FitDistributions(subset.data)
    
    dir.create("SecondDataChecks")
    setwd("SecondDataChecks/")
    PlotDistributions(subset.data, list.of.dists)
    setwd("..")
    
    rm(subset.data, list.of.dists)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### PRELIMINARY MDS PLOT ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



MDScolors <- gsub(pattern = "FF", replacement = "", x = arg.colors)


if (arg.plotmds == TRUE && arg.databatch == TRUE){
    mdsplot <- MDSPlot(data.batch, arg.group, "", MDScolors)
    ggsave(paste0(arg.filename, "_MDSplot_batchcorr.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

} else if (arg.plotmds == TRUE && arg.databatch == FALSE){
    if (arg.variant[1] == "seq") {
        mdsplot <- MDSPlot(data.frame(arg.data$E), arg.group, "", MDScolors)
    } else {
        mdsplot <- MDSPlot(arg.data, arg.group, "", MDScolors)
    }
    ggsave(paste0(arg.filename, "_MDSplot.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)
    
    rm(mdsplot)
    
} else {
    cat("\n- No preliminary plot requested.\n")
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                 ### Kmeans ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





if (arg.kmeans == TRUE) {
    
    if (arg.variant[1] == "seq") {
        k.data <- arg.data$E
    } else {
        k.data <- arg.data
    }
    
    # Number of sample sets to generate
    nsets <- 1:ceiling(nrow(k.data)/1000)
    if(length(nsets) > 10) {
        nsets <- 1:10
    }
    
    # Number of variables in each sample set
    setsize <- nrow(k.data)
    if (setsize > 2000) {
        setsize <- 2000
    }

    
    # Number of kmeans to try
    if(ncol(k.data) <= 100) {
        nks <- 2:6
    } else if (ncol(k.data) > 100 && ncol(k.data) <= 500) {
        nks <- 2:11
    } else {
        nks <- 2:16
    }
    
    
    cat(paste0("Based on size of dataset, ", length(nsets), " sample sets will be generated of size ", setsize, " and ", length(nks), " clusters will be tested - N.B This may take up to 10 min!\nRunning......"))
    
    dir.create("KmeansResults")
    setwd("KmeansResults/")
    
    
    list.of.dfs <- list()
    
    if (arg.databatch == TRUE) {
        for (idx in 1:length(nsets)) {
            df <- t(data.batch[sample(nrow(data.batch), setsize), ])
            list.of.dfs[[idx]] <- df
        }
        Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
        Kmeans.Out <- PlotKmeans(data.batch, Kmeans.list, nks, labels.kmeans, arg.filename)
    } else {
        for (idx in 1:length(nsets)) {
            df <- t(k.data[sample(nrow(k.data), setsize), ])
            list.of.dfs[[idx]] <- df
        }
        Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
        Kmeans.Out <- PlotKmeans(k.data, Kmeans.list, nks, labels.kmeans, arg.filename)
    }
    
    
    out <- cbind(arg.metadata, Kmeans.Out)
    
    write.table(out, paste0(arg.filename,"_Metadata_Kmeans.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    setwd("..")
    
    rm(k.data, Kmeans.list, Kmeans.Out, out)
}






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CLEAN UP AND SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 2
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rm(SplitList, ReadMyFile, ReplaceNAs, ReplaceZero, NormalizeData, FitDistributions, PlotDistributions, MDSPlot, EstimateKmeans)
gc(full = TRUE)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### DIFFERENTIAL EXPRESSION ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Differential Expression Analysis with Limma


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# First dataset

dir.create("DEAAResults")
setwd("DEAAResults/")
    

# Make design matrix
if (arg.databatch == "FALSE") {
    arg.design.str <- "model.matrix(~0+arg.group"
    out.name <- "_DE"
} else if (length(arg.batch) > 0) {
    arg.design.str <- "model.matrix(~0+arg.group+arg.batch"
    out.name <- "_databatch_DE"
} else {
    stop("Batch correction selected but no batches column found!")
}
    
if(is.null(arg.covarD)) {
    arg.design <- eval(parse(text=paste0(arg.design.str, ")")))
} else {
    if (length(arg.covarD) == 1) {
        df <- data.frame(arg.metadata[,colnames(arg.metadata) %in% arg.covarD])
        colnames(df) <- arg.covarD
    } else {
        df <- arg.metadata[,colnames(arg.metadata) %in% arg.covarD]
    }
        
    s <- lapply(split(as.matrix(df), col(df)), factor)
    my.names <- paste0("arg.", colnames(df))
    list2env(setNames(s, my.names), envir=.GlobalEnv)
    my.names <- paste0(my.names, collapse = "+")
    arg.design <- eval(parse(text=paste0(arg.design.str,"+",my.names,")")))
    rm(s)
}



# Making group contrasts
combinations <- data.frame(t(combn(paste0("arg.group", levels(arg.group)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(arg.design))))

    
# Apply DE_limma function to all comparisons
res.DE <- DAFeatureApply(contrast.matrix, arg.data, arg.design, arg.logFC, arg.FDR, NULL, FALSE)


# Write results out as excel file
if (!is.null(res.DE)) {
    #DE.out <- ExcelOutput(res.DE, paste0(arg.filename, out.name))
    DE.out <- TextOutput(res.DE, paste0(arg.filename, out.name))
    rownames(DE.out) <- NULL
    res.DE.names <- unique(DE.out$name)
    rm(res.DE)
} else {
    cat("No signficant DE/DA hits found. Check output file from differential expression analysis. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
}

setwd("..")




if (arg.variant[1] == "seq") {
    arg.data <- data.frame(arg.data$E)
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Second dataset


if(!is.null(arg.PPI) & !is.null(arg.GmiRI)) {
    
    combinations <- data.frame(t(combn(paste0("arg.sgroup", levels(arg.sgroup)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    
    
    if (arg.sdatabatch == FALSE) {
        arg.sdesign.str <- "model.matrix(~0+arg.sgroup"
        out.name <- "_second_DE"
    } else if (length(arg.sbatch) > 0) {
        arg.sdesign.str <- "model.matrix(~0+arg.sgroup+arg.sbatch"
        out.name <- "_second_databatch_DE"
    } else {
        stop("Batch correction selected but no batches column found!")
    }
    
    if(is.null(arg.scovarD)) {
        arg.sdesign <- eval(parse(text=paste0(arg.sdesign.str, ")")))
    } else {
        if (length(arg.scovarD) == 1) {
            df <- data.frame(arg.smetadata[,colnames(arg.smetadata) %in% arg.scovarD])
            colnames(df) <- arg.scovarD
        } else {
            df <- arg.smetadata[,colnames(arg.smetadata) %in% arg.scovarD]
        }
        
        s <- lapply(split(as.matrix(df), col(df)), factor)
        my.names <- paste0("arg.s", colnames(df))
        list2env(setNames(s, my.names), envir=.GlobalEnv)
        my.names <- paste0(my.names, collapse = "+")
        arg.sdesign <- eval(parse(text=paste0(arg.sdesign.str,"+",my.names,")")))
        rm(s)
    }
    
    # Making group contrasts
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(arg.sdesign))))
    
    # Apply DE_limma function to all comparisons
    res.sDE <- DAFeatureApply(contrast.matrix, arg.sdata, arg.sdesign, arg.slogFC, arg.sFDR, NULL, FALSE)
    
    # Write results out as excel file
    if (!is.null(res.sDE)) {
        #sDE.out <- ExcelOutput(res.sDE, paste0(arg.filename, out.name))
        sDE.out <- TextOutput(res.sDE, paste0(arg.filename, out.name))
        rownames(sDE.out) <- NULL
        res.sDE.names <- unique(sDE.out$name)
    } else {
        cat("No signficant DE/DA hits found for analysis of second dataset. Check output file from differential expression analysis. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }
}




if (exists("arg.sdata") && arg.variant[2] == "seq") {
    arg.sdata <- data.frame(arg.sdata$E)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                                       ## LASSO Regression ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (!is.null(arg.lasso)) {
    if(arg.lasso <= 0.0 || arg.lasso > 1.0 ) {
        stop("\n- The input for flag -l denotes hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression. Re-run the pipeline again with correct -l input or remove -l all together.\n")
    }
    
    
    # Length of each group
    len <- as.numeric(table(arg.group))
    test.train <- unique(len >= 19)
    too.few <- unique(len < 9)
    
    
    # Stop Lasso if too few samples
    if (TRUE %in% too.few) {
        stop("\n- LASSO cannot be performed, too few samples per group, minimum is 10!\n")
    }
    
    
    dir.create("LASSOResults")
    setwd("LASSOResults/")
    
    group.LASSO <- arg.group
    seeds <- sample(1:1000, 10)
    LASSO.res <- list()
    
    cat("Cross-validation for grouped multinomial LASSO is running with 10 random seeds, this will take some minutes...")
    
    if(FALSE %in% test.train) {
        if (arg.databatch == TRUE) {
            if (length(levels(as.factor(group.LASSO))) > 2) {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, arg.lasso, FALSE ,TRUE)
                    LASSO.res[[idx]] <-  LR
                }
            } else {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, arg.lasso, FALSE, FALSE)
                    LASSO.res[[idx]] <-  LR
                }
            }
        } else {
            if (length(levels(as.factor(group.LASSO))) > 2) {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], arg.data, group.LASSO, arg.lasso, FALSE ,TRUE)
                    LASSO.res[[idx]] <-  LR
                }
            } else {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], arg.data, group.LASSO, arg.lasso, FALSE ,FALSE)
                    LASSO.res[[idx]] <-  LR
                }
            }
        }
    } else {
        if (arg.databatch == TRUE) {
            if (length(levels(as.factor(group.LASSO))) > 2) {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, arg.lasso, TRUE ,TRUE)
                    LASSO.res[[idx]] <-  LR
                }
            } else {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, arg.lasso, TRUE, FALSE)
                    LASSO.res[[idx]] <-  LR
                }
            }
        } else {
            if (length(levels(as.factor(group.LASSO))) > 2) {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], arg.data, group.LASSO, arg.lasso, TRUE ,TRUE)
                    LASSO.res[[idx]] <-  LR
                }
            } else {
                for (idx in 1:length(seeds)) {
                    LR <- LASSOFeature(seeds[[idx]], arg.data, group.LASSO, arg.lasso, TRUE ,FALSE)
                    LASSO.res[[idx]] <-  LR
                }
            }
        }
    }
    
    
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Extract results of 10 runs - Write out and plot results
    VarsSelect <- Reduce(intersect, lapply(LASSO.res, '[[', 1))
    
    if (length(VarsSelect) < 2) {
        stop("\n- There is no overlap in 10 elastic net runs. If you ran LASSO (-l was et to 1.0) you can try and relax alpha and perform elastic net instead (0.0 < -l < 1.0). Otherwise you data may have to high of a noise ratio to sample size, LASSO should not be performed.\n")
    }
    
    
    VarsSelect <- data.frame(VarsSelect[-1])
    colnames(VarsSelect) <- c("LASSO.Var.Select")
    

    # Write out LASSO/EN results
    write.table(VarsSelect, paste0(arg.filename,"_LASSO.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Consensus DEA and LASSO
    consensus <- DE.out[DE.out$name %in% VarsSelect$LASSO.Var.Select,]
    
    if (nrow(consensus) > 0) {
        write.table(consensus, paste0(arg.filename,"_DEA_LASSO_Consensus.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        pdf(paste0(arg.filename, "_overlap_DEAA_LASSO_EN.pdf"), height=8, width=12)
        if (length(levels(arg.group)) == 2) {
            venn <- venn.diagram(list(A=unique(as.character(DE.out[DE.out$dir =="up",]$name)), B=unique(as.character(DE.out[DE.out$dir =="down",]$name)), C=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis Up", "D(EA) Analysis Down", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(3, begin = 0.2, end = 0.8, option="cividis"))
            
        } else {
            venn <- venn.diagram(list(A=unique(as.character(DE.out$name)), B=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis All", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(2, begin = 0.2, end = 0.8, option="cividis"))
            
        }
        grid.draw(venn)
        dev.off()
    } else {
        cat("There is no consensus between LASSO regression and DEA/DAA.")
    }
    
    
    
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Cross Validation errors
    LassoRun <- paste0(rep("Run", 10), 1:10)
    CrossValErrormean <- round(unlist(lapply(LASSO.res, '[[', 2)), digits = 4)
    cat(paste0("\nThe average leave-one-out cross validation error for LASSO/elastic-net was: ", mean(CrossValErrormean), "% and the higest error returned from any of the 10 runs was: ", max(CrossValErrormean),"%. Generally the cross validation error should be low ~ 5.0 %, as large errors may indicate a poor model and/or very heterogeneous data. On the other hand, an error of 0 might indicate over-fitting. See CAMPP manual for specifics.\n\n"))
    pCVEM <- data.frame(cbind(CrossValErrormean, LassoRun))
    pCVEM <- ggplot(data=pCVEM, aes(x=LassoRun, y=CrossValErrormean)) + geom_bar(aes(fill = as.factor(LassoRun)), stat="identity") + theme_minimal() + scale_x_discrete(limits=c(LassoRun)) + scale_fill_viridis(begin = 0.0, end = 0.0, discrete=TRUE, option="cividis" ) + theme(legend.position="none") + ylab("CrossValErrormean in %") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16))
    ggsave(paste0(arg.filename, "_CrossValidationPlot.pdf"), plot = pCVEM)
    
    
    
    
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Area under the curve AUC
    if(TRUE %in% test.train) {
        
        ll <- list()
        llev <- levels(as.factor(group.LASSO))
        
        for (idx in 1:length(llev)) {
            pos <- which(group.LASSO == as.character(llev[idx]))
            ll[[idx]] <- pos
        }
        
        my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))
        
        
        if (arg.databatch == TRUE) {
            LASSO.data <- data.batch
        } else {
            LASSO.data <- arg.data
        }
        
        
        testD <- data.frame(t(LASSO.data[rownames(LASSO.data) %in% as.character(VarsSelect$LASSO.Var.Select), my.samp]))
        testG <- as.integer(group.LASSO[my.samp])
        
        trainD <- data.frame(t(LASSO.data[rownames(LASSO.data) %in% as.character(VarsSelect$LASSO.Var.Select), -my.samp]))
        trainG <- as.integer(group.LASSO[-my.samp])
        
        mn.net <- nnet::multinom(trainG ~ ., data=trainD)
        mn.pred <- predict(mn.net, newdata=testD, type="prob")
        roc.res <- multiclass.roc(testG, mn.pred)
        roc.res <- data.frame(round(as.numeric(sub(".*: ", "", roc.res$auc)), digits = 2))
        colnames(roc.res) <- "AUC"
        cat(paste0("Are under the curve (AUC) for variables selected from 10 LASSO/EN runs was: ", roc.res$AUC))
        write.table(roc.res, paste0(arg.filename,"_AUC.txt"), row.names=FALSE, col.names = TRUE, quote = FALSE)
    }
    

    setwd("..")
    try(rm(VarsSelect, LR, venn, CrossValErrormean, mn.net, mn.pred, roc.res), silent=T)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                    ### Plotting Results Heatmap ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (!is.null(arg.plotheatmap)) {
    # color scheme
    arg.colors.hm <- GetColors(as.character(arg.group), arg.colors)
    
    if (arg.databatch == TRUE){
        arg.hm <- data.batch
        name <- "_batchcorr"
    } else {
        arg.hm <- arg.data
        name <- ""
    }
    
    
    if (arg.plotheatmap == "ALL") {
        stop("\nOption ALL is not allowed for heatmap, too heavy! Pick either 'DE', 'DA', 'LASSO', 'EN' or 'Consensus'.\n")
    } else if (arg.plotheatmap %in% c("DA", "DE")) {
        arg.hm <- arg.hm[rownames(arg.hm) %in% res.DE.names,]
    } else {
        if(!is.null(arg.lasso)) {
            if (arg.plotheatmap == "Consensus") {
                arg.hm <- arg.hm[rownames(arg.hm) %in% as.character(consensus$name),]
            }
            if (arg.plotheatmap %in% c("EN", "LASSO")) {
                arg.hm <- arg.hm[rownames(arg.hm) %in% as.character(VarsSelect$LASSO.Var.Select),]
            }
        } else {
            stop("You have specified argument -a which will produce a heatmap with results from either DA/DE analysis, LASSO/Elastic-Net Regression or the Consensus of these. Input to argument -a must be a string specifying which results to plot, options are: DA, DE, LASSO, EN, Consensus")
        }
    }
    

    # heatmap colors in blue
    arg.hm.gradient <- viridis(300, option="cividis")
    arg.range <- c(round(min(DE.out$logFC)), round(max(DE.out$logFC)))
    
    # Heatmap as pdf
    MakeHeatmap(arg.hm, arg.hm.gradient, arg.colors.hm, arg.colors, arg.group, arg.filename, arg.range)
    
    rm(arg.hm)
    
} else {
    cat("\n- No heatmap requested.\n")
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### CORRELATION ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Data (TIF) and sdata correlations

if (exists("arg.sdata") & !is.null(arg.corrby)) {
    
    dir.create("CorrelationResults")
    setwd("CorrelationResults/")
    
    
    retainedSamp <- intersect(colnames(arg.data), colnames(arg.sdata))
    
    if (arg.corrby == "ALL") {
        retainedvar <- intersect(rownames(arg.data), rownames(arg.sdata))
    } else if (arg.corrby %in% c("DA", "DE")) {
        retainedvar <- intersect(res.DE.names, rownames(arg.sdata))
    } else {
        if(!is.null(arg.lasso)) {
            if (arg.survival == "Consensus" ) {
                retainedvar <- intersect(as.character(consensus$name), rownames(arg.sdata))
            } else {
                retainedvar <- intersect(as.character(VarsSelect$LASSO.Var.Select), rownames(arg.sdata))
            }
        } else {
            stop("You have chosen to perform correlation analysis with results from LASSO/EN BUT this analysis has not been performed. Please re-run pipeline with parameter -l (see Manual or help (-h))")
        }
    }
    
    
    if (arg.databatch == TRUE) {
        data.corr <- data.batch[rownames(data.batch) %in% retainedvar, colnames(data.batch) %in% retainedSamp]
    } else {
        data.corr <- arg.data[rownames(arg.data) %in% retainedvar, colnames(arg.data) %in% retainedSamp]
    }
    
    if (arg.sdatabatch == TRUE) {
        corr.out <- sdata.batch[rownames(sdata.batch) %in% retainedvar, colnames(sdata.batch) %in% retainedSamp]
    } else {
        corr.out <- arg.sdata[rownames(arg.sdata) %in% retainedvar, colnames(arg.sdata) %in% retainedSamp]
    }
    
    # Perform correction analysis and generate overall correlation plot
    res.corr <- CorrAnalysis(data.corr, corr.out, arg.filename)
    
    
    # print out significant hits in Excel
    res.corr$sig.corr <- ifelse(res.corr$fdr <= arg.FDR, "yes", "no")
    
    write.table(res.corr, paste0(arg.filename,"_corr.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


    # Individual correlation plots
    corr.features <- as.character(res.corr[which(res.corr$sig.corr == "yes"),]$name)
    if (length(corr.features) > 1) {
        CorrelationPlots(data.corr, corr.out, corr.features, arg.filename)
    } else {
        cat("\n- No significant correlations, no plotting.\n")
    }
    setwd("..")
    try(rm(res.corr,corr.out,data.corr), silent=T)
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### SURVIVAL ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






### SETTING UP DATA ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (!is.null(arg.survival)) {
    
    if (length(intersect(colnames(arg.metadata), must.contain)) < 3) {
        stop("Survival analysis requested, but one or more columns; 'survival', 'outcome', 'outcome.time' are missing from metadata file!")
    }
    
    
    dir.create("SurvivalResults")
    setwd("SurvivalResults/")
    
    # Setting up dataset for survival analysis, matching samples either batch or no batch. Only DA/DE featrues are used.
    metadata.surv <- arg.metadata[which(arg.metadata$survival == 1),]
    arg.samples <- metadata.surv$ids
    
    if (arg.databatch == TRUE) {
        data.surv <- data.batch[,colnames(data.batch) %in% arg.samples]
    } else {
        data.surv <- arg.data[,colnames(arg.data) %in% arg.samples]
    }
    
    
    
    if (arg.survival == "ALL") {
        cat("\nYou have chosen to use all varibles for survival analysis, this is NOT advisable! See options for -u in with the help flag -h.\n")
    } else if (arg.survival %in% c("DA", "DE")) {
        data.surv <- data.surv[rownames(data.surv) %in% res.DE.names,]
    } else {
        if(!is.null(arg.lasso)) {
            if (arg.survival == "Consensus" ) {
                data.surv <- data.surv[rownames(data.surv) %in% as.character(consensus$name),]
            } else {
                data.surv <- data.surv[rownames(data.surv) %in% as.character(VarsSelect$LASSO.Var.Select),]
            }
        } else {
            stop("You have chosen to perform survival analysis with results from LASSO/EN but this analysis has not been performed. Please re-run pipeline with parameter -l (see Manual or help (-h))")
        }
    }
    
    
    
    
    
    # If user has specified covariates for survivalanalysis, extract these.
    if(is.null(arg.covarS)) {
        surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$outcome)
        colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "outcome")
    } else {
        my.covars <- metadata.surv[,colnames(metadata.surv) %in% arg.covarS]
        surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$outcome, my.covars)
        colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "outcome", arg.covarS)
    }
    
   
    # datadist format for cph function with cubic splines.
    features <- rownames(data.surv)
    dd <- datadist(surv_object); options(datadist='dd')
    
    mylist <- list(features, dd, surv_object)
    
    
    
    
    ### User Specified Covariates ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    # Evaluating class of covariate, continious or categorical.
    arg.covarS.original <- arg.covarS
    
    if(!is.null(arg.covarS)) {
        for (i in 1:length(arg.covarS)) {
            if(class(eval(parse(text=paste0("surv_object$", arg.covarS[i])))) %in% c("integer", "numeric")) {
                arg.covarS[i] <- paste0("rcs(", arg.covarS[i], ")")
            }
        }
        
        arg.covarS <- paste(arg.covarS,collapse="+")
    }
    
    ### Liniarity of Continious Covariates ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    # Survival models with or without user specified covariates.
    covariate_linearity <- list()
    
    if(!is.null(arg.covarS)) {
        for (f in features) {
            acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + ", arg.covarS, ", data = surv_object, x=TRUE,y=TRUE)"))
            eval(acall)
            covariate_linearity[[as.character(f)]] <- result
        }
    } else {
        for (f in features) {
            acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),"), data = surv_object, x=TRUE,y=TRUE)"))
            eval(acall)
            covariate_linearity[[as.character(f)]] <- result
        }
    }
    
    
    # Check for liniarity violation:
    check1 <- which(as.numeric(unlist(lapply(covariate_linearity, function(x) length(x)))) == 1)
    
    if (length(check1) > 0) {
        surv_object <- surv_object[,-check1]
        covariate_linearity <- covariate_linearity[-check1]
    }
    
    # Check for error, other.
    check2 <- lapply(covariate_linearity, function(i) tryCatch({ anova(i) }, error=identity))
    check2 <- which(as.character(vapply(check2, is, logical(1), "error")) == TRUE)
    
    
    if (length(check2) > 0) {
        surv_object <- surv_object[,-check2]
        covariate_linearity <- covariate_linearity[-check2]
    }
    
    
    # Remove checks to save memory
    rm(check1, check2)
    
    # Re-assign features
    features <- colnames(surv_object)
    
    # Anova comparing linear covariate and non-linear covariate.
    covariate_linearity  <- lapply(covariate_linearity, function(x) anova(x))
    
    # Unlisting and extracting pvalues for non-linear fit.
    covariate_linearity <- data.frame(do.call(rbind, lapply(covariate_linearity, function(x) x[grep("Nonlinear|NONLINEAR", rownames(x)),3])))
    
    # FDR correction for multiple testing
    covariate_linearity <- apply(covariate_linearity, 2, function(x) p.adjust(x, method = "fdr", n=nrow(covariate_linearity)))
    
 
    # Colnames with and without user-specified covariates.
    if(is.null(arg.covarS)) {
        colnames(covariate_linearity) <- c("Feature")
    } else {
        covar.names <- unlist(strsplit(arg.covarS, split="[+]"))
        if(length(grep("rcs", covar.names)) > 0) {
            covar.names <- c("Feature", covar.names[grep("rcs", covar.names)], "Total")
            covar.names <- gsub("rcs\\(|\\)", "", covar.names)
            colnames(covariate_linearity) <- covar.names
        } else {
            colnames(covariate_linearity) <- c("Feature")
        }
    }
    
 
    
    # Extracting p-values and filtering for significance.
    covariate_nonlinear <- colnames(covariate_linearity)[apply(covariate_linearity, 2, function(x) any(x < 0.05))]
    
    # Updating covariates with cubic splines.
    if (length(covariate_nonlinear) > 0) {
        cat(paste0("\n- The following continious covariate(s) may be violating the assumption of linearity:  ", covariate_nonlinear,".\nCubic splines will be added.\n"))
        if(!is.null(arg.covarS.original)) {
            for (i in 1:length(arg.covarS.original)) {
                if (arg.covarS.original[i] %in% covariate_nonlinear) {
                    arg.covarS.original[i] <- paste0("rcs(", arg.covarS.original[i], ")")
                }
            }
        }
    }
    
    rm(covariate_linearity)
    
   
    
    ### Testing Cox Proportional Hazard Assumption ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    # Picking model based on result of linearity check.
    pha_check <- list()
    
    
    for (f in features) {
        if ("feature" %in% covariate_nonlinear) {
            if(is.null(arg.covarS.original)) {
                acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f), "), data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                eval(acall)
                pha_check[[as.character(f)]] <- result
            } else {
                acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + ", arg.covarS.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                eval(acall)
                pha_check[[as.character(f)]] <- result
            }
        } else {
            if(is.null(arg.covarS.original)) {
                acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ ", as.character(f), ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                eval(acall)
                pha_check[[as.character(f)]] <- result
                
            } else {
                acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ ", as.character(f), "+ ", arg.covarS.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                eval(acall)
                pha_check[[as.character(f)]] <- result
            }
        }
    }
    
    

    
    # Filtering p-values for significance.
    if(length(pha_check[1]) == 1) {
         pha.fail.test <- unlist(lapply(pha_check, function(x) any(x < 0.05)))
         pha.fail.test <- names(pha.fail.test)[pha.fail.test == TRUE]
    } else {
        pha.fail.test <- as.character(unlist(lapply(pha_check, function(x) rownames(x)[apply(x, 1, function(u) any(u < 0.05))])))
        if ("GLOBAL" %in% pha.fail.test) {
            pha.fail.test <- pha.fail.test[-which(pha.fail.test == "GLOBAL")]
        }
    }
    
    

    cat("\nWARNING: The following features and/or covariates failed the test of proportional hazard: ", pha.fail.test, "\nIF the covariates that failed are categorical you may use strata by re-running the pipline adding flag -y followed by the names of the categorical covariates to stratify (if multiple then separate by comma). \nN.B, this pipeline does not handle continuous variables that violate the proportional hazard assumption, if any of these failed PH test, the hazard ratios of these should NOT be evaluated.\n")
    
    
    
    # User stratification of categorical covariates which violate the proportional hazard assumption.
    if(!is.null(arg.stratify)) {
        for (i in 1:length(arg.covarS.original)) {
            if (arg.covarS.original[i] %in% arg.stratify) {
                arg.covarS.original[i] <- paste0("strat(", arg.covarS.original[i], ")")
            }
        }
    }
    
    
    
    
    
    
    
    ### Survival Analysis with updated covariates ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    

    
    # Survival analysis, cox-regression
    
 
    if(!is.null(arg.covarS)) {
        features.surv <- features[-c((length(features)-length(arg.covar)):length(features))]
    } else {
        features.surv <- features[-c((length(features)-1):length(features))]
    }
    
    # Picking model based on result of linearity check and proportional hazard text.
    survival.results <- list()
    

    
    for (f in features.surv) {
        if ("feature" %in% covariate_nonlinear) {
            if(is.null(arg.covarS.original)) {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),"), data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                survival.results[[as.character(f)]] <- result
            } else {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") +", arg.covarS.original, ", data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                survival.results[[as.character(f)]] <- result
            }
        } else {
            if(is.null(arg.covarS.original)) {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), ", data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                survival.results[[as.character(f)]] <- result
                
            } else {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), " +", arg.covarS.original, ", data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                survival.results[[as.character(f)]] <- result
            }
        }
    }
   
   
    # Setting up data and writing out excel file with results and making HR stemplot
    survival.data <- SurvivalCOX(survival.results, arg.filename, arg.survplot)
    
    write.table(survival.data, paste0(arg.filename,"_survival.txt"), sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    setwd("..")
    rm(dd, pha_check, pha.fail.test, survival.results, features.surv, survival.data)
}






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 3
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rm(DAFeature, DAFeatureApply, LASSOFeature, TextOutput, MakeHeatmap, SurvivalCOX, CorrAnalysis, CorrelationPlots)
gc(full = TRUE)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                               ## Weighed Gene Co-Expression Network Analysis ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (!is.null(arg.WGCNA)) {
    
    dir.create("WGCNAResults")
    setwd("WGCNAResults/")
    
    if (arg.databatch == TRUE){
        data.WGCNA <- t(data.batch)
    } else {
        data.WGCNA <- t(arg.data)
    }
    
    if (arg.WGCNA %in% c("DA", "DE")) {
        data.WGCNA <- data.WGCNA[,colnames(data.WGCNA) %in% res.DE.names]
    }

    dir.create("WGCNAPlots")
    setwd("WGCNAPlots/")
    # Check data
    gsg <- goodSamplesGenes(data.WGCNA)
    cat(paste0("- Data set is OK for WGCNA - ", gsg$allOK,".\n"))
    WGCNAres <- WGCNAAnalysis(data.WGCNA, arg.cutoffWGCNA ,arg.filename)
    setwd("..")
    
    if (arg.WGCNA %in% c("DA", "DE")) {
        logFCs <- DE.out
        logFCs$gene <- DE.out$name
        WGCNAres <- lapply(WGCNAres, function(x) merge(x, logFCs, by = "gene"))
        WGCNAres <- lapply(WGCNAres, function(x) x[with(x, order(comparison, Module)),])
    }
    
    mod.names <- names(WGCNAres)
    for (idx in 1:length(WGCNAres)) {
        write.table(WGCNAres[[idx]], paste0(mod.names[idx],"_moduleRes.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    
    setwd("..")
    rm(WGCNAres)
}

gc()





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                               ## Protein-Protein and MiRNA-Gene Network Analysis ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







if (!is.null(arg.PPI)) {
    
    DB <- DownloadPPInt(arg.PPI[[1]])
    
    dir.create("InteractionResults")
    setwd("InteractionResults/")
    dir.create("InteractionPlots")
    
    DE.out$name <- gsub("_", "-", DE.out$name)
    
    PPIs <- GetDeaPPInt(DB, DE.out)
    

   if (!is.null(arg.GmiRI)) {
       
       sDE.out$name <- gsub("_", "-", sDE.out$name)
       
       cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
       GmiRIs <- GetGeneMiRNAInt(arg.GmiRI[[2]], sDE.out, DE.out)
       
       PPGmiRIs <- MatchPPGmiRInt(GmiRIs, PPIs)
       PPGmiRTrim <- TrimWriteInt(PPGmiRIs)
       
       setwd("InteractionPlots/")
       cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
       PlotInt(PPGmiRTrim)
       
       setwd("../..")
        
    } else {
    
    PPTrim <- TrimWriteInt(PPIs)
    
    
    setwd("InteractionPlots/")
    cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
    PlotInt(PPTrim)
        
    setwd("../..")
    
    }
    

} else if (!is.null(arg.GmiRI)) {
    
    dir.create("InteractionResults")
    setwd("InteractionResults/")
    dir.create("InteractionPlots")
    
    DE.out$name <- gsub("_", "-", DE.out$name)
    
    cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
    GmiRIs <- GetGeneMiRNAInt(arg.GmiRI[[2]], DE.out)
    GmiRTrim <- TrimWriteInt(GmiRIs)
    
    setwd("InteractionPlots/")
    cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
    PlotInt(GmiRTrim)
    
    setwd("../..")

} else {
    cat("\nNo Interaction Networks requested\n.")
}


cat("\nCAMPP RUN DONE!\n")








# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Close Renv library.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (is.rlib == TRUE){
    renv::deactivate()
}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# End log file
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gc(full = TRUE)

sink()

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



