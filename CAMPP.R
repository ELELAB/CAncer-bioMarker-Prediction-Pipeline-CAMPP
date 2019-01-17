
#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                            # USER ARGUMENTS
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


sink("./CAMPPlog.txt", append=TRUE, split=TRUE)


library("getopt")


# flag specification # ref: https://tgmstat.wordpress.com/2014/05/21/r-scripts/
spec = matrix(c(
  "help", "h", 0, "logical",
  "data", "d", 1, "character",
  "metadata", "m", 1, "character",
  "datacheck", "k", 2, "logical",
  "serumdata", "s", 2, "character",
  "survival", "u", 2, "logical",
  "survplot", "q", 2, "numeric",
  "transform", "t", 2, "character",
  "databatch", "b", 2, "logical",
  "serumbatch", "e", 2, "logical",
  "filename", "n", 2, "character",
  "logFC", "f", 2, "numeric",
  "FDR", "r", 2, "numeric",
  "survcovar", "p", 2, "character",
  "stratify", "y", 2, "character",
  "plotmds", "o", 2, "logical",
  "colors", "c", 2, "character",
  "plotheatmap", "a", 2, "character",
  "lasso", "l", 2, "logical",
  "WGCNA", "w", 2, "logical",
  "cutoffWGCNA", "x", 2, "numeric",
  "DElist", "i", 2, "character"), byrow=TRUE, ncol=4)


opt = getopt(spec)





# Help
if(!is.null(opt$help)) {
    cat("\nFlags:\n-d --data: excelsheet with normalized counts, rows as features (genes, miRNAs, proteins, N-features) and columns as samples.\n-m --metadata: excelsheet with metadata, minimum two columns named 'ids' (sample names matching those in the object above) and 'group' (diagnosis, tumor stage, ect.).\n(I) If the data comes from experimental batches and you want to correct for this, a column named 'batch' specifying which batch each sample belongs to (A,B,C,D, time1, time2, time3 ect) should also be included in the metadata. N.B specifying batches by numbers alone is not allowed.\n(II) If you want to perform correlation analysis a column named 'serum' must be included in the metadata specifying (in a binary way) which samples have a matched serum samples (1) and which that do not (0). N.B. if paired samples are used the column 'serum' should only have the value 1 for those samples (either tumours or normals, A or B ect.) you choose to test for - not both.\n(III) If you are interested in performing survival analysis a column named 'survival' must be included specifying (in a binary way) which samples have survival information (1) and which do not (0). N.B. if you have paired cancer and normal samples the column 'survival' should only have the value 1/0 for tumour samples (NA or other character values should be used for normal samples.\n(IV) If you want to include covariates in your survival analysis these should be included in the metadata sheet as a column(s).\n-s --serumdata: excelsheet of normalized counts from serum with the same order and as the count matrix (option -d).\n-u --survival: survival info must be included in the metadata excel sheet. The metadata file must contain at least four columns named; 'ids'(sample identifiers), 'age' (age in years at diagnosis, surgery or entry into trail), 'outcome.time' (time until end of follow-up in weeks, months or years, censuring, death) and 'outcome' (numeric 0 = censuring, 1=death). N.B. if you have (paired) normal samples the columns with survival information for these samples should contain NA values.\n-q --survplot: Arguments which specifies number of features to include per survival plot, e.g. many features requireds splitting of the plot, default features per plot is 50.\n-t --transform: should data be transformed? Current options are 'log2', 'log10' 'logit' or 'voom'. If argument is left out, no transformation of data will occur.\n-b --databatch: TRUE or FALSE specifies if you want to correct for experimental sample (tissue/interstitial fluid) batches. Batch information should be included in the metadata in a column named 'batch'.\n-e --serumbatch: TRUE or FALSE specifies if you want to correct for experimental serum sample batches. Batch information should be included in the metadata in a column named 'sbatch'.\n-n --filename: Name of result files from analysis.\n-f --logFC: Log fold change cut-off defining significant hits (proteins, genes, miRNAs or N-features). If omitted logFC cutoff will be set to 1.\n-r --FDR: Alpha level for significance.\n-o --plotmds: TRUE or FALSE specifies if a preliminary MDSplot should be made for data overview.\n-p --survcovar: Covariates to include in survival analysis. If multiple of these, they should be specified with commas as separator (e.g. Covar1,Covar2,Covar3), names should match the desired columns in the metadata sheet.\n-y --stratify: This flag may be used if some of the categorical (NOT continious) covariates violate the cox proportional assumption. The pipline checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the pipeline with this flag followed by the names of the categorical covariates which failed and these will be stratified.\n-c --colors: Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) must be separted by commas, example: green,red,blue. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.\n-l --lasso: Argument specifying in a logical way if LASSO regression should be performed. Default setting is FALSE.\n-w --WGCNA: Argument specifying in a logical way if Weighed Gene Co-expression Network Analysis should be performed. Default setting is FALSE. \n-x --cutoffWGCNA: Argument specifying the cutoff value, in %, for top most interconnected genes (or other features) from each modules identified in the Weighed Gene Co-expression Network Analysis.\n-i --DElist: Personal list of DE targets, one ID (name) per line. IDs (names) much match at least one of those found in the count data rows.\n-a --plotheatmap: TRUE or FALSE specifies if heatmap of DE/DA features should be made.\n\n")
    
        stop("Argument -h (help) selected, exiting script.")
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                         ### LOAD PACKAGES AND SET ARGUMENTS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD PACKAGES
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


list.of.packages <- c("limma", "sva", "glmnet", "openxlsx", "xlsx", "ggplot2", "heatmap.plus", "plyr", "data.table", "viridis", "squash", "survcomp", "survminer", "scales", "rms", "stackoverflow", "WGCNA", "fitdistrplus", "impute")

lapply(list.of.packages, library, character.only=T)

source("CAMPPFunctions.R")




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ARGUMENTS SPECIFYING DATASETS
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Data
if (is.null(opt$data)){
    stop("Argument data (-d) is missing. Data provided must be an excel sheet with features (genes, proteins ect.) as rows and samples as columns.")
} else {
    arg.data <- opt$data
    arg.data <- openxlsx::read.xlsx(arg.data, colNames = TRUE, rowNames = TRUE)
    data.names <- rownames(arg.data)
    arg.data <- as.matrix(as.data.frame(lapply(arg.data, as.numeric)))
    rownames(arg.data) <- data.names
}


# Serumdata
if (is.null(opt$serumdata)){
    arg.serumdata <- NULL
} else {
    arg.serumdata <- opt$serumdata
    arg.serumdata <- openxlsx::read.xlsx(arg.serumdata, colNames = TRUE, rowNames = TRUE)
    serumdata.names <- rownames(arg.serumdata)
    arg.serumdata <- as.matrix(as.data.frame(lapply(arg.serumdata, as.numeric)))
    rownames(arg.serumdata) <- serumdata.names
}



# Missing Values (NAs)

if (TRUE %in% unique(as.vector(is.na(arg.data)))) {
    arg.data <- ReplaceNAs(arg.data)
    
}

if (!is.null(arg.serumdata)) {
    if (TRUE %in% unique(as.vector(is.na(arg.serumdata)))) {
        arg.serumdata <- ReplaceNAs(arg.serumdata)
    }
}



# Metadata
if (is.null(opt$metadata)){
    stop("Argument metadata (-m) is missing. Metadata provided must be an excel sheet with minimum two columns named 'ids' (sample names matching those in the object above) and 'group' (diagnosis, tumor stage, ect.).")
} else {
    arg.metadata <- opt$metadata
    arg.metadata <- openxlsx::read.xlsx(arg.metadata, colNames = TRUE, rowNames = FALSE)
    
    # Match data and metadata
    arg.metadata <- arg.metadata[arg.metadata$ids %in% colnames(arg.data),]

    # Vector with group and batch for DE and MDS
    try(arg.group <- as.factor(as.character(arg.metadata$group)))
    if (length(arg.group)<1) {
        stop("Column 'group' is missing from the metadata. This column is mandatory, errors will arise. See -h for help.")
    }
    try(arg.batch <- as.factor(as.character(arg.metadata$batch)))
    try(arg.sbatch <- as.factor(as.character(arg.metadata$sbatch)))
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### ADDITIONAL ARGUMENTS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Data transformation
#if (is.null(opt$transform)) {
#    print("Argument transform (-t) is missing. Data will not be transformed.")
#} else if (opt$transform == "log2") {
#    arg.data <- log2(arg.data)
#    if(!is.null(arg.serumdata)) {
#        arg.serumdata <- log2(arg.serumdata)
#    }
#} else if (opt$transform == "logit") {
#    arg.data <- logit(arg.data)
#    if(!is.null(arg.serumdata)) {
#        arg.serumdata <- logit(arg.serumdata)
#    }
#} else if (opt$transform == "log10"){
#    arg.data <- log10(arg.data)
#   if(!is.null(arg.serumdata)) {
#        arg.serumdata <- log10(arg.serumdata)
#    }
#} else if (opt$transform == "voom"){
#    design_voom <- model.matrix(~0+arg.group)
#    arg.data <- voom(arg.data, design_voom, plot=FALSE)
#    if(!is.null(arg.serumdata)) {
#        arg.serumdata <- voom(arg.serumdata, design_voom, plot=FALSE)
#    }
#} else {
#    print("Argument t is abused. Data will not be transformed.")
#}



# Data transformation
if (is.null(opt$transform)) {
    print("Argument transform (-t) is missing. Data will not be transformed.")
} else {
    hasZero <- unique(as.vector(arg.data == 0))
    if (TRUE %in% hasZero) {
        print("Log tranformation cannot be performed with values == 0, these will be replaced by either variable minimum or mean - See user manual.")
        if (opt$transform == "log2") {
            arg.data <- log2(ReplaceZero(arg.data))
        } else if (opt$transform == "logit") {
            arg.data <- logit(ReplaceZero(arg.data))
        } else if (opt$transform == "log10"){
            arg.data <- log10(ReplaceZero(arg.data))
        } else if (opt$transform == "voom"){
            design_voom <- model.matrix(~0+arg.group)
            arg.data <- voom(arg.data, design_voom, plot=FALSE)
        } else {
            print("Argument t is abused. Data will not be transformed.")
        }
    } else {
        if (opt$transform == "log2") {
            arg.data <- log2(arg.data)
        } else if (opt$transform == "logit") {
            arg.data <- logit(arg.data)
        } else if (opt$transform == "log10"){
            arg.data <- log10(arg.data)
        } else if (opt$transform == "voom"){
            design_voom <- model.matrix(~0+arg.group)
            arg.data <- voom(arg.data, design_voom, plot=FALSE)
        } else {
            print("Argument t is abused. Data will not be transformed.")
        }
    }
}



if (!is.null(arg.serumdata) & !is.null(opt$transform)) {
    hasZero <- unique(as.vector(arg.serumdata == 0))
    if (TRUE %in% hasZero) {
        print("Log tranformation cannot be performed with values == 0, these will be replaced by either variable minimum or mean - See user manual.")
        if (opt$transform == "log2") {
            arg.serumdata <- log2(ReplaceZero(arg.serumdata))
        } else if (opt$transform == "logit") {
            arg.serumdata <- logit(ReplaceZero(arg.serumdata))
        } else if (opt$transform == "log10"){
            arg.serumdata <- log10(ReplaceZero(arg.serumdata))
        } else if (opt$transform == "voom"){
            design_voom <- model.matrix(~0+arg.serumdata)
            arg.serumdata <- voom(arg.serumdata, design_voom, plot=FALSE)
        } else {
            print("Argument t is abused. Data will not be transformed.")
        }
    } else {
        if (opt$transform == "log2") {
            arg.serumdata <- log2(arg.serumdata)
        } else if (opt$transform == "logit") {
            arg.serumdata <- logit(arg.serumdata)
        } else if (opt$transform == "log10"){
            arg.serumdata <- log10(arg.serumdata)
        } else if (opt$transform == "voom"){
            design_voom <- model.matrix(~0+arg.group)
            arg.serumdata <- voom(arg.serumdata  , design_voom, plot=FALSE)
        } else {
            print("Argument t is abused. Data will not be transformed.")
        }
    }
}




# Data batch correction
if (is.null(opt$databatch)){
    arg.databatch <- FALSE
} else {
    arg.databatch <- opt$databatch
}


# Serum data batch correction
if (is.null(opt$serumbatch)){
    arg.serumbatch <- FALSE
} else {
    arg.serumbatch <- opt$serumbatch
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


# LogFC
if (is.null(opt$logFC)){
   print("No log fold change cut-off for significant hits has been chosen. For log2 transformed data a logFC cut-off of +/- 1 is most commenly used. For for untransformed data a fold change of +/- 2 is equivalent")
   arg.logFC <- 0
} else {
    arg.logFC <- opt$logFC
}


# FDR
if (is.null(opt$FDR)){
    arg.FDR <- 0.05
} else {
    arg.FDR <- opt$logFDR
}


# Colors
if (is.null(opt$colors)){
    arg.colors <- as.vector(viridis(length(levels(as.factor(as.character(arg.metadata$group)))), begin = 0.2, end = 0.8, option="cividis"))
    #if (length(levels(as.factor(as.character(arg.metadata$group)))) < 3) {
    #    arg.colors <- arg.colors[1:2]
    #}
} else {
    arg.colors <- as.list(strsplit(opt$colors, ",")[[1]])
}


# Filename
if (is.null(opt$filename)){
    arg.filename <- "my_results"
} else {
    arg.filename <- opt$filename
}



# Custum DElist
if (is.null(opt$DElist)){
    arg.DElist <- NULL
} else {
    arg.DElist <- openxlsx::read.xlsx(opt$DElist, colNames = FALSE)$X1
}


# Heatmap
if (is.null(opt$plotheatmap)){
    arg.plotheatmap <- FALSE
} else {
    arg.plotheatmap <- opt$plotheatmap
}


# LASSO
if (is.null(opt$lasso)){
    arg.lasso <- FALSE
} else {
    arg.lasso <- opt$lasso
}



# WGCNA
if (is.null(opt$WGCNA)){
    arg.WGCNA <- FALSE
} else {
    arg.WGCNA <- opt$WGCNA
}


# CutoffWGCNA
if (is.null(opt$cutoffWGCNA)){
    arg.cutoffWGCNA <- 25
} else {
    arg.cutoffWGCNA <- opt$cutoffWGCNA
}



# Survival Analysis
if (is.null(opt$survival)){
    arg.survival <- FALSE
} else {
    arg.survival <- opt$survival
}



# Survival covariates
if (is.null(opt$survcovar)){
    arg.survcovar <- NULL
} else {
    arg.survcovar <- unlist(strsplit(opt$survcovar, split=","))
}



# Stratify
if (is.null(opt$stratify)){
    arg.stratify <- NULL
} else {
    arg.stratify <- unlist(strsplit(opt$stratify, split=","))
}



# Survival Plot
if (is.null(opt$survplot)){
    arg.survplot <- 50
} else {
    arg.survplot <- opt$survplot
}




# Create directory Results
dir.create("Results")
setwd("Results/")


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### BATCH CORRECTION ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




if (arg.databatch == TRUE){
    if (length(arg.batch) > 0) {
        arg.design <-  model.matrix(~arg.group)
        data.batch <- ComBat(as.matrix(arg.data), arg.batch, arg.design, par.prior=TRUE,prior.plots=FALSE)
    } else {
        data.batch <- arg.data
        print("No column named batch in the metadata file. Continuing without batch correction.")
    }
} else {
    print("No batch correction requested")
}





if (arg.serumbatch == TRUE){
    if (length(arg.sbatch) > 0) {
        arg.design <- model.matrix(~arg.group)
        serum.data.batch <- ComBat(as.matrix(arg.serumdata), arg.sbatch, arg.design, par.prior=TRUE,prior.plots=FALSE)
    } else {
        serum.data.batch <- arg.serumdata
        print("No column named sbatch in the metadata file. Continuing without serum batch correction.")
    }
} else {
    print("No serum batch correction requested")
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                        ### Distributional Checks ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (arg.datacheck == TRUE) {
    if (arg.databatch == TRUE) {
        subset.data <- data.batch[sample(nrow(data.batch), 6),]
    } else {
        subset.data <- arg.data[sample(nrow(arg.data), 6),]
    }
    
    list.of.dists <- fit_distributions(subset.data)
    
    dir.create("DataChecks")
    setwd("DataChecks/")
    plot_distributions(subset.data, list.of.dists)
    setwd("..")
}




if (arg.datacheck == TRUE & !is.null(arg.serumdata)) {
    if (arg.serumbatch == TRUE) {
        subset.data <- serum.data.batch[sample(nrow(serum.data.batch), 6),]
    } else {
        subset.data <- arg.serumdata[sample(nrow(arg.serumdata), 6),]
    }
    
    list.of.dists <- fit_distributions(subset.data)
    
    dir.create("SerumdataChecks")
    setwd("SerumdataChecks/")
    plot_distributions(subset.data, list.of.dists)
    setwd("..")
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### PRELIMINARY MDS PLOT ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



if (arg.plotmds == TRUE && arg.databatch == TRUE){
    mdsplot <- myMDSplot(data.batch, arg.group, "", arg.colors[1:length(levels(as.factor(as.character(arg.metadata$group))))])
    ggsave(paste0(arg.filename, "_MDSplot_batchcorr.pdf"), bg = "transparent", plot = mdsplot)

} else if (arg.plotmds == TRUE && arg.databatch == FALSE){
    mdsplot <- myMDSplot(arg.data, arg.group, "", arg.colors[1:length(levels(as.factor(as.character(arg.metadata$group))))])
    ggsave(paste0(arg.filename, "_MDSplot.pdf"), bg = "transparent", plot = mdsplot)
    
} else {
    print("No preliminary plot requested.")
}










# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### DIFFERENTIAL EXPRESSION ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# Differential Expression Analysis with Limma


res.DE <- NULL


if (is.null(arg.DElist)){
    
    combinations <- data.frame(t(combn(paste0("arg.group", levels(arg.group)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    
    if (arg.databatch == "FALSE") {
        # Design incorporating group
        arg.design <-  model.matrix(~0+arg.group)
        
        # Making group contrasts
        contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(arg.design))))
        
        # Apply DE_limma function to all comparisons
        res.DE <- DA_feature_apply(contrast.matrix, arg.data, arg.design, arg.logFC, arg.FDR, NULL, FALSE)
        
        # Write results out as excel file
        DE.out <- excel_output(res.DE, paste0(arg.filename,"_DE"))
    
    } else {
        
        if (length(arg.batch) > 0) {
            
            # Design incorporating group and batch
            arg.design <-  model.matrix(~0+arg.group+arg.batch)
            
            # Making group contrasts
            contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(arg.design))))
            
            # Apply DE_limma function to all comparisons
            res.DE <- DA_feature_apply(contrast.matrix, arg.data, arg.design, arg.logFC, arg.FDR, NULL, FALSE)
            
            # Write results out as excel file
            DE.out <- excel_output(res.DE, paste0(arg.filename,"_databatch_DE"))
        }
    }
} else {
    print("You have provided a custom list of DE/DA features. Differential expression/abundance analysis with limma will be skipped.")

}










# DE feature names

if (!is.null(res.DE)) {
    res.DE.names <- do.call(rbind, unlist(res.DE, recursive=FALSE))$name
} else if (!is.null(arg.DElist)) {
    res.DE.names <- arg.DElist
} else {
    stop("No signficant DE/DA hits found. Check output file from differential expression analysis. You can also provide a custom file of DE/DA features. See argument -l.")
}






# Plotting results of differential expression analysis

if (arg.plotheatmap == TRUE) {
    # color scheme
    arg.colors.hm <- get_colors(as.character(arg.group), arg.colors[1:length(levels(as.factor(as.character(arg.metadata$group))))])
    
    if (arg.databatch == TRUE){
        arg.DE.hm <- data.batch[rownames(data.batch) %in% res.DE.names,]
        name <- "_batchcorr"
    } else {
        arg.DE.hm <- arg.data[rownames(arg.data) %in% res.DE.names,]
        name <- ""
    }
    
    # heatmap colors in blue
    arg.hm.gradient <- viridis(300)
    
    # Heatmap as pdf
    my_heatmap(arg.DE.hm, arg.hm.gradient, arg.colors.hm, arg.group, arg.filename)

} else {
    print("No heatmap requested.")
}






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                                       ## LASSO Regression ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



if (arg.lasso == TRUE) {
    
    group.LASSO <- as.integer(arg.group)
    seeds <- sample(1:1000, 10)
    LASSO.res <- list()
    
    if (arg.databatch == TRUE) {
        if (length(levels(as.factor(group.LASSO))) > 2) {
            for (idx in 1:length(seeds)) {
                LR <- LASSO_feature(seeds[[idx]], data.batch, group.LASSO, TRUE)
                LASSO.res[[idx]] <-  LR
            }
        } else {
            for (idx in 1:length(seeds)) {
                LR <- LASSO_feature(seeds[[idx]], data.batch, group.LASSO, FALSE)
                LASSO.res[[idx]] <-  LR
            }
        }
    } else {
        if (length(levels(as.factor(group.LASSO))) > 2) {
            for (idx in 1:length(seeds)) {
                LR <- LASSO_feature(seeds[[idx]], arg.data, group.LASSO, TRUE)
                LASSO.res[[idx]] <-  LR
            }
        } else {
            for (idx in 1:length(seeds)) {
                LR <- LASSO_feature(seeds[[idx]], arg.data, group.LASSO, FALSE)
                LASSO.res[[idx]] <-  LR
            }
        }
    }
    LASSO.res1 <- Reduce(intersect, lapply(LASSO.res, '[[', 1))
    LASSO.res1 <- data.frame(LASSO.res1[-1])
    colnames(LASSO.res1) <- c("LASSO.Var.Select")
    xlsx::write.xlsx(LASSO.res1, file=paste0(arg.filename,"_LASSO.xlsx"), row.names=FALSE)
    consensus <- DE.out[DE.out$name %in% LASSO.res1$LASSO.Var.Select,]
    xlsx::write.xlsx(consensus, file=paste0(arg.filename,"_DEA_LASSO_Consensus.xlsx"), row.names=FALSE)
    
    LASSO.res2 <- round(unlist(lapply(LASSO.res, '[[', 2)), digits = 4)
    cat(paste0("Mean cross validation error (cv.glmnet) = ", LASSO.res2,".\n"))
}







# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### SERUM-TIF CORRELATION ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Data (TIF) and Serumdata correlations

if (!is.null(arg.serumdata)) {
    # Extract sample ID IF sample has matched serum.
    arg.samples <- arg.metadata[arg.metadata$serum == 1,]$ids
    
    if (arg.databatch == TRUE) {
        data.corr <- data.batch[,colnames(data.batch) %in% arg.samples]
        if (arg.serumbatch == TRUE) {
            serum.corr <- serum.data.batch
        }  else {
            serum.corr <- arg.serumdata
        }
    } else {
        data.corr <- arg.data[,colnames(arg.data) %in% arg.samples]
        if (arg.serumbatch == TRUE) {
            serum.corr <- serum.data.batch
        }  else {
            serum.corr <- arg.serumdata
        }
    }
    
    # Perform correction analysis and generate overall correlation plot
    res.corr <- my_correlation(data.corr, serum.corr, res.DE.names, arg.filename)
    
    
    # Print out significant hits in Excel
    res.corr$sig.corr <- ifelse(res.corr$fdr <= arg.FDR, "yes", "no")
    xlsx::write.xlsx(res.corr, file=paste0(arg.filename,"_corr_serum.xlsx"), row.names=FALSE)


    # Individual correlation plots
    #corr.features <- as.character(res.corr[res.corr$sig.corr == "yes",]$name)
    corr.features <- as.character(res.corr[which(res.corr$sig.corr == "yes"),]$name)
    if (length(corr.features) > 1) {
        my_correlation_plots(data.corr, serum.corr, corr.features, arg.filename)
    } else {
        print("No significant correlations, no plotting.")
    }
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### SURVIVAL ANALYSIS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






                                                                        ### SETTING UP DATA ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (arg.survival == TRUE){
   
   # Setting up dataset for survival analysis, matching TIF and serum samples either batch or no batch. Only DA/DE featrues are used.
   metadata.surv <- arg.metadata[which(arg.metadata$survival == 1),]
   arg.samples <- metadata.surv$ids

   if (arg.databatch == TRUE) {
       data.surv <- data.batch[,colnames(data.batch) %in% arg.samples]
   } else {
       data.surv <- arg.data[,colnames(arg.data) %in% arg.samples]
   }

   data.surv <- data.surv[rownames(data.surv) %in% res.DE.names,]
   
   
   
   # If user has specified covariates for survivalanalysis, extract these.
   if(is.null(arg.survcovar)) {
       surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$age, metadata.surv$outcome)
       colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "age", "outcome")
   } else {
       my.covars <- metadata.surv[,colnames(metadata.surv) %in% arg.survcovar]
       surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$age, metadata.surv$outcome, my.covars)
       colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "age", "outcome", arg.survcovar)
   }
   
   
   # datadist format for cph function with cubic splines.
   features <- rownames(data.surv)
   dd <- datadist(surv_object); options(datadist='dd')


   
   
                                                                        ### User Specified Covariates ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
   
   
   
   # Evaluating class of covariate, continious or categorical.
   arg.survcovar.original <- arg.survcovar
   
   if(!is.null(arg.survcovar)) {
       for (i in 1:length(arg.survcovar)) {
           if(class(eval(parse(text=paste0("surv_object$", arg.survcovar[i])))) %in% c("integer", "numeric")) {
               arg.survcovar[i] <- paste0("rcs(", arg.survcovar[i], ")")
           }
       }
       
       arg.survcovar <- paste(arg.survcovar,collapse="+")
   }
   
   
   
   
                                                                        ### Liniarity of Continious Covariates ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
   
   
   
   # Survival models with or without user specified covariates.
   covariate_linearity <- list()
   
   if(!is.null(arg.survcovar)) {
       for (f in features) {
           acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(age) + rcs(", as.character(f),") + ", arg.survcovar, ", data = surv_object, x=TRUE,y=TRUE)"))
           eval(acall)
           covariate_linearity[[as.character(f)]] <- result
       }
   } else {
       for (f in features) {
           acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(age) + rcs(", as.character(f),"), data = surv_object, x=TRUE,y=TRUE)"))
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


   # Re-assign features
   features <- colnames(surv_object)
   

   # Anova comparing linear covariate and non-linear covariate.
   covariate_linearity  <- lapply(covariate_linearity, function(x) anova(x))
   
   # Unlisting and extracting pvalues for non-linear fit.
   covariate_linearity <- data.frame(do.call(rbind, lapply(covariate_linearity, function(x) x[grep("Nonlinear|NONLINEAR", rownames(x)),3])))
   
   # FDR correction for multiple testing
   covariate_linearity <- apply(covariate_linearity, 2, function(x) p.adjust(x, method = "fdr", n=nrow(covariate_linearity)))
   
   
   
   
   
   # Colnames with and without user-specified covariates.
   if(is.null(arg.survcovar)) {
       colnames(covariate_linearity) <- c("Age", "Feature", "Total")
   } else {
       covar.names <- unlist(strsplit(arg.survcovar, split="[+]"))
       if(length(grep("rcs", covar.names)) > 0) {
           covar.names <- c("Age", "Feature", covar.names[grep("rcs", covar.names)] ,"Total")
           covar.names <- gsub("rcs\\(|\\)", "", covar.names)
           colnames(covariate_linearity) <- covar.names
       } else {
           colnames(covariate_linearity) <- c("Age", "Feature", "Total")
       }
   }
   
   
   # Extracting p-values and filtering for significance.
   covariate_nonlinear <- colnames(covariate_linearity)[apply(covariate_linearity, 2, function(x) any(x < 0.05))]
   
   
   # Updating covariates with cubic splines.
   if (length(covariate_nonlinear) > 0) {
       cat(paste0("The following continious covariate(s) may be violating the assumption of linearity:  ", covariate_nonlinear,".\nCubic splines will be added.\n"))
       for (i in 1:length(arg.survcovar.original)) {
           if (arg.survcovar.original[i] %in% covariate_nonlinear) {
               arg.survcovar.original[i] <- paste0("rcs(", arg.survcovar.original[i], ")")
           }
       }
   }

   






                                                                ### Testing Cox Proportional Hazard Assumption ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



   # Picking model based on result of linearity check.
   pha_check <- list()

   
   for (f in features) {
       if("age" %in% covariate_nonlinear && "feature" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)){
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(age) + rcs(", as.character(f), "), data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(age) + rcs(", as.character(f),") + ", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           }
       }
       else if ("feature" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ age + rcs(", as.character(f), "), data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ age + rcs(", as.character(f),") + ", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           }
       }
       else if ("age" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(age) + ", as.character(f), ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(age) + ", as.character(f), " + ", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           }
       }
       else {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ age +", as.character(f), ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ age +", as.character(f), "+ ", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
               eval(acall)
               pha_check[[as.character(f)]] <- result 
           }
       }
   }
   
   



   # Filtering p-values for significance.
   pha.fail.test <- as.character(unlist(lapply(pha_check, function(x) rownames(x)[apply(x, 1, function(u) any(u < 0.05))])))


   cat("WARNING: The following features and/or covariates failed the test of proportional hazard: ", pha.fail.test, "\nIF the covariates that failed are categorical you may use strata by re-running the pipline adding flag -y followed by the names of the categorical covariates to stratify (if multiple then separate by comma). \nN.B, this pipeline does not handle continuous variables that violate the proportional hazard assumption, if any of these failed PH test, the hazard ratios of these should NOT be evaluated.\n")
   
   
   
   # User stratification of categorical covariates which violate the proportional hazard assumption.
   if(!is.null(arg.stratify)) {
       for (i in 1:length(arg.survcovar.original)) {
           if (arg.survcovar.original[i] %in% arg.stratify) {
               arg.survcovar.original[i] <- paste0("strat(", arg.survcovar.original[i], ")")
           }
       }
   }


   





                                                                ### Survival Analysis with updated covariates ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




   # Survival analysis, cox-regression
   
   # Picking model based on result of linearity check and proportional hazard text.
   survival.results <- list()
   
   
   for (f in features) {
       if("age" %in% covariate_nonlinear && "feature" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)){
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + rcs(age), data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + rcs(age) +", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           }
       }
       else if ("feature" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + age, data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + age +", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           }
       }
       else if ("age" %in% covariate_nonlinear) {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), "+ rcs(age), data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), "+ rcs(age) +", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           }
       }
       else {
           if(is.null(arg.survcovar.original)) {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), "+ age, data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           } else {
               acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), "+ age +", arg.survcovar.original, ", data = surv_object, x=TRUE,y=TRUE)"))
               eval(acall)
               survival.results[[as.character(f)]] <- result
           }
       }
   }
   
   
   survival.results <- survival.results[-c((length(survival.results)-2):length(survival.results))]

   # Setting up data and writing out excel sheet with results and making HR stemplot
   survival.data <- my_survival(survival.results, arg.filename, arg.survplot)
   xlsx::write.xlsx(survival.data, file=paste0(arg.filename,"_survival.xlsx"), row.names=TRUE)
   

}


sink()



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                               ## Weighed Gene Co-Expression Network Analysis ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


if (arg.WGCNA == TRUE) {
    
    if (arg.databatch == TRUE){
        data.WGCNA <- t(data.batch)
    } else {
        data.WGCNA <- t(arg.data)
    }
    
    # Check data
    gsg <- goodSamplesGenes(data.WGCNA)
    cat(paste0("Data set is OK for WGCNA - ", gsg$allOK,".\n"))
    
    # Plot power estimates of data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    powers <-  c(c(1:10), seq(from = 12, to=20, by=2));
    sft <- pickSoftThreshold(data.WGCNA, dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")
    
    pdf(paste0(arg.filename,"_WGCNA_softpowerplot.pdf"))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n")
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],col="red")
    dev.off()
    
    # Set power to best softpower estimate
    softPower <- sft$powerEstimate
    softPower
    
    # Construct adjecancy and topology matrix
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    # Adjacency matrix and Topology Matrix
    adj <-adjacency(data.WGCNA, power = softPower);
    
    TOM <- TOMsimilarity(adj)
    colnames(TOM) <- colnames(data.WGCNA)
    rownames(TOM) <- colnames(data.WGCNA)
    
    dissTOM <- (1-TOM)
    
    # Dendogram
    geneTree <- hclust(as.dist(dissTOM),method="ward.D2")
    
    
    # Construct Modules
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    minModuleSize <- 10
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
    table(dynamicMods)
    
    # Extract module color labels
    dynamicColors <- labels2colors(dynamicMods)
    table(dynamicColors)
    
    # Eigen features in each module
    MEList <- moduleEigengenes(data.WGCNA, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- (1-cor(MEs))
    
    # Cluster module eigengenes
    METree <- hclust(dist(MEDiss), method = "ward.D2")

    MEDissThres <- 0.25
    merge <- mergeCloseModules(data.WGCNA, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    mergedColors <- merge$colors
    mergedMEs <- merge$newMEs
    
    # Plot merged modules
    pdf(paste0(arg.filename,"_WGCNA_ModuleTree.pdf"))
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()
    
    moduleColors <- mergedColors
    colorOrder <- c("grey", standardColors(50))
    moduleLabels <- match(moduleColors, colorOrder)-1
    
    # Interconnectivity of features in modules
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    IC <- intramodularConnectivity(adj,  moduleColors)
    mod.cols <- levels(as.factor(moduleColors))
    
    WGCNAres <- ModuleIC(mod.cols, moduleColors, IC, data.WGCNA, arg.cutoffWGCNA, arg.filename)
    WGCNAres <- do.call("rbind", WGCNAres)
    rownames(WGCNAres) <- gsub(".*[.]", "", rownames(WGCNAres))
    xlsx::write.xlsx(WGCNAres, file=paste0(arg.filename,"_WGCNAres.xlsx"), row.names=TRUE)
}





