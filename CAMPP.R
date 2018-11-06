
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
  "DElist", "l", 2, "character"), byrow=TRUE, ncol=4)


opt = getopt(spec)





# Help
if(!is.null(opt$help)) {
    cat("\nFlags:\n-d --data: excelsheet with normalized counts, rows as features (genes, miRNAs, proteins, N-features) and columns as samples.\n-m --metadata: excelsheet with metadata, minimum two columns named 'ids' (sample names matching those in the object above) and 'group' (diagnosis, tumor stage, ect.).\n(I) If the data comes from experimental batches and you want to correct for this, a column named 'batch' specifying which batch each sample belongs to (A,B,C,D, time1, time2, time3 ect) should also be included in the metadata. N.B specifying batches by numbers alone is not allowed.\n(II) If you want to perform correlation analysis a column named 'serum' must be included in the metadata specifying (in a binary way) which samples have a matched serum samples (1) and which that do not (0). N.B. if paired samples are used the column 'serum' should only have the value 1 for those samples (either tumours or normals, A or B ect.) you choose to test for - not both.\n(III) If you are interested in performing survival analysis a column named 'survival' must be included specifying (in a binary way) which samples have survival information (1) and which do not (0). N.B. if you have paired cancer and normal samples the column 'survival' should only have the value 1/0 for tumour samples (NA or other character values should be used for normal samples.\n(IV) If you want to include covariates in your survival analysis these should be included in the metadata sheet as a column(s).\n-s --serumdata: excelsheet of normalized counts from serum with the same order and as the count matrix (option -d).\n-u --survival: survival info must be included in the metadata excel sheet. The metadata file must contain at least four columns named; 'ids'(sample identifiers), 'age' (age in years at diagnosis, surgery or entry into trail), 'outcome.time' (time until end of follow-up in weeks, months or years, censuring, death) and 'outcome' (numeric 0 = censuring, 1=death). N.B. if you have (paired) normal samples the columns with survival information for these samples should contain NA values.\n-q --survplot: Arguments which specifies number of features to include per survival plot, e.g. many features requireds splitting of the plot, default features per plot is 50.\n-t --transform: should data be transformed? Current options are 'log2', 'log10' 'logit' or 'voom'. If argument is left out, no transformation of data will occur.\n-b --databatch: TRUE or FALSE specifies if you want to correct for experimental sample (tissue/interstitial fluid) batches. Batch information should be included in the metadata in a column named 'batch'.\n-e --serumbatch: TRUE or FALSE specifies if you want to correct for experimental serum sample batches. Batch information should be included in the metadata in a column named 'sbatch'.\n-n --filename: Name of result files from analysis.\n-f --logFC: Log fold change cut-off defining significant hits (proteins, genes, miRNAs or N-features). If omitted logFC cutoff will be set to 1.\n-r --FDR: Alpha level for significance.\n-o --plotmds: TRUE or FALSE specifies if a preliminary MDSplot should be made for data overview.\n-p --survcovar: Covariates to include in survival analysis. If multiple of these, they should be specified with commas as separator (e.g. Covar1,Covar2,Covar3), names should match the desired columns in the metadata sheet.\n-y --stratify: This flag may be used if some of the categorical (NOT continious) covariates violate the cox proportional assumption. The pipline checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the pipeline with this flag followed by the names of the categorical covariates which failed and these will be stratified.\n-c --colors: Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) must be separted by commas, example: green,red,blue. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.\n-l --DElist: Personal list of DE targets, one ID (name) per line. IDs (names) much match at least one of those found in the count data rows.\n-a --plotheatmap: TRUE or FALSE specifies if heatmap of DE/DA features should be made.\n\n")
    
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


list.of.packages <- c("limma", "sva", "openxlsx", "xlsx", "ggplot2", "heatmap.plus", "plyr", "data.table", "RColorBrewer", "squash", "survcomp", "survminer", "scales", "rms", "stackoverflow")

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



# Metadata
if (is.null(opt$metadata)){
    stop("Argument metadata (-m) is missing. Metadata provided must be an excel sheet with minimum two columns named 'ids' (sample names matching those in the object above) and 'group' (diagnosis, tumor stage, ect.).")
} else {
    arg.metadata <- opt$metadata
    arg.metadata <- openxlsx::read.xlsx(arg.metadata, colNames = TRUE, rowNames = FALSE)
    # Vector with group and batch for DE and MDS
    try(arg.group <- as.factor(as.character(arg.metadata$group)))
    if (length(arg.group)<1) {
        stop("Column 'group' is missing from the metadata. This column is mandatory, errors will arise. See -h for help.")
    }
    try(arg.batch <- as.factor(as.character(arg.metadata$batch)))
    try(arg.sbatch <- as.factor(as.character(arg.metadata$sbatch)))
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






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### ADDITIONAL ARGUMENTS ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Data transformation
if (is.null(opt$transform)) {
    print("Argument transform (-t) is missing. Data will not be transformed.")
} else if (opt$transform == "log2") {
    arg.data <- log2(arg.data)
    if(!is.null(arg.serumdata)) {
        arg.serumdata <- log2(arg.serumdata)
    }
} else if (opt$transform == "logit") {
    arg.data <- logit(arg.data)
    if(!is.null(arg.serumdata)) {
        arg.serumdata <- logit(arg.serumdata)
    }
} else if (opt$transform == "log10"){
    arg.data <- log10(arg.data)
    if(!is.null(arg.serumdata)) {
        arg.serumdata <- log10(arg.serumdata)
    }
} else if (opt$transform == "voom"){
    design_voom <- model.matrix(~0+arg.group)
    arg.data <- voom(arg.data, design_voom, plot=FALSE)
    if(!is.null(arg.serumdata)) {
        arg.serumdata <- voom(arg.serumdata, design_voom, plot=FALSE)
    }
} else {
    print("Argument t is abused. Data will not be transformed.")
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
    arg.colors <- as.vector(brewer.pal(length(levels(as.factor(as.character(arg.metadata$group)))),"Dark2"))
    if (length(levels(as.factor(as.character(arg.metadata$group)))) < 3) {
        arg.colors <- arg.colors[1:2]
    }
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



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                        ### BATCH CORRECTION ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




if (arg.databatch == TRUE){
    if (length(arg.batch) > 0) {
        arg.design <-  model.matrix(~arg.group)
        data.batch <- ComBat(arg.data, arg.batch, arg.design, par.prior=TRUE,prior.plots=FALSE)
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
        serum.data.batch <- ComBat(arg.serumdata, arg.sbatch, arg.design, par.prior=TRUE,prior.plots=FALSE)
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
                                                                        ### PRELIMINARY MDS PLOT ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



if (arg.plotmds == TRUE && arg.databatch == TRUE){
    mdsplot <- myMDSplot(data.batch, arg.group, "", arg.colors[1:length(levels(as.factor(as.character(arg.metadata$group))))])
    ggsave(paste0(arg.filename, "_MDSplot_batchcorr.pdf"), plot = mdsplot)

} else if (arg.plotmds == TRUE && arg.databatch == FALSE){
    mdsplot <- myMDSplot(arg.data, arg.group, "", arg.colors[1:length(levels(as.factor(as.character(arg.metadata$group))))])
    ggsave(paste0(arg.filename, "_MDSplot.pdf"), plot = mdsplot)
    
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
        excel_output(res.DE, paste0(arg.filename,"_DE"))

    
    } else {
        
        if (length(arg.batch) > 0) {
            
            # Design incorporating group and batch
            arg.design <-  model.matrix(~0+arg.group+arg.batch)
            
            # Making group contrasts
            contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(arg.design))))
            
            # Apply DE_limma function to all comparisons
            res.DE <- DA_feature_apply(contrast.matrix, arg.data, arg.design, arg.logFC, arg.FDR, NULL, FALSE)
            
            # Write results out as excel file
            excel_output(res.DE, paste0(arg.filename,"_databatch_DE"))
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
    arg.hm.gradient <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)
    
    # Heatmap as pdf
    my_heatmap(arg.DE.hm, arg.hm.gradient, arg.colors.hm, arg.group, arg.filename)

} else {
    print("No heatmap requested.")
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
   


   # Setting up data and writing out excel sheet with results and making HR stemplot
   survival.data <- my_survival(survival.results, arg.survplot)
   xlsx::write.xlsx(survival.data, file=paste0(arg.filename,"_survival.xlsx"), row.names=TRUE)
   

}


sink()
