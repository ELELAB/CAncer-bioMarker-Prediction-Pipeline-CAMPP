#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




													### FUNCTIONS ###



													
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BASE FUNCTION FOR INSTALLATION OF R-PACKAGES
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




install.packages.auto <- function(x, y) {
    if(isTRUE(x %in% .packages(all.available=TRUE))) {
        eval(parse(text = paste0("require('",x,"')")))
    } else {
        eval(parse(text = paste0("install.packages('", x, "', repos = '", y, "', dependencies = TRUE)")))
    }
    if(isTRUE(x %in% .packages(all.available=TRUE))) {
        eval(parse(text = paste0("require('",x,"')")))
    } else {
        eval(parse(text = paste0("BiocManager::install('",x,"', dependencies = TRUE)")))
    }
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading in files
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

SplitList <- function(my.list) {
    my.list <- as.character(unlist(strsplit(my.list, split=",")))
    return(my.list)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading in files
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



ReadMyFile <- function(my.data, my.expr) {
    if(my.expr == TRUE) {
        file <- try(my.data <- openxlsx::read.xlsx(my.data, colNames = TRUE, rowNames = FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Data file is not .xlsx, trying .txt\n")
            file <- try(my.data <- read.delim(my.data, header = TRUE), silent = TRUE)
            if (class(file) == "try-error") {
                file <- try(my.data <- read.delim(my.data, header = TRUE, sep=";"), silent = TRUE)
                if (class(file) == "try-error") {
                    stop("\n- Data file must be .xlsx or .txt\n")
                }
            }
        }
        
        
        # Average duplicates and get IDs
        colnames(my.data)[1] <- "IDs"
        IDs <- my.data$IDs
        my.data$IDs <- NULL
        
        my.data <- as.data.frame(lapply(my.data, as.numeric))
        my.data$IDs <- IDs
        my.data <-  data.frame(data.table(my.data)[, lapply(.SD, mean), by=IDs])
        
        my.names <- my.data$IDs
        my.data$IDs <- NULL
        my.data <- as.matrix(my.data)
        rownames(my.data) <- gsub("-", "_", my.names)
        
    } else {
        file <- try(my.data <- openxlsx::read.xlsx(my.data, colNames = TRUE, rowNames = FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Metadata file is not .xlsx, trying .txt\n")
            file <- try(my.data <- read.delim(my.data, header = TRUE), silent = TRUE)
            if (class(file) == "try-error") {
                file <- try(my.data <- read.delim(my.data, header = TRUE, sep =";"), silent = TRUE)
                if (class(file) == "try-error"){
                    stop("\n- Metadata file must be .xlsx or .txt\n")
                }
            }
        }
    }
    rm(file)
    return(my.data)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DATA CHECKS
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking Data for NA Values.
# Takes as arguments;
    # my.data = a dataframe of expression/abundance counts.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



ReplaceNAs <- function(my.data) {
    na_row <- apply(my.data, 1, function(x) (sum(is.na(x))/ncol(my.data))*100)
    
    cat(paste0("\n- The input data has between " , round(min(na_row), digits = 2), "% - ", round(max(na_row), digits = 2),"%", " missing values per row.\n- Variables (rows) with more than 70% missing values will be removed.\n"))
    
    
    removeNA <- which(as.vector(na_row) > 70)
    if(length(removeNA) > 0) {
        my.data <- my.data[-removeNA,]
        
    }
    na_col <- apply(my.data, 2, function(x) (sum(is.na(x))/nrow(my.data))*100)
    cat(paste0("\n- The input data has between " , round(min(na_col), digits = 2), "% - ", round(max(na_col), digits = 2),"%", " missing values per column.\n- Samples (rows) with more than 80% missing values will be removed. Variables (rows) with more than 50% missing values are imputed using the overall mean per sample.\n... Performing missing value imputation. N.B Uncertainty increases with number of missing values!.\n"))
    
    removeNA <- which(as.vector(na_col) > 80)
    if(length(removeNA) > 0) {
        my.data <- my.data[,-removeNA]
    }
    
    still.NA <- unique(as.vector(is.na(my.data)))
    
    if (TRUE %in% still.NA) {
        varnames <-  rownames(my.data)
        
        if (checkData(as.matrix(my.data))[1] == FALSE) {
            my.data <- as.data.frame(lapply(my.data, as.numeric))
        }
        
        file <- try(my.data.lls <- data.frame(completeObs(llsImpute(as.matrix(my.data), k = 10, correlation="spearman", allVariables=TRUE))), silent =TRUE)
        hasNegB <- unique(as.vector(my.data < 0))
        hasNegA <- unique(as.vector(my.data.lls < 0))
        
        if (class(file) == "try-error" || TRUE %in% hasNegA & hasNegB == FALSE) {
            my.data <- impute.knn(as.matrix(my.data), rowmax = 0.7)
            my.data <- data.frame(my.data$data)
        } else {
            my.data <- my.data.lls
            rm(file)
        }
    }
    rownames(my.data) <- varnames
    return(my.data)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replace zeros:
# Takes as arguments;
    # my.data = a dataframe of expression/abundance counts.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



ReplaceZero <- function(my.data, my.group) {
    
    smallestGr <- min(as.numeric(table(my.group)))
    
    greaterthanBG <- apply(my.data, 1, function(x) sum(x > 0))
    lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr)

    if (length(lessthanBG) > 0) {
        my.data <- my.data[-lessthanBG,]
    }
    
    min_per_row <- as.vector(apply(my.data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(my.data)){
        my.data[i, my.data[i,] == 0] <- min_per_row[i]
    }
    return(my.data)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fitting Data Distributions:
# Takes as arguments;
# my.data = a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


NormalizeData <- function(my.variant, my.data, my.group, my.transform, my.standardize, my.data.original = NULL) {
    if (my.variant == "seq") {
        if (!is.null(my.data.original)) {
            my.data <- my.data.original
        }
        my.data <- DGEList(counts=my.data)
        design <- model.matrix(~0+my.group)
        keep <- filterByExpr(my.data, design)
        my.data <- my.data[keep,,keep.lib.sizes=FALSE]
        my.data <- calcNormFactors(my.data, method = "TMM")
        my.data <- voom(my.data, design, plot=TRUE)
        cat("\n-v = seq. Data will be filtered for lowly expressed variables, normalized and voom transformed.\n")
        
    } else if (my.variant %in% c("array", "ms", "other")) {
        if (my.transform == "log2") {
            my.data <- log2(my.data)
        } else if (my.transform == "logit") {
            my.data <- logit(my.data)
        } else if (my.transform == "log10"){
            my.data <- log10(my.data)
        } else {
            cat("\n-t is not specified for data, log transformation will NOT be performed.\n")
            my.data <- my.data
        }
        if (my.standardize == "mean") {
            my.data <- scale(my.data, scale = FALSE)
            cat(paste0("\n-v = array and -z = mean. Data will be mean centered.", NB, "\n"))
            
        } else if (my.standardize == "median") {
            rowmed <- apply(my.data,1,median)
            my.data <- my.data - rowmed
            cat(paste0("\n-v = array and -z = median. Data will be median centered.", NB, "\n"))
        } else if (!(my.standardize %in% c("mean", "median")) & my.variant == "array") {
            my.data <- normalizeBetweenArrays(my.data, method="quantile")
            cat(paste0("\n-v = array. Data will be normalized on the quantiles.",NB, "\n"))
        } else {
            cat("\n- No standardization requested. If argument -v is 'array', data will be normalized on quantile (NormalizeBetweenArrays), otherwise no normalization will be performed.\n")
        }
    } else {
        stop("\n- Option -v is mandatory and specifies data type (variant). Options are; array, seq, ms or other. See user manual for specifics.\n")
    }
    return(my.data)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fitting Data Distributions:
# Takes as arguments;
    # my.data = a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



FitDistributions <- function(my.data) {
    discretetype <- unique(as.vector(apply(my.data, 1, function(x) x%%1==0)))
    hasNeg <- unique(as.vector(my.data < 0))
    list.of.lists <- list()
    for(idx in 1:nrow(my.data)) {
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
                l <- list()
                if(exists("fit_n")) {
                    l[[length(l)+1]] <- fit_n
                }
                list.of.lists[[idx]] <-  l
            } else {
                try(fit_w  <- fitdist(as.numeric(my.data[idx,]), "weibull"), silent = TRUE)
                try(fit_g  <- fitdist(as.numeric(my.data[idx,]), "gamma"), silent = TRUE)
                try(fit_ln <- fitdist(as.numeric(my.data[idx,]), "lnorm"), silent = TRUE)
                try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
                l <- list()
                if(exists("fit_w")) {
                    l[[length(l)+1]] <- fit_w
                }
                if(exists("fit_g")) {
                    l[[length(l)+1]] <- fit_g
                }
                if(exists("fit_ln")) {
                    l[[length(l)+1]] <- fit_ln
                }
                if(exists("fit_n")) {
                    l[[length(l)+1]] <- fit_n
                }
                list.of.lists[[idx]] <-  l
            }
        }
        if (discretetype == TRUE) {
            try(fit_p <- fitdist(as.numeric(my.data[idx,]), "pois"), silent = TRUE)
            try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
            l <- list()
            if(exists("fit_p")) {
                l[[length(l)+1]] <- fit_p
            }
            if(exists("fit_n")) {
                l[[length(l)+1]] <- fit_n
            }
            list.of.lists[[idx]] <-  l
        }
    }
    names(list.of.lists) <- rownames(my.data)
    return(list.of.lists)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting Distributions:
# Takes as arguments;
    # my.data = a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!(See "Fitting Data Distributions" above)
    # list.of.lists = Output from the function "Fitting Data Distributions" (see above). The my.data and list.of.list should have the same dimentions, e.g. length of list == nrows of dataframe.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


PlotDistributions <- function(my.data, list.of.lists) {
    discretetype <- unique(as.vector(apply(my.data, 1, function(x) x%%1==0)))
    hasNeg <- unique(as.vector(my.data < 0))
    for(idx in 1:length(list.of.lists)) {
        pdf(paste0(names(list.of.lists)[idx], ".pdf"), height = 8, width = 12)
        par(mfrow=c(2,3))
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                descdist(as.numeric(my.data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("norm")
                plot.colors <- viridis(1)
            } else {
                descdist(as.numeric(my.data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("Weibull", "lognormal", "gamma", "norm")
                plot.colors <- viridis(4)
            }
        }
        if (!FALSE %in% discretetype) {
            descdist(as.vector(my.data[idx,]), discrete = TRUE,  boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
            plot.legend <- c("poisson", "norm")
            plot.colors <- viridis(2)
        }
        denscomp(list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors)
        cdfcomp (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, datapch=16)
        qqcomp  (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        ppcomp  (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        dev.off()
    }
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multidimensional Scaling Plot:
# Takes as arguments; 
	# my.data = a dataframe of expression/abundance counts 
	# my.group and my.lables = a vector of IDs for coloring and labeling (may be the same or different, length should be equal to ncol(dataframe))
	# my.cols = a vector of colors (one color for each group)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


MDSPlot <- function(my.data, my.group, my.labels, my.cols) {
    d<-dist(t(my.data))
    fit <- cmdscale(d,eig=TRUE, k=2)
    res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
    ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = my.group, colour = my.group), shape=21, size = 5, stroke = 0.1) + geom_text(aes(x=M1,y=M2, label= my.labels)) + guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=my.cols) + scale_colour_manual(values=rep("white", length(my.cols))) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function for Estimating Kmeans
# Takes as arguments;
# my.data = a dataframe of expression/abundance counts
# n = number of sample subsets to generate
# l = size of sample subsets
# k = number of groups to test
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


EstimateKmeans <- function(df, n) {
    BIC <- mclustBIC(df)
    mod1 <- Mclust(df, x = BIC)
    Ks <- as.numeric(summary(mod1, parameters = TRUE)$G)
    cat(paste0("\nCluster run complete - out of ", length(n), " in total..."))
    return(Ks)
}

PlotKmeans <- function(my.data, clus.list, k, my.labels, my.name) {
    res.list <- list()
    nclus <- unique(unlist(clus.list))
    if (unique(is.na(nclus)) == TRUE) {
        cat("\nNo 'best' ks could be determined. There may be little or poor clustering of samples. Ks 1:5 will be returned.\n")
        nclus <- k
    }
    nclus <- sort(nclus)
    for (idx in 1:length(nclus)) {
        set.seed(10)
        nth <- detectCores(logical = TRUE)
        Kclus <- Kmeans(t(my.data), nclus[[idx]], nthread=nth)
        Clusters <- as.factor(paste0("C",data.frame(Kclus$cluster)$Kclus.cluster))
        d<-dist(t(my.data))
        fit <- cmdscale(d,eig=TRUE, k=2)
        res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
        res$Clusters <- Clusters
        res.list[[idx]] <- Clusters
        Pal <- viridisLite::viridis(length(levels(as.factor(res$Clusters))))
        p <- ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = Clusters, colour = Clusters), shape=21, size = 3, stroke = 0.1) + guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=Pal) + scale_colour_manual(values=rep("white", length(Pal))) + theme_bw() + geom_text(label = my.labels, colour = "grey50", nudge_x = 0.25, nudge_y = 0.25, size =3) + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 14), axis.title=element_text(size=14)) + theme(legend.position = "top") + theme(axis.text=element_text(size=14)) + theme(axis.text = element_text(colour = "black"))
        ggsave(file=paste0("BestKmeans_C", as.character(nclus[[idx]]), ".pdf"), p, width = 10, height = 8)
    }
    names(res.list) <- paste0("Clus", nclus)
    res.df <- as.data.frame(res.list)
    return(res.df)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT FEATURES:
# Takes as arguments; 
	# my.contrast = a contrast between groups of interest
	# my.data = a dataframe of expression/abundance counts
        # my.design = a design matrix with all comparisons
	# my.coLFC and my.coFDR = cutoffs for logFC and FDR 
	# if blocking than a vector of patient IDs
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DAFeature <- function(my.contrast, my.data, my.design, coLFC, coFDR, my.block=NULL) {
    if(is.null(my.block)) {
        fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
    }
    else {
        corfit <- duplicateCorrelation(my.data, my.design, block=my.block)
        fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design, block = my.block, correlation=corfit$consensus), my.contrast))
    }
    tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
    
    up <- tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]
    down <- tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]
    
    up$name <- rownames(up)
    down$name <- rownames(down)
    
    up$dir <- rep("up", nrow(up))
    down$dir <- rep("down", nrow(down))
    
    final <- list(up, down)
    return(final)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE: 
# Takes as arguments; 
	# my.contrasts = all contrasts between groups of interest
	# my.data = a dataframe of expression/abundance counts
	# my.design = a design matrix with all comparisons
	# my.coLFC and my.coFDR = cutoffs for logFC and FDR
	# if blocking than a vector of patient IDs
	# TRUE/FALSE statment specifying output format, if TRUE the function return a vector of feature IDs only
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DAFeatureApply <- function(my.contrasts, my.data, my.design, coLFC, coFDR, my.block=NULL, my.vector) {
  my.features.l <- apply(my.contrasts, 2, function(x) DAFeature(x, my.data, my.design, coLFC, coFDR, my.block))
  if(my.vector == TRUE) {
    my.features <- do.call(rbind, lapply(my.features.l, function(x) do.call(rbind, x)))
    my.features <- unique(do.call(rbind, strsplit(rownames(my.features), "[.]"))[,2])
    return(my.features)
  }
  else {
    return(my.features.l)
  }
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LASSO FUNCTION
# Takes arguments:
# my.data = data marix of values
# my.group = vector of integers specifying group.
# IF my.multinorm=TRUE, then analysis will be multinomial, IF my.multinorm=FALSE, then then analysis will be binomial
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



LASSOFeature <- function(my.seed, my.data, my.group, my.LAorEN, my.validation=FALSE, my.multinorm=TRUE) {
    
    if (my.validation == TRUE) {
        
        ll <- list()
        llev <- levels(as.factor(my.group))
        
        for (idx in 1:length(llev)) {
            pos <- which(my.group == as.character(llev[idx]))
            ll[[idx]] <- pos
        }
        
        my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))
        
        
        my.data.val <- my.data[,my.samp]
        my.group.val <- as.integer(my.group[my.samp])
        
        my.data <- my.data[,-my.samp]
        my.group <- as.integer(my.group[-my.samp])
    }
    
    if(my.multinorm == TRUE) {
        set.seed(my.seed)
        nth <- detectCores(logical = TRUE)
        registerDoMC(cores=nth)
        my.fit <- cv.glmnet(x = t(my.data), y = my.group, family="multinomial", type.multinomial = "grouped", nfolds = 10, alpha = my.LAorEN, parallel=TRUE)
        my.coef <- coef(my.fit, s=my.fit$lambda.1se)
        my.ma <- as(my.coef[[1]], "matrix")
    } else {
        set.seed(my.seed)
        my.fit <- cv.glmnet(x = t(my.data), y = my.group, family = "binomial", type.measure = "class", nfolds = 10, alpha = my.LAorEN)
        my.coef <- coef(my.fit, s=my.fit$lambda.1se)
        my.ma <- as(my.coef, "matrix")
    }
    
    
    if (my.validation == TRUE) {
        meanerror <- round(as.numeric(mean(predict(my.fit, t(my.data.val), s=my.fit$lambda.1se, type="class") != my.group.val))*100, digits = 2)
        cat(paste0("\nOne LASSO/EN run out completed. Mean class error was = ", meanerror, "%\n"))
        
    } else {
        meanerror <- round(as.numeric(mean(predict(my.fit, t(my.data), s=my.fit$lambda.1se, type="class") != my.group))*100, digits = 2)
        cat(paste0("\nOne LASSO/EN run out completed. Mean cross-validation error was = ", meanerror, "%\n"))
    }
    
    
    my.ma <- names(my.ma[my.ma[,1] != 0, ])
    
    return(list(my.ma, meanerror))
    
    rm(my.fit)
    rm(my.coef)
}







# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING OUTPUT AS EXCEL FILE
# Takes as arguments:
# my.list= a list of dataframes from DA_feature_apply with p-values, FDRs, logFC ect.
# my.sheetname = name of sheet to save.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




TextOutput <- function(my.list, my.filename) {
    if (is.null(my.list)) {
        cat("\nDifferential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict?\n")
    } else {
        my.list <- do.call(rbind, unlist(my.list, recursive=FALSE))
        my.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(my.list))
        my.names <- gsub("arg.group|arg.sgroup", "", my.names)
        my.list$comparison <- my.names
        file <- try(write.table(my.list, paste0(my.filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            stop("\n- Differential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict? Also, check you metadata file has the correct minimum required columns for analysis.")
        }
    }
    return(my.list)
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments: 
	# my.truestatus = a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
	# my.colors = a vector with colors to use (a character vector with the length of the number of groups/levels).
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


GetColors <- function(my.truestatus, my.colors) {
    hm_col <- data.frame(status = levels(as.factor(my.truestatus)), mycolor = my.colors)
    true_status <- data.frame(status = my.truestatus)
    col <- join(true_status, hm_col)
    col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor))
    return(as.matrix(col$mycolor))
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR MAKING HEATMAP
# Takes as arguments:
        # my.DE = dataframe with counts for differential expressed/abundant features
        # my.gradient = color gradient to use for heatmap
        # my.colorshm = color pallet for groups full length
        # my.colors = color pallet for groups at levels
        # my.group = a vector of groups to color by
        # my.filename = name of output plot
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



MakeHeatmap <-  function(my.DE, my.gradient, my.colorshm, my.colors, my.group, my.filename, my.range) {
    pdf(paste0(my.filename,"_heatmap.pdf"), height = 13, width=11)
    heatmap.plus(as.matrix(scale(my.DE, scale = FALSE)), col=my.gradient, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(my.DE), labCol='', ColSideColors=cbind(my.colorshm, rep("white", length(my.group))), margins = c(14,8), cexCol=1, cexRow = 0.8)
    map <- makecmap(my.range[1]:my.range[2])
    map$colors <- viridis((length(map$breaks)-1), option="cividis")
    hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
    legend("topleft", legend = as.character(levels(as.factor(my.group))), fill=my.colors, cex = 0.7)
    dev.off()
}






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GENERATING SURVIVAL PLOT
# Takes as arguments;
	# my.survres = A survival object in the form of a list with a cox regression for each feature.
        # my.filename = name of output plot	
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
number_ticks <- function(n) {function(limits) pretty(limits, n)}


SurvivalCOX <- function(my.survres, my.filename, my.n) {
    survival_conf  <- lapply(my.survres, function(x) summary(x)[2, c(4,6,7)])
    survival_conf  <- data.frame(do.call(rbind, survival_conf))
    colnames(survival_conf) <- c("HR", "lower", "upper")
    survival_pval <-  as.numeric(unlist(lapply(my.survres, function(x) anova(x)[1,3])))
    survival_conf$pval <- survival_pval
    survival_conf$fdr <- p.adjust(survival_pval, method = "fdr", n=length(survival_pval))
    survival_conf$feature <- as.factor(names(my.survres))
    survival_conf$Significant <- as.factor(ifelse(survival_conf$fdr <=0.05, 1, 0))
    survival_conf$InverseFDR <- 1/survival_conf$fdr
    
    
    if(nrow(survival_conf) > 100) {
        n.splits <- floor(nrow(survival_conf)/my.n)
        n.splits <- chunk2(1:nrow(survival_conf),n.splits)
        features <- lapply(n.splits, function(x) survival_conf[x,])
        
        p1 <- lapply(features, function(x) ggplot(x, aes(feature, HR)) + geom_point(aes(colour = Significant, size = InverseFDR)) + scale_color_manual(values=c("grey70", "black")) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=number_ticks(10)) + xlab("Features") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x))))
        for (i in 1:length(p1)) {
            pdf(paste0(as.character(names(p1)[i]),"_individual_corrplots.pdf"), height=6, width=12)
            print(p1[i])
            dev.off()
        }

    } else {
        pdf(paste0(my.filename, "_corrplot.pdf"), height=6, width=12)
        p1 <- ggplot(survival_conf, aes(feature, HR)) + geom_point(aes(colour = Significant, size = InverseFDR)) + scale_color_manual(values=c("grey70", "black")) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=number_ticks(10)) + xlab("Features") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))
        print(p1)
        dev.off()
    }
    return(survival_conf)
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION ANALYSIS
# Takes as arguments: 
	# d1 and d2 = two dataframes with values to be correlated. These must have the same dimensions and rownames must the same.
	# my.feature = list of features to be correlated (e.g. a set of features), if all features are to be used then set feature to  rownames(df1) 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


CorrAnalysis <- function(d1, d2, my.filename) {
    
    my.names <- rownames(d1)
    
    d1 <- as.matrix(d1)
    d2 <- as.matrix(d2)
    
    pear_corr <- sapply(1:nrow(d1), function(i) cor(d1[i,], d2[i,], method = "pearson"))
    
    # pearson correlation p-values
    pear_p_val <- sapply(1:nrow(d1), function(i) cor.test(d1[i,], d2[i,], method = "pearson")$p.value)
    
    # correction for multiple testing with fdr
    fdr <- p.adjust(pear_p_val, method = "fdr")
    
    # make dataframe
    pear_corr_full <- data.frame(my.names, pear_corr, pear_p_val, fdr, log2(1/fdr))
    colnames(pear_corr_full) <- c("name", "cor_coef", "pval", "fdr", "Inverse_Scaled_FDR")
    
    corr.plot <- ggplot(pear_corr_full, aes(name, cor_coef)) +  geom_point(aes(colour = Inverse_Scaled_FDR), size=7, stroke = 0, shape = 16) + scale_color_viridis(begin = 0.2, end = 0.8, direction = -1, option="cividis") + scale_y_continuous(breaks=number_ticks(10)) + theme_bw() +  theme(panel.grid.major.x = element_blank()) + geom_text(aes(label=name), size=4, hjust = 0.8, vjust=-0.2, color="grey30") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16), legend.title = element_text(size=14)) + xlab("") + ylab("Correlation Coefficient") + geom_hline(yintercept=c(0.0, 0.5, -0.5), color=rep("grey30", 3))
    ggsave(paste0(my.filename, "_corrplot.pdf"), bg = "transparent", width=12, height=6, plot = corr.plot)
    return(pear_corr_full)
}



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION SCATTER PLOTS
# Takes as arguments;
	# my.data = a dataframe with expression/abundance counts for tissue or TIF
	# my.serumdata = a dataframe with expression/abundance counts for serum
	# my.filename = name of output plot
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

CorrelationPlots <- function(my.data, my.serumdata, my.features, my.filename) {
  my.data <- my.data[rownames(my.data) %in% my.features,]
  my.serumdata <- my.serumdata[rownames(my.serumdata) %in% my.features,]
  features <- lapply(1:nrow(my.data), function(x) data.frame(as.numeric(my.data[x,]), as.numeric(my.serumdata[x,])))
  features <- lapply(features, setNames, c("TIF_Tissue", "Serum"))
  p1 <- lapply(features, function(x) ggplot(x, aes(TIF_Tissue, Serum)) + geom_point(shape=1, size=2.5) + theme_bw() + geom_smooth(method=lm, colour="black") +  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16)))
  names(p1) <- my.features
  for (i in 1:length(p1)) {
      p1[[i]] <- p1[[i]] + ggtitle(names(p1)[i])
  }
  nc <- ceiling(length(features)/2)
  pdf(paste0(my.filename,"_individual_corrplots.pdf"), height=6, width=12)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR WGCNA - Most Interconnected Features.
# Takes as arguments;
# vec.of.modulecolors = vector of module colors
# my.moduleColors = module color object from WGCNA
# my.IC = interconnectivity object from WGCNA
# my.ExpData = dataset, features (genes, proteins ect) as columns and samples as rows.
# my.n = Number indicating top n% most interconnected features in dataset
# my.name = name of output plot(s)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ModuleIC <- function(vec.of.modulecolors, my.moduleColors, my.IC, my.ExpData, my.n, my.softPower, my.name) {
    my.modules <- list()
    for (idx in 1:length(vec.of.modulecolors)) {
        moduleC <- colnames(my.ExpData)[which(my.moduleColors == vec.of.modulecolors[[idx]])]
        if (length(moduleC) <= 100) {
            pdf(paste0(my.name, "_module", vec.of.modulecolors[[idx]], "_moduleHM.pdf"), height=14, width=14)
            plotNetworkHeatmap(my.ExpData, moduleC, weights = NULL, useTOM = TRUE, power = my.softPower, networkType = "unsigned", main = "Heatmap of the network")
            dev.off()
        }
        moduleCIC <- my.IC[rownames(my.IC) %in% moduleC, ]
        moduleCIC <- moduleCIC[order(moduleCIC$kWithin,decreasing = TRUE), ]
        modTOP <- ceiling((length(moduleC)/100) * my.n)
        modTOP <- moduleCIC[1:modTOP,]
        modTOP$Module <- paste0("module", rep(vec.of.modulecolors[[idx]],nrow(modTOP)))

        modTOP$gene <- as.factor(rownames(modTOP))
        bp <- ggplot(modTOP, aes(x=gene, y=kWithin))+ geom_bar(stat="identity", fill=as.character(vec.of.modulecolors[[idx]]), width = 0.4) + theme_minimal() + geom_text(aes(label=gene), color = "grey30", size = 2.5, angle = 90, hjust = -.05) + theme(axis.text.x=element_blank())
        ggsave(paste0(my.name,"_", vec.of.modulecolors[[idx]], "_moduleIC.pdf"), plot = bp, width = 14, height = 10)
        
        my.modules[[idx]] <- modTOP
        
    }
    names(my.modules) <- vec.of.modulecolors
    return(my.modules)
}




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR WGCNA
# Takes as arguments;
# my.data = dataframe with expression/abundance values
# my.thresholds = a vector of length two. First argument specifies the minimum size of a module and the second argument specifies the dissimilarity cutoff for merging of modules.
# my.name = name of output plot(s)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WGCNAAnalysis <- function(my.data, my.thresholds, my.name) {
    
    enableWGCNAThreads()
    powers <-  c(c(1:10), seq(from = 12, to=20, by=2));
    sft <- pickSoftThreshold(my.data, dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
    
    pdf(paste0(my.name,"_WGCNA_softpowerplot.pdf"))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n")
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],col="red")
    dev.off()
    
    # Set power to best softpower estimate
    softPower <- sft$powerEstimate
    
    if (is.na(softPower) | softPower > 9) {
        if (nrow(my.data) < 20) {
            softPower <- 9
        } else if (nrow(my.data) >= 20 & nrow(my.data) < 30) {
            softPower <- 8
        } else if (nrow(my.data) >= 30 & nrow(my.data) < 40) {
            softPower <- 7
        } else {
            softPower <- 6
        }
    }
    
    # Construct adjecancy and topology matrix
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    # Adjacency matrix and Topology Matrix
    adj <-adjacency(my.data, power = softPower)
    
    if(ncol(my.data) > 8000) {
        
        cat("Dataset too large for classic WGCNA, running blockwiseModules intead.\n")
        
        WGCNAres <- blockwiseModules(
        my.data,
        weights = NULL,
        checkMissingData = TRUE,
        blocks = NULL,
        maxBlockSize = 8000,
        blockSizePenaltyPower = 5,
        nPreclusteringCenters = as.integer(min(ncol(my.data)/20, 100*ncol(my.data)/5000)),
        randomSeed = 12345,
        corType = "pearson",
        maxPOutliers = 1,
        quickCor = 0,
        pearsonFallback = "individual",
        power = softPower,
        networkType = "unsigned",
        replaceMissingAdjacencies = FALSE,
        suppressTOMForZeroAdjacencies = FALSE,
        TOMType = "unsigned",
        TOMDenom = "min",
        getTOMs = NULL,
        saveTOMs = FALSE,
        saveTOMFileBase = "blockwiseTOM",
        deepSplit = 2,
        detectCutHeight = 0.995,
        minModuleSize = my.thresholds[1],
        minBranchEigennodeDissim = mergeCutHeight,
        tabilityCriterion = c("Individual fraction", "Common fraction"),
        reassignThreshold = 1e-6,
        minCoreKME = 0.5,
        minCoreKMESize = my.thresholds[1]/3,
        minKMEtoStay = 0.3,
        mergeCutHeight = 0.25)
        
        moduleColors <- WGCNAres$colors
        
    } else {
        TOM <- TOMsimilarity(adj)
        colnames(TOM) <- colnames(my.data)
        rownames(TOM) <- colnames(my.data)
        
        cat("Done with topology calculations.\n")
        
        dissTOM <- (1-TOM)
        
        # Dendogram
        geneTree <- hclust(as.dist(dissTOM),method="ward.D2")
        
        cat("Done with hclust.\n")
        
        # Construct Modules
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = my.thresholds[1])
        
        
        # Extract module color labels
        dynamicColors <- labels2colors(dynamicMods)
        
        # Change colors to viridis
        my.cols <- data.frame(dynamicColors)
        colnames(my.cols) <- c("oldcols")
        
        n.colors <- length(levels(as.factor(dynamicColors)))

        WGCNA.cols <- c("#5C6D70", "#E88873", "#B5BD89", "#E8AE68", "#104F55", "#4062BB", "#72195A", "#8ACB88", "#A997DF", "#4B296B", "#BC69AA", "#1C3144", "#C3D898", "#C2C1A5", "#B23A48", "#3A7D44", "#B6F9C9", "#FF8360", "#296EB4", "#E8C1C5", "#2E282A", "#AAA95A", "#34D1BF", "#3454D1", "#F2F4CB", "#898980", "#ADF5FF", "#372549", "#7DD181", "#F56476")
        
        if(n.colors > 30) {
            l <- n.colors-30
            WGCNA.cols <- c(WGCNA.cols, viridis(l, option="inferno"))
        }
        
        #colortransform <- data.frame(levels(as.factor(dynamicColors)), viridis(length(levels(as.factor(dynamicColors))), option="inferno"))
        colortransform <- data.frame(levels(as.factor(dynamicColors)), WGCNA.cols[1:n.colors])
        colnames(colortransform) <- c("oldcols", "colortransform")
        dynamicColors <- as.character(join(my.cols, colortransform)$colortransform)
        
        # Eigen features in each module
        MEList <- moduleEigengenes(my.data, colors = dynamicColors)
        MEs <- MEList$eigengenes
        MEDiss <- (1-cor(MEs))
        
        # Cluster module eigengenes
        if (ncol(MEDiss) <= 1) {
            print("N.B. You only have one module in your tree, please consider changing parameter -x. Specify a smaller minimum module size (default is 10) or decrease % dissimilarity for mergining modules (default is 25%).")
        } else {
            METree <- hclust(dist(MEDiss), method = "ward.D2")
        }
        
        cat("Done with hclust of MEDiss.\n")
        
        merge <- mergeCloseModules(my.data, dynamicColors, cutHeight = my.thresholds[2]/100, verbose = 3)
        mergedColors <- merge$colors
        mergedMEs <- merge$newMEs
        
        # Plot merged modules
        pdf(paste0(my.name,"_WGCNA_moduleTree.pdf"), height = 10, width = 12)
        plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleColors <- mergedColors
        colorOrder <- c("grey", standardColors(50))
        moduleLabels <- match(moduleColors, colorOrder)-1
    }
    
    
    # Interconnectivity of features in modules
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    IC <- intramodularConnectivity(adj,  moduleColors)
    mod.cols <- levels(as.factor(moduleColors))
    
    cat("Done with IC.\n")
    
    WGCNAres <- ModuleIC(mod.cols, moduleColors, IC, my.data, my.thresholds[3], softPower, my.name)
    return(WGCNAres)
}







# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which downloads the STRING database.
# Take arguments:
# my.geneIDs = String specifying type of gene ID to use to match IDs in STRING database. Options are: ensembl_peptide_id, hgnc_symbol, ensembl_gene_id, ensembl_transcript_id or uniprotswissprot.
# my.version = Version of STRING database to use. Default is version 11.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DownloadPPInt <- function(my.geneIDs, my.version = "11.0") {
    approvedGeneIDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
    setwd("..")
    file <- try(load("stringDB.Rdata"))
    setwd("Results/")
    
    if (class(file) == "try-error") {
        print("\nDownloading and preparing string database for protein-protein interactions - this may take a couple of minutes!\n")
        download.file(paste0("https://stringdb-static.org/download/protein.links.v",my.version,"/9606.protein.links.v",my.version,".txt.gz"), paste0("9606.protein.links.v", my.version, ".txt.gz"))
        stringDB <- data.table::fread(paste0("9606.protein.links.v", my.version, ".txt.gz"))
        
        stringDB$protein1 <- gsub("9606.", "", stringDB$protein1)
        stringDB$protein2 <- gsub("9606.", "", stringDB$protein2)
      
        if (my.geneIDs %in% approvedGeneIDs[-1]) {
            uqprot <- unique(c(stringDB$protein1, stringDB$protein2))
            uqprot <- split(uqprot, ceiling(seq_along(uqprot)/10000))
            print("Calling IDs from BiomaRt")
            ensemblMT <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
            #map <- unique(getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = uqprot, mart = ensemblMT))
            map <- lapply(uqprot, function(x) getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = x, mart = ensemblMT))
            map <- do.call("rbind", map)
            colnames(map) <- c("protein1", "ID1")
            stringDB <- merge(stringDB, map, by = "protein1")
            colnames(map) <- c("protein2", "ID2")
            stringDB <- merge(stringDB, map, by = "protein2")
            stringDB <- stringDB[,-c(1,2)]
        }
        
        stringDB <- stringDB[stringDB$combined_score > as.numeric(quantile(stringDB$combined_score)[2]),]
        setwd("..")
        save(stringDB, file = "stringDB.Rdata")
        setwd("Results/")
    }
    
    return(stringDB)
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which subsets a dataframe of differential expression analysis into a list of dataframes by comparison.
# Take arguments:
# my.data = a dataframe with results of differential expression analysis.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DFtoList <- function(my.data) {
    
    my.data$name <- gsub("_", "-", my.data$name)
    comparisons <- levels(as.factor(my.data$comparison))
    df.list <- list()
    
    for (idx in 1:length(comparisons)) {
        df <- my.data[my.data$comparison == as.character(comparisons[idx]),]
        df.list[[idx]] <- df[,-c(2:4,6)]
    }
    
    names(df.list) <- comparisons
    
    df.lengths <- as.numeric(unlist(lapply(df.list, function(x) nrow(x))))
    df.lengths.min <- which(df.lengths < 2)
    if (length(df.lengths.min) > 0) {
        df.list <- df.list[-df.lengths.min]
    }
    return(df.list)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which extracts the protein-protein interactions where each protein (gene) is differentially abundant (expressed).
# Take arguments:
# my.DB = output of the function "DownloadPPInt", e.g. a curated string database.
# my.Geneset = a dataframe with results of differential expression analysis.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


GetDeaPPInt <- function(my.DB, my.Genedat) {
    
    df.list <- DFtoList(my.Genedat)
    nodes.list <- list()
    
    for (idx in 1:length(df.list)) {
        # Merging StringDB and datasets
        df <- df.list[[idx]]
        df$ID1 <- df$name
        nodes <- merge(my.DB, df, by = "ID1")
        df$ID2 <- df$name
        nodes <- unique(merge(nodes, df, by = "ID2"))
        nodes <- nodes[order(nodes$combined_score, decreasing = TRUE),]
        df <- as.data.frame(nodes[, c(2,1,3,4,5,7,9,10,12,13)])
        colnames(df) <- c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "logFC.node2", "fdr.node2", "dir.node2", "comparison")
        df <- df[!duplicated(apply(df,1,function(x) paste(sort(x),collapse=''))),]
        nodes.list[[idx]] <- df
    }
    
    names(nodes.list) <- names(df.list)
    return(nodes.list)
}





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which extracts the miRNA-gene interactions where miRNAs (and potentially genes, if this data is available,) are differentially expressed.
# Take arguments:
# arg.miRset = string specifying what miRNA database to use, options are: mirtarbase (validated), targetscan (predicted) or tarscanbase (validated and predicted).
# my.miRdat = a dataframe with results of differential expression analysis (miRNAs).
# my.Geneset = a dataframe with results of differential expression analysis (genes)- by default this argument is set to NULL.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


GetGeneMiRNAInt <- function(arg.miRset, my.miRdat, my.Genedat = NULL) {
    
    miR.list <- DFtoList(my.miRdat)
    nodes.list <- list()
    
    for (idx in 1:length(miR.list)) {
        
        df <- miR.list[[idx]]
        df$node1 <- df$name
        
        
        if (arg.miRset == "mirtarbase") {
            mirs <- get_multimir(mirna = df$node1, table = "mirtarbase")@data
            if(length(mirs) > 0) {
                mirs <- unique(mirs[,3:4])
                mirs$score <- rep(1, nrow(mirs))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. Try re-running the pipeline with a lower cutoff for significance (-f) or change argument -i to either targetscan or tarscanbase")
            }
        } else if (arg.miRset == "targetscan") {
            mirs <- get_multimir(mirna = df$node1, table = "targetscan")@data
            if(length(mirs) > 0) {
                mirs <- mirs[,c(3,4,7)]
                mirs$score <- abs(as.numeric(mirs$score))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. Try re-running the pipeline with a lower cutoff for significance (-f) or change argument -i to either mirtarbase or tarscanbase")
            }
        } else if (arg.miRset == "tarscanbase") {
            mirsV <- get_multimir(mirna = df$node1, table = "mirtarbase")@data
            if(length(mirsV) > 0) {
                mirsV <- mirsV[, 3:4]
                mirsV$score <- rep(1, nrow(mirsV))
            }
            mirsP <- get_multimir(mirna = df$node1, table = "targetscan")@data
            if(length(mirsP) > 0) {
                mirsP <- mirsP[,c(3,4,7)]
            }
            if (length(mirsV) > 0 | length(mirsP) > 0) {
                mirs <- unique(rbind(mirsV, mirsP))
                mirs$score <- abs(as.numeric(mirs$score))
                
                mirs$IDs <- paste0(mirs[,1], "|", mirs[,2])
                mirs <- data.table(mirs[,-c(1,2)])
                mirs <-  data.frame(mirs[, lapply(.SD, mean), by=IDs])
                mirs <- data.frame(do.call(rbind, strsplit(mirs$IDs, "[|]")), as.numeric(mirs[,2]))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. No network can be generated!")
            }
        } else {
            stop("If the argument -i is specified, the second part of this argument it must be set to either mirtarbase, targetscan or tarscanbase!")
        }
        
        mirs <- mirs[!duplicated(apply(mirs,1,function(x) paste(sort(x),collapse=''))),]
        mirs <- mirs[mirs$score > as.numeric(quantile(mirs$score)[3]),]
        dfmiR <- merge(df, mirs, by = "node1")
        dfmiR <- dfmiR[, c(1,7,8,2,3,5,6)]
        colnames(dfmiR) <- c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "comparison")
        dfmiR$score <- dfmiR$score * 1000
        dfmiR <- dfmiR[-grep("^$", dfmiR$node2),]
        
        if(is.null(my.Genedat)) {
            dfmiR$logFC.node2 <- as.numeric(rep(0, nrow(dfmiR)))
            dfmiR$fdr.node2 <- rep(NA, nrow(dfmiR))
            dfmiR$dir.node2 <- rep(NA, nrow(dfmiR))
        }
        nodes.list[[idx]] <- dfmiR
    }
    
    names(nodes.list) <- names(miR.list)
    
    if(!is.null(my.Genedat)) {
        
        miRgene.list  <- list()
        
        gene.list <- DFtoList(my.Genedat)
        
        overlap <- intersect(names(gene.list), names(nodes.list))
        
        nodes.list <- nodes.list[names(nodes.list) %in% overlap]
        gene.list <- gene.list[names(gene.list) %in% overlap]
        
        for (idx in 1:length(gene.list)) {
            
            df <- gene.list[[idx]]
            df$node2 <- df$name
            
            dfmiRgene <- merge(df, nodes.list[[idx]], by = "node2")
            dfmiRgene <- dfmiRgene[, c(7,1,8,9,10,11,2,3,5,6)]
            colnames(dfmiRgene) <-c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "logFC.node2", "fdr.node2", "dir.node2", "comparison")
            dfmiRgene <- dfmiRgene[dfmiRgene$dir.node1 != dfmiRgene$dir.node2, ]
            
            miRgene.list[[idx]] <- dfmiRgene
        }
        nodes.list <- miRgene.list
        names(nodes.list) <- overlap
    }
    return(nodes.list)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function for combining miRNA-gene interactions and protein-protein (gene-gene) interactions.
# Take arguments:
# nodes.list1 = List of networks, i.e. output of the function "GetGeneMiRNAInt".
# nodes.list2 = List of networks, i.e. output of the function "GetDeaPPInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



MatchPPGmiRInt <- function(nodes.list1, nodes.list2) {
    
    nodes.list <- list()
    
    overlap <- intersect(names(nodes.list1), names(nodes.list2))
    
    nodes.list1 <- nodes.list1[names(nodes.list1) %in% overlap]
    nodes.list2 <- nodes.list2[names(nodes.list2) %in% overlap]
    
    for (idx in 1:length(nodes.list1)) {
        
        df1 <- nodes.list1[[idx]]
        df2 <- nodes.list2[[idx]]
        
        nodes <- rbind(df1[df1$dir.node1 == "up",], df1[df1$dir.node1 == "down",], df2[df2$dir.node1 == "up",], df2[df2$dir.node1 == "down",])
        nodes.list[[idx]] <- nodes
    }
    
    names(nodes.list) <- overlap
    return(nodes.list)
}








# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to write out interaction lists and trim interactions for plotting
# Take arguments:
# my.nodes.list = List of networks, i.e. output of the function "GetGeneMiRNAInt" or ourput of "GetDeaPPInt" or output of "MatchPPGmiRInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


TrimWriteInt <- function(my.nodes.list) {
    
    node.lengths <- as.numeric(unlist(lapply(my.nodes.list, function(x) nrow(x))))
    node.lengths.min <- which(node.lengths < 2)
    
    if (length(node.lengths.min) > 0) {
        my.nodes.list <- my.nodes.list[-node.lengths.min]
    }
    
    set.names <- names(my.nodes.list)
    trimlist <- list()
    
    for (idx in 1:length(my.nodes.list)) {
        p.nodes <- my.nodes.list[[idx]]
        p.nodes$myrank <- rowMeans(apply(cbind(abs(p.nodes$logFC.node1), abs(p.nodes$logFC.node2)), 2, rank))
        p.nodes <- p.nodes[order(p.nodes$myrank, decreasing = TRUE),]
        write.table(p.nodes, paste0(set.names[[idx]],"_AllInteractions.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
        
        if (nrow(p.nodes) > 100) {
            miRidx <- grep("miR|let", p.nodes$node1)
            if (length(miRidx) > 0) {
                miR.nodes <- p.nodes[miRidx,]
                
                if (sum(miR.nodes$logFC.node2) == 0) {
                    miRidx <- which(miR.nodes$score > as.numeric(quantile(miR.nodes$score)[4]))
                    miR.nodes <- miR.nodes[miRidx,]
                }
                
                p.nodes <- p.nodes[-miRidx,]
                
                if(nrow(miR.nodes) >= 100) {
                    p.nodes <- miR.nodes[1:100,]
                } else {
                    p.nodes <- p.nodes[1:(100-nrow(miR.nodes)),]
                    p.nodes <- rbind(p.nodes, miR.nodes)
                }
            } else {
                p.nodes <- p.nodes[1:100,]
            }
        }
        
        p.info <- data.table(c(as.character(p.nodes$node1), as.character(p.nodes$node2)), c(p.nodes$logFC.node1, p.nodes$logFC.node2))
        colnames(p.info) <- c("name", "logFC")
        p.info<- data.frame(p.info[, .(Freq = .N), by = .(name, logFC)])
        p.info <- p.info[order(p.info$Freq, decreasing = TRUE),]
        p.info$group <- 1:nrow(p.info)
        vircls <- viridis(2, direction = -1 , end = 0.9, option = "cividis")
        #vircls <- viridis(2, end = 0.6, direction = -1, option = "magma")
        p.info$calfb <- ifelse(p.info$logFC > 0, vircls[1], "grey50")
        p.info$calfb <- ifelse(p.info$logFC < 0, vircls[2], as.character(p.info$calfb))
        
        myorder <- unique(as.vector(t(data.frame(p.nodes[,1:2]))))
        p.info <- p.info[match(myorder, p.info$name),]
        
        trimmed <- list(p.nodes, p.info)
        names(trimmed) <- c("p.nodes", "p.info")
        trimlist[[idx]] <- trimmed
    }
    names(trimlist) <- set.names
    return(trimlist)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function plotting interaction networks.
# Take arguments:
# my.trimmed.list = List of trimmed networks, i.e. output of the function "TrimWriteInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



PlotInt <- function(my.trimmed.list) {
    
    trim.lengths <- as.numeric(unlist(lapply(my.trimmed.list, function(x) nrow(x))))
    trim.lengths.min <- which(trim.lengths < 10)
    
    if (length(trim.lengths.min) > 0) {
        my.trimmed.list <- my.trimmed.list[-trim.lengths.min]
    }
    
    trimmed.names <- names(my.trimmed.list)
    
    for (idx in 1:length(my.trimmed.list)) {
        
        p.nodes <- my.trimmed.list[[idx]]$p.nodes
        p.info <- my.trimmed.list[[idx]]$p.info
        
        final.nodes <- as.matrix(p.nodes[, 1:2])
        degrees <- as.numeric(p.info$Freq)
        names(degrees) <- p.info$name
        meta <- data.frame(p.info$group, degrees, p.info$name, ind=1:nrow(p.info))
        my.order <- as.integer(meta[order(meta$degrees, decreasing = TRUE),]$ind)
        
        tiff(paste0(trimmed.names[[idx]],"_TopInteracions.tiff"), width = 14, height = 10, units = 'in', res = 600)

        #cex.nodes = 2^(log2(abs(p.info$logFC))/10)
        
        arcplot(final.nodes, ordering=my.order, labels=p.info$name, cex.labels=0.6,
        show.nodes=TRUE, col.nodes=p.info$calfb, bg.nodes=p.info$calfb,
        cex.nodes = (log(degrees)/2)+0.2, pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = log2(p.nodes$score/100)*1.5)
        dev.off()
    }
}






