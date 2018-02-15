#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




													### FUNCTIONS ###
													
													
													
													

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BOXPLOTS AND VIOLINPLOTS
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function for generating boxplots: 
# Takes as arguments; 
	# my.data = a dataframe of expression/abundance counts
	# my.name = a name for the plot (given as a sting)
	# my.group = a vector of groups
	# my.sn = a vector of numbers indicating number of samples in each group
	# my.cols = a vector of colors for each box (one color for each group)
	
my.boxplots <- function(my.data, my.name, my.group, my.sn, my.cols) {
  png(paste0(my.name,".png"), height = 800, width = 1200)
  p1 <- apply(my.data, 1, function(x) ggplot(,aes(my.group, as.numeric(x))) + geom_boxplot(aes(fill = my.group)) + scale_fill_manual(values=my.cols) + theme_bw() + theme(legend.position="none") + geom_text(data=data.frame(), aes(x=names(c(by(as.numeric(x), my.group, median))), y=c(by(as.numeric(x), my.group, median)), label=my.sn), col='black', size=6) + labs(x = "BC Subtypes", y = "feature Abundance"))
  nc <- ceiling(nrow(my.data)/4)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}

# Function for generating violin plots: 
# Takes as arguments; 
	# my.data = a dataframe of expression/abundance counts
	# my.name = a name for the plot (given as a sting) 
	# my.group = a vector of groups
	# my.cols = a vector of colors for each box (one color for each group)
	
my.violin <- function(my.data, my.name, my.group, my.cols) {
  pdf(paste0(my.name,".pdf"), height = 6, width = 10)
  p1 <- apply(my.data, 1, function(x) ggplot(,aes(my.group, as.numeric(x))) + geom_violin(aes(fill = my.group), trim = FALSE) + stat_summary(fun.y=median, geom="point", size=2, color="black") + scale_fill_manual(values=my.cols) + theme_bw() + theme(legend.position="none") + theme(axis.text = element_text(size=13, colour = "black"), axis.title = element_text(size=13)) + labs(x = "feature", y = "Log abundance"))
  nc <- ceiling(nrow(my.data)/2)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multidimensional Scaling Plot: 
# Takes as arguments; 
	# my.data = a dataframe of expression/abundance counts 
	# my.group and my.lables = a vector of IDs for coloring and labeling (may be the same or different, length should be equal to ncol(dataframe))
	# my.cols = a vector of colors (one color for each group)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res) 
  p + geom_point(aes(x=M1,y=M2,color=my.group)) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
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



DA_feature <- function(my.contrast, my.data, my.design, coLFC, coFDR, my.block=NULL) {
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


DA_feature_apply <- function(my.contrasts, my.data, my.design, coLFC, coFDR, my.block=NULL, my.vector) {
  my.features.l <- apply(my.contrasts, 2, function(x) DA_feature(x, my.data, my.design, coLFC, coFDR, my.block)) 
  if(my.vector == TRUE) {
    my.features <- do.call(rbind, lapply(my.features.l, function(x) do.call(rbind, x)))
    my.features <- unique(do.call(rbind, strsplit(rownames(my.features), "[.]"))[,2])
    return(my.features)
  }
  else {
    return(my.features.l)
  }
}


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS THE "DA_feature_apply" FROM ABOVE.
# Takes as arguments; 
	# my.data = a dataframe of expression/abundance counts
	# my.group = a vector of groups do perform contrasts on (same length as ncol(dataframe))
	# my.coLFC and my.coFDR = a cutoff for logFC and FDR
	# if remove is different from NULL, a vector of indices to remove must be supplied

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DA_all_contrasts <- function(my.data, my.group, my.logFC, my.FDR, my.levels, my.remove=NULL) {
  if (!is.null(my.remove)) {
    my.data <- my.data[, -my.remove]
    my.group <- my.group[-my.remove]
  }
    G <- factor(as.character(my.group), levels=my.levels)
    combinations<- data.frame(t(combn(paste0("G", levels(G)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    mod_design <-  model.matrix(~0+G)
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))
    my.DA <- DA_feature_apply(contrast.matrix, my.data, mod_design, my.logFC, my.FDR, NULL, FALSE)
    return(my.DA)
  }





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING OUTPUT AS EXCEL FILE
# Takes as arguments:
# my.list= a list of dataframes from DA_feature_apply with p-values, FDRs, logFC ect.
# my.sheetname = name of sheet to save.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




excel_output <- function(my.list, my.sheetname) {
    if (is.null(my.list)) {
        print("Differential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict?")
    } else {
        my.list <- do.call(rbind, unlist(my.list, recursive=FALSE))
        my.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(my.list))
        my.names <- gsub("arg.group", "", my.names)
        my.list$comparison <- my.names
        xlsx::write.xlsx(my.list, file=paste0(my.sheetname,".xlsx"), row.names = FALSE)
    }
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments: 
	# my.truestatus = a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
	# my.colors = a vector with colors to use (a character vector with the length of the number of groups/levels).
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_colors <- function(my.truestatus, my.colors) {
    if (length(my.colors) > length(levels(as.factor(my.truestatus)))) {
    my.colors <- my.colors[-1]
    }
    hm_col <- data.frame(levels(as.factor(my.truestatus)), my.colors)
    colnames(hm_col) <- c("status", "mycolor")
    true_status <- data.frame(my.truestatus)
    myorder <- 1:nrow(true_status)
    true_status$order <- myorder
    colnames(true_status) <- c("status", "order")
    col <- merge(true_status, hm_col, by="status", all.x =TRUE)
    col <- col[order(col$order),]
    col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor))
    return(as.matrix(col$mycolor))
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------$
# FUNCTION FOR MAKING HEATMAP
# Takes as arguments:
        # my.DE = dataframe with counts for differential expressed/abundant features
        # my.gradient = color gradient to use for heatmap
        # my.colors = color pallet for groups
        # my.group = a vector of groups to color by
        # my.filename = name of output plot
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------$



my_heatmap <-  function(my.DE, my.gradient, my.colors, my.group, my.filename) {
    pdf(paste0(my.filename,"_heatmap.pdf"))
    heatmap.plus(as.matrix(scale(my.DE, scale = FALSE)), col=my.gradient, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(my.DE), labCol='', ColSideColors=cbind(my.colors, rep("white", length(my.group))), margins = c(14,8), cexCol=1.2, cexRow = 1.3)
    map <- makecmap(-3:4)
    map$colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = length(map$breaks)-1)
    hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
    dev.off()
}






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GENERATING SURVIVAL PLOT
# Takes as arguments;
	# my.survres = A survival object in the form of a list with a cox regression for each feature.
        # my.filename = name of output plot	
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
number_ticks <- function(n) {function(limits) pretty(limits, n)}

my_survival <- function(my.survres, my.filename) {
  survival_conf  <- lapply(my.survres, function(x) summary(x)[2, c(4,6,7)])
  survival_conf  <- data.frame(do.call(rbind, survival_conf))
  colnames(survival_conf) <- c("HR", "lower", "upper")
  survival_pval <-  as.numeric(unlist(lapply(my.survres, function(x) anova(x)[1,3])))
  survival_conf$pval <- survival_pval
  survival_conf$fdr <- p.adjust(survival_pval, method = "fdr", n=length(survival_pval))
  survival_conf$feature <- as.factor(names(my.survres))
  survival_conf$sig <- as.factor(ifelse(survival_conf$fdr <=0.05, 1, 0))
  survival_conf$InverseFDR <- 1/survival_conf$fdr
  survival.plot <- ggplot(survival_conf, aes(x=feature, y=HR)) + geom_point(aes(colour = sig, size = InverseFDR)) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + scale_x_discrete(limits=as.character(names(survival.results))) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=number_ticks(10)) + xlab("Features") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))
  ggsave(paste0(my.filename, "_survivalplot.pdf"), plot = survival.plot)
  return(survival_conf)
}




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GENERATING SURVIVAL ggplot CURVES
# Takes as arguments;
	# my.data = a dataframe with expression/abundance counts
	# my.survivaldata = a dataframe with results of cox reagression
	# my. index = indices of features (features) to use
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my_survival_curve <- function(my.data, my.survivaldata, my.index) {
  features <- lapply(my.index, function(x) data.frame(as.numeric(my.data[x,]), my.survivaldata$time_to_Outcome_years, my.survivaldata$Outcome_status, my.survivaldata$Age_at_surgery))
  features <- lapply(features, setNames, c("feature", "years", "status", "age"))
  features <- lapply(features, function(x) coxph(Surv(years, status) ~ age + feature, data = x))
  return(features)
}






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION ANALYSIS
# Takes as arguments: 
	# d1 and d2 = two dataframes with values to be correlated. These must have the same dimensions and rownames must the same.
	# my.feature = list of features to be correlated (e.g. a set of features), if all features are to be used then set feature to  rownames(df1) 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


my_correlation <- function(d1, d2, my.features, my.filename) {
    
    d1 <- d1[rownames(d1) %in% my.features,]
    d2 <- d2[rownames(d2) %in% my.features,]
    
    my.names <- rownames(d1)
    
    pear_corr <- sapply(1:nrow(d1), function(i) cor(d1[i,], d2[i,], method = "pearson"))
    
    # pearson correlation p-values
    pear_p_val <- sapply(1:nrow(d1), function(i) cor.test(d1[i,], d2[i,], method = "pearson")$p.value)
    
    # correction for multiple testing with fdr
    fdr <- p.adjust(pear_p_val, method = "fdr")
    
    # make dataframe
    pear_corr_full <- data.frame(my.names, pear_corr, pear_p_val, fdr, log2(1/fdr))
    colnames(pear_corr_full) <- c("name", "cor_coef", "pval", "fdr", "Inverse_Scaled_FDR")
    
    corr.plot <- ggplot(pear_corr_full, aes(name, cor_coef)) +  geom_point(aes(colour = Inverse_Scaled_FDR), size=7) + scale_colour_gradient(low="lightskyblue1", high="navyblue") + scale_y_continuous(breaks=number_ticks(10)) + theme_bw() +  theme(panel.grid.major.x = element_blank()) + geom_text(aes(label=name), size=6, hjust = 0.8, vjust=-0.2, color="grey30") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16), legend.title = element_text(size=14)) + xlab("") + ylab("Correlation Coefficient") + geom_hline(yintercept=c(0.0, 0.5, -0.5), color=c("grey30","maroon3", "maroon3"))
    ggsave(paste0(my.filename, "_corrplot.pdf"), width=12, height=6, plot = corr.plot)
    return(pear_corr_full)
}



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION SCATTER PLOTS
# Takes as arguments;
	# my.data = a dataframe with expression/abundance counts for tissue or TIF
	# my.serumdata = a dataframe with expression/abundance counts for serum
	# my.filename = name of output plot
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my_correlation_plots <- function(my.data, my.serumdata, my.features, my.filename) {
  my.data <- my.data[rownames(my.data) %in% my.features,]
  my.serumdata <- my.serumdata[rownames(my.serumdata) %in% my.features,]
  features <- lapply(1:nrow(my.data), function(x) data.frame(as.numeric(my.data[x,]), as.numeric(my.serumdata[x,])))
  features <- lapply(features, setNames, c("TIF_Tissue", "Serum"))
  p1 <- lapply(features, function(x) ggplot(x, aes(TIF_Tissue, Serum)) + geom_point(shape=1, size=2.5) + theme_bw() + geom_smooth(method=lm) + ggtitle("") +  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16)))
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
# FUNCTION FOR MAKING MULTIPLE ggplotS IN ONE FIGURE - FUNCTION WAS OBAINED FROM
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
