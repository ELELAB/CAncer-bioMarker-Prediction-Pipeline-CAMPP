#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




													### Install Packages ###





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list.of.packages <- c("BiocManager", "getopt","limma", "edgeR", "glmnet", "sva", "openxlsx", "ggplot2", "gridExtra", "heatmap.plus", "plyr", "data.table", "viridis", "squash", "survcomp", "survminer", "scales", "rms", "stackoverflow", "WGCNA", "fitdistrplus", "impute", "pcaMethods", "pROC", "VennDiagram", "mclust", "multiMiR", "biomaRt", "devtools", "arcdiagram")


missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (length(missing.packages) > 0) {
    cat('\nThe following package(s) is(are) not installed: ', '"', missing.packages, '"', '\nAn R script named "CAMPPmissingpackages.R" has been created to help you install the packages needed.\nTo install these packages do the following: Open R on the server (done by simply writing R at the prompt). In R source the script written out by entering:\n\nsource("CAMPPmissingpackages.R")\n\nYou will be asked for a CRAN mirror, pick the one closest to your geographical location.\n\n')
    cat('\n#Adapted from the original work by Joshua Wiley, http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html\ninstall.packages.auto <- function(x) {\n\tx <- as.character(substitute(x))\n\tif(isTRUE(x %in% .packages(all.available=TRUE))) {\n\t\teval(parse(text = sprintf("require(\\"%s\\")", x)))\n\t} else {\n\t\teval(parse(text = sprintf("install.packages(\\"%s\\", dependencies = TRUE)", x)))\n\t}\n\tif(isTRUE(x %in% .packages(all.available=TRUE))) {\n\t\teval(parse(text = sprintf("require(\\"%s\\")", x)))\n\t} else {\n\t\teval(parse(text = sprintf("BiocManager::install(\\"%s\\")", x)))\n\t}\n}', file="CAMPPmissingpackages.R",append=TRUE)
    
    
    if ("arcdiagram" %in% missing.packages) {
        missing.packages <- missing.packages[-c(grep("arcdiagram", missing.packages))]
        arc <- TRUE
    } else {
        arc <- FALSE
    }
    
    cat(paste0('\ninstall.packages.auto(','"', missing.packages, '"', ')'), file="CAMPPmissingpackages.R",append=TRUE)
    
    if (arc == TRUE) {
        cat(paste0('\nlibrary(devtools)\ninstall_github("https://github.com/gastonstat/arcdiagram")'),'\nquit()',file="CAMPPmissingpackages.R",append=TRUE)
    }
    
} else {
    cat('All packages are installed, ready to run analysis.\n')
}
