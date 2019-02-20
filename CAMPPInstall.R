#Author: Thilde Bagger Terkelsen
#Contact: thilde@cancer.dk
#Place of employment: Danish Cancer Society Research Center
#Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




													### Install Packages ###





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list.of.packages <- c("getopt","limma", "edgeR", "glmnet", "sva", "openxlsx", "xlsx", "ggplot2", "heatmap.plus", "plyr", "data.table", "viridis", "squash", "survcomp", "survminer", "scales", "rms", "stackoverflow", "WGCNA", "fitdistrplus", "impute", "pcaMethods", "pROC", "VennDiagram")


missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (length(missing.packages) > 0) {
    cat('\nThe following package(s) is(are) not installed: ', '"', missing.packages, '"', '\nAn R script named "missing_packages.R" has been created to help you install the packages needed.\nTo install these packages do the following: Open R on the server (done by simply writing R at the prompt). In R source the script written out by entering:\nsource("missing_packages.R")\nYou will be asked for a CRAN mirror, pick the one closest to your geographical location.\n\n')
    cat('\n#Adapted from the original work by Joshua Wiley, http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html\ninstall.packages.auto <- function(x) {\n\tx <- as.character(substitute(x))\n\tif(isTRUE(x %in% .packages(all.available=TRUE))) {\n\t\teval(parse(text = sprintf("require(\\"%s\\")", x)))\n\t} else {\n\t\teval(parse(text = sprintf("install.packages(\\"%s\\", dependencies = TRUE)", x)))\n\t}\n\tif(isTRUE(x %in% .packages(all.available=TRUE))) {\n\t\teval(parse(text = sprintf("require(\\"%s\\")", x)))\n\t} else {\n\t\tsource("http://bioconductor.org/biocLite.R")\n\t\teval(parse(text = sprintf("biocLite(\\"%s\\")", x)))\n\t\teval(parse(text = sprintf("require(\\"%s\\")", x)))\n\t}\n}', file="CAMPPmissingpackages.R",append=TRUE)
    cat(paste0('\ninstall.packages.auto(','"', missing.packages, '"', ')'),'\nquit()', file="CAMPPmissingpackages.R",append=TRUE)
} else {
    cat('All packages are installed, ready to run analysis.\n')
}



