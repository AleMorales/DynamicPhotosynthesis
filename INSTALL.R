# This script will install the dependencies required to run the R script inside the Code folder
# An internet connection is required
install.packages(c("dplyr", "ggplot2", "broom", "readr" ,"Rcpp", "RcppArmadillo", 
                   "rjson", "GetoptLong", "R6", "plyr", "Hmisc", 
                   "doParallel", "matrixStats", "abind", "Deriv"))
curdir = getwd()
setwd("Packages")
system("R CMD INSTALL RcppSundials --preclean --clean")
system("R CMD INSTALL SimulationModels --preclean --clean")
system("R CMD INSTALL fad --preclean --clean")
system("R CMD INSTALL diff --preclean --clean")
system("R CMD INSTALL RcppNLopt --preclean --clean")
system("R CMD INSTALL nlfit --preclean --clean")
system("R CMD INSTALL MiniModel --preclean --clean")
setwd(curdir)