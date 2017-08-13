################################################################################
# Check that the currently-installed version of R
# is at least the minimum required version.
################################################################################
R_min_version = "3.4.0"
R_version = paste0(R.Version()$major, ".", R.Version()$minor)
if(compareVersion(R_version, R_min_version) < 0){
  stop("You do not have the latest required version of R installed.\n", 
       "Launch should fail.\n",
       "See: http://cran.r-project.org/ and update your version of R.")
} else {
  message("R version looks okay:\n", R.version$major, ".", R.version$minor)
}
################################################################################
# Make sure devtools is installed...
################################################################################
if(!require("devtools")){
  install.packages("devtools")
}
if(
  compareVersion(a = as.character(packageVersion("devtools")),
                 b = "1.13.2") < 0
){
  message("\nUpdating devtools package to latest release version...\n")
  install.packages("devtools")
}
################################################################################
# Install basic required packages if not available/installed.
################################################################################
install_missing_packages = function(pkg, version = NULL, verbose = TRUE){
  availpacks = .packages(all.available = TRUE)
  source("http://bioconductor.org/biocLite.R")
  missingPackage = FALSE
  if(!any(pkg %in% availpacks)){
    if(verbose){
      message("The following package is missing.\n",
              pkg, "\n",
              "Installation will be attempted...")
    }
    missingPackage <- TRUE
  }
  if(!is.null(version) & !missingPackage){
    # version provided and package not missing, so compare.
    if( compareVersion(a = as.character(packageVersion(pkg)),
                       b = version) < 0 ){
      if(verbose){
        message("Current version of package\n", 
                pkg, "\t", 
                packageVersion(pkg), "\n",
                "is less than required.
                Update will be attempted.")
      }
      missingPackage <- TRUE
    }
  }
  if(missingPackage){
    biocLite(pkg, suppressUpdates = TRUE)
  }
}
################################################################################
# Special snapshot version of shinyFiles = "0.6.3"
################################################################################
if(
  compareVersion(a = as.character(packageVersion("shinyFiles")),
                 b = "0.6.3") < 0
){
  message("\nAttempting to update/install shinyFiles package from joey711/shinyFiles on GitHub...\n")
  devtools::install_github("joey711/shinyFiles")
}

################################################################################
# Define list of package names and required versions.
################################################################################
deppkgs = c(
  data.table = "1.10.4",
  dada2 = "1.5.2",
  DT = "0.2",
  ggplot2 = "2.2.1",
  jsonlite = "1.5",
  magrittr = "1.5",
  markdown = "0.8",
  rmarkdown = "1.6",
  RColorBrewer = "1.1.2",
  scales = "0.4.1",
  ShortRead = "1.34.0",
  shiny = "1.0.3",
  shinyFiles = "0.6.3")
# Loop on package check, install, update
pkg1 = mapply(install_missing_packages,
              pkg = names(deppkgs), 
              version = deppkgs,
              MoreArgs = list(verbose = TRUE), 
              SIMPLIFY = FALSE,
              USE.NAMES = TRUE)
################################################################################
# Load packages that must be fully-loaded 
################################################################################
packagesToLoad = c(
  "data.table",
  "dada2",
  "DT",
  "magrittr",
  "ggplot2",
  "shiny",
  "shinyFiles"
)
for(i in packagesToLoad){
  library(i, character.only = TRUE)
  message(i, " package version:\n", packageVersion(i))
}
################################################################################
