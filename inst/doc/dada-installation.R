## ----bioc-install, eval=FALSE--------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(suppressUpdates = FALSE)
#  biocLite("ShortRead", suppressUpdates = FALSE)

## ----install-packages, eval=FALSE----------------------------------------
#  install.packages("path/to/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----install-github-example, eval=FALSE----------------------------------
#  install.packages("~/github/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----packageVersion------------------------------------------------------
packageVersion("dada2")

## ----bioc-install-missing, eval=FALSE------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("missing_package_1")
#  biocLite("missing_package_2")
#  # ... and so on

## ----install-packages-rev2, eval=FALSE-----------------------------------
#  install.packages("path/to/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----load-dada2, message=FALSE-------------------------------------------
library("dada2")

## ----documentation-example, eval=FALSE, message=FALSE--------------------
#  help(package="dada2")
#  ?derepFastq
#  ?dada

