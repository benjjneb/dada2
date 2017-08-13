################################################################################
#' Launch Shiny-DADA2.
#' 
#' @return Launches a shiny app
#' 
#' @importFrom shiny runApp
#' 
#' @export
#' 
shinyDADA2 = function(){
  # find and launch the app
  appDir <- system.file("shiny", package = "dada2")
  runApp(appDir, display.mode = "normal")
}
################################################################################
