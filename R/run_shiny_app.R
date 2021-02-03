#' Run Shiny App
#'
#' @description
#' A simple function which runs the GUPPI Shiny application in a local browser.
#'
#' @export
#'

run_shiny_app <- 
   function(
   ) {
      shiny::runApp(
       appDir = system.file("shiny", package = "GUPPI"),
       launch.browser = TRUE
      )
   }