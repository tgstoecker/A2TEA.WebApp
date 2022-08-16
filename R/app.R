#source('R/global.R')
#source('R/ui.R', local = TRUE)
#source('R/server.R', local = TRUE)

#################
#Wrapping the App inside a function
#################
#' Start the A2TEA.WebApp
#' @export
A2TEA_App <- function(...) {


#start the App
#shinyApp(ui = ui,
#         server = server)

#credits to DeanAttali -
#https://stackoverflow.com/a/49623819
appDir <- system.file("./inst/", package = "A2TEA.WebApp")
  if (appDir == "") {
    stop("Could not find myapp. Try re-installing `mypackage`.", call. = FALSE)
  }
shiny::runApp(appDir, display.mode = "normal")
  
#######################
#ending App function call wrap
#######################
}