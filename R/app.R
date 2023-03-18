#################
#Wrapping the App inside a function
#################
#' Start the A2TEA.WebApp
#' @param ... just leave empty
#' @export
A2TEA_App <- function(...) {

#start the App
#shinyApp(ui = ui,
#         server = server)

#credits to DeanAttali -
#https://stackoverflow.com/a/49623819
appDir <- system.file("webapp", package = "A2TEA.WebApp")
  if (appDir == "") {
    stop("Could not find webapp Try re-installing `A2TEA.WebApp`.", call. = FALSE)
  }
shiny::runApp(appDir, display.mode = "normal")

#######################
#ending App function call wrap
#######################
}
