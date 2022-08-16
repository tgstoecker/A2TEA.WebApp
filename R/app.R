source('R/global.R')
source('R/ui.R', local = TRUE)
source('R/server.R')

#################
#Wrapping the App inside a function
#################
#' Start the A2TEA.WebApp
#' @export
A2TEA_App <- function(...) {


#start the App
shinyApp(ui = ui,
         server = server)

#######################
#ending App function call wrap
#######################
}