library(shiny)
library(data.table)
library(ggplot2) 
library(plotly) 
library(preprocessCore)

# Loading ui.R and server.R
source("ui.R")
source("server.R")

# Run the Shiny app
shinyApp(ui = ui, server = server)

