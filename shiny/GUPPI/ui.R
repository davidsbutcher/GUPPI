
library(magrittr)
library(GUPPI)
library(dplyr)
library(assertthat)
library(tictoc)
library(RPushbullet)
library(glue)
library(feather)
library(shiny)
library(shinydashboard)
library(purrr)
library(readxl)
library(Peptides)
library(stringr)
library(forcats)
library(rmarkdown)
library(flexdashboard)
library(tibble)
library(knitr)
library(UpSetR)
library(pander)
library(ggplot2)
library(DT)
library(furrr)
library(future)
library(plotly)
library(ggthemes)
library(viridis)
library(waffle)
library(RSQLite)
library(DBI)
library(GO.db)
library(fs)
library(dbplyr)
library(tidyr)
library(Biobase)
library(UniProt.ws)
library(AnnotationDbi)

options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Define UI for application that draws a histogram
shinyUI(
    fixedPage(

        titlePanel("GUPPI beta shiny app"),

        # SIDEBAR

        sidebarLayout(
            sidebarPanel(
                fileInput(
                    "tdrep",
                    "Upload .tdReport file",
                    accept = c(".tdReport"),
                    multiple = TRUE
                ),
                br(),
                numericInput(
                    "taxon",
                    "Taxon number",
                    value = "9606"
                ),
                br(),
                actionButton(
                    "start",
                    "Make it happen"
                ),
                br(),
                downloadButton("downloadReport", label = "Download Report"),
            ),

            # MAIN PANEL

            mainPanel(
                textOutput("confirm")
            )
        )
    )
)
