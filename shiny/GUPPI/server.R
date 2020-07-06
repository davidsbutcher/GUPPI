
library(assertthat)
library(tictoc)
library(RPushbullet)
library(glue)
library(feather)
library(shiny)
library(purrr)
library(readxl)
library(Peptides)
library(stringr)
library(rmarkdown)
library(dplyr)
library(magrittr)
library(RSQLite)
library(DBI)

# Define server logic required to draw a histogram
shinyServer(
    function(input, output) {

        observeEvent(
            input$start,
            {
                outputDir <- tempdir()
                tempReport <- tempfile(fileext = ".html")

                GUPPI::guppi(
                    dirname(input$tdrep$datapath),
                    basename(input$tdrep$datapath),
                    input$taxon,
                    outputDir,
                    fdr = 0.01,
                    makeDashboard = T,
                    dashboardPath = tempReport,
                    saveOutput = F,
                    usePB = F
                )

                output$downloadReport <-
                    downloadHandler(
                        filename = glue::glue(
                            "{format(Sys.time(), '%Y%m%d_%H%M%S')}_GUPPIreport.html"
                        ),
                        content = function(file) {
                            file.copy(tempReport, file)
                        }
                    )

                output$confirm <-
                    renderText(
                        {
                            "tdReport analyzed succesfully"
                        }
                    )

            }
        )

    }
)
