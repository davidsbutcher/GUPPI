
library(shiny)

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
