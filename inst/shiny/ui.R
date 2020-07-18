
library(magrittr)
library(GUPPI)
library(viztools)
library(dplyr)
library(assertthat)
library(tictoc)
library(glue)
library(feather)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
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
library(writexl)
library(readxl)
library(sessioninfo)


options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Hidden tabs -------------------------------------------------------------

VT_parameter_tabs <-
    tabsetPanel(
        id = "params",
        type = "hidden",
        tabPanel(
            "blank",
            br(),
            br()
        ),
        tabPanel(
            "upset",
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                selectInput(
                    "upset_name",
                    "UpSet type",
                    choices = c("Protein", "Proteoform")
                )
            ),
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                textInput(
                    "upset_barcolor",
                    "Bar color",
                    "#4C4184"
                )
            )
        ),
        tabPanel(
            "intdeg",
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                selectInput(
                    "intdeg_name",
                    "Int. Deg. type",
                    choices = c("Protein", "Proteoform")
                )
            ),
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                textInput(
                    "intdeg_fillcolor",
                    "Fill color",
                    "#4C4184"
                )
            ),
            sliderInput(
                "intdeg_yrange",
                "Y range",
                0,
                100,
                100,
                step = 1
            )
        ),
        tabPanel(
            "heatmap",
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                selectInput(
                    "heatmap_name",
                    "Heatmap type",
                    choices = c("Protein", "Proteoform")
                )
            ),
            div(
                style="display: inline-block;vertical-align:top; width: 150px;",
                selectInput(
                    "heatmap_orientation",
                    "Orientation",
                    choices = c("h", "v")
                )
            ),
            sliderInput(
                "heatmap_binsize",
                "Bin size (Da)",
                500,
                5000,
                1000,
                step = 500
            ),
            hr(),
            h5("Leave ranges blank for automatic sizing"),
            numericRangeInput(
                "heatmap_axisrange",
                "Mass axis range (kDa)",
                NULL
            ),
            numericRangeInput(
                "heatmap_countrange",
                "Count range",
                NULL
            )
        ),
        tabPanel(
            "waffle",
            "No options for now"
        )
    )

# UI ----------------------------------------------------------------------

shinyUI(
    fixedPage(

        useShinyjs(),

        tags$head(
            tags$style(HTML("hr {border-top: 1px solid #A9A9A9;}")),
            tags$style(
                HTML(
                ".rank-list-container.custom-sortable {
                    background-color: #6699cc;
                    overflow-y: scroll;
                    height: 150px;
                    width: 150px;
                 }
                 .custom-sortable .rank-list-item {
                    font-size: 0.75em;
                }"
                )
            )
        ),

        titlePanel("GUPPI beta shiny app"),

        # SIDEBAR

        sidebarLayout(
            sidebarPanel(
                tabsetPanel(
                    id = "mainpanel",
                    type = "pills",
                    tabPanel(
                        "GUPPI",
                        hr(),
                        radioGroupButtons(
                            inputId = "tdrep_fileinput",
                            label = "File location",
                            choices = c(
                                "Upload",
                                "Local Filesystem"
                            )
                        ),
                        br(),
                        div(
                            id = "input_server",
                            fileInput(
                                "tdrep",
                                "Upload .tdReport file",
                                accept = c(".tdReport"),
                                multiple = TRUE
                            ),
                            "NOTE: Large reports can take a long time to upload!"
                        ),
                        div(
                            id = "input_local",
                            shinyFilesButton(
                                "tdrep_local",
                                "Select .tdReport file(s)",
                                "Select one or more tdReport files",
                                multiple = TRUE
                            ),
                            br()
                        ),
                        br(), br(),
                        radioGroupButtons(
                            inputId = "tdrep_fracs",
                            label = "Fraction assignment",
                            choices = c(
                                "Automatic",
                                "Manual"
                            )
                        ),
                        br(),
                        selectInput(
                            "taxon",
                            "Taxon number",
                            choices = c("83333","2097","6239")
                        ),
                        br(),
                        actionButton(
                            "GUPPIstart",
                            "Analyze tdReport"
                        ),
                        br(), br(),
                        downloadButton("downloadReport", label = "Download Report")
                    ),
                    tabPanel(
                        "viztools",
                        hr(),
                        br(),
                        tabsetPanel(
                            id = "viztoolspanel",
                            type = "tabs",
                            tabPanel(
                                "Make Plot",
                                br(),
                                selectInput(
                                    "file1",
                                    "Choose a tdReport",
                                    choices = "Analyze w/ GUPPI first"
                                ),
                                div(
                                    id = "input_VT",
                                    selectInput(
                                        inputId = "plot_type",
                                        "Plot type",
                                        choices = c(
                                            "UpSet" = "upset",
                                            "Int. Degree" = "intdeg",
                                            "Heatmap" = "heatmap",
                                            "Waffle" = "waffle"
                                        )
                                    ),
                                    hr(),
                                    VT_parameter_tabs,
                                    hr(),
                                    actionBttn(
                                        "VTstart",
                                        "Update Preview"
                                    ),
                                    br(),
                                    br(),
                                    downloadButton("downloadPDF", label = "Download PDF"),
                                    downloadButton("downloadSVG", label = "Download SVG"),
                                    downloadButton("downloadPNG", label = "Download PNG")
                                )
                            ),
                            tabPanel(
                                "Image Settings",
                                br(),
                                div(
                                    style="display: inline-block;vertical-align:top; width: 150px;",
                                    selectInput(
                                        "download_font",
                                        "Plot font",
                                        choices = c(
                                            "sans",
                                            "serif",
                                            "mono"
                                        )
                                    )
                                ),
                                div(
                                    style="display: inline-block;vertical-align:top; width: 150px;",
                                    radioGroupButtons(
                                        inputId = "download_unit",
                                        label = "Size unit",
                                        choices =
                                            c("cm", "inch"),
                                        justified = TRUE
                                    )
                                ),
                                numericInput(
                                    "download_width",
                                    "Image width",
                                    value = 20,
                                    step = 0.5
                                ),
                                numericInput(
                                    "download_height",
                                    "Image height",
                                    value = 12.5,
                                    step = 0.5
                                )
                                ,
                                sliderInput(
                                    "download_dpi",
                                    "Image DPI (PNG only)",
                                    min = 50,
                                    max = 600,
                                    value = 300,
                                    step = 50
                                )
                            )
                        ),
                        width = 4
                    ),
                    tabPanel(
                        "About",
                        br(),
                        includeMarkdown("about.md")
                    )
                )
            ),


            # MAIN PANEL

            mainPanel(
                htmlOutput("ULconfirm"),
                uiOutput(
                    "tdrep_fracassign"
                ),
                br(),
                textOutput("confirm"),
                br(),
                textOutput("error"),
                br(),
                plotOutput("outputPlot")

            )
        )
    )
)
