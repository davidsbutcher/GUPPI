
library(magrittr)
library(GUPPI)
library(viztools)
library(dplyr)
library(assertthat)
library(tictoc)
library(glue)
library(tippy)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(colourpicker)
library(purrr)
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
library(xfun)
library(htmltools)
library(mime)


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
            colourInput(
               "upset_barcolor",
               "Bar color",
               value = "#4C4184",
               showColour = c("both"),
               palette = c("square"),
               allowTransparent = FALSE
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
            colourInput(
               "intdeg_fillcolor",
               "Fill color",
               value = "#4C4184",
               showColour = c("both"),
               palette = c("square"),
               allowTransparent = FALSE
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
               choices = c(
                  "Vertical" = "v",
                  "Horizontal" = "h"
               )
            )
         ),
         div(
            style="display: inline-block;vertical-align:top; width: 150px;",
            selectInput(
               "heatmap_fillScale",
               "Color palette",
               choices = c(
                  "Plasma" = "C",
                  "Magma" = "A", 
                  "Inferno" = "B",
                  "Viridis" = "D",
                  "Cividis" = "E"
               )
            )
         ),
         div(
            style="display: inline-block;vertical-align:top; width: 150px;",
            sliderInput(
               "heatmap_binsize",
               "Bin size (Da)",
               500,
               5000,
               1000,
               step = 500
            )
         ),
         hr(),
         h5("Leave ranges blank for automatic sizing"),
         div(
            style="display: inline-block;vertical-align:top; width: 150px;",
            numericRangeInput(
               "heatmap_axisrange",
               "Mass axis range (kDa)",
               NULL
            )
         ),
         div(
            style="display: inline-block;vertical-align:top; width: 150px;",
            numericRangeInput(
               "heatmap_countrange",
               "Count range",
               NULL
            )
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
      # shinythemes::themeSelector(),
      titlePanel("GUPPI"),
      theme = "maglab_theme_old.css",
      
      useShinyjs(),
      
      setBackgroundColor(
         color = c("#FFFFFF"),
         gradient = "linear",
         direction = "bottom"
      ),
      
      tags$head(
         tags$style(HTML("hr {border-top: 1px solid #A9A9A9;}")),
         tags$style(
            HTML(
               ".shiny-notification {
                        position:fixed;
                        top: calc(20%);
                        left: calc(40%);
                    }"
            )
         )
      ),
      
      sidebarLayout(
         
         # Sidebar -----------------------------------------------------------------
         
         sidebarPanel(
            tabsetPanel(
               id = "mainpanel",
               type = "pills",
               tabPanel(
                  "Analyze",
                  hr(),
                  tabsetPanel(
                     id = "guppipanel",
                     type = "tabs",
                     tabPanel(
                        "Upload File",
                        br(),
                        div(
                           id = "input_server",
                           fileInput(
                              "tdrep",
                              tippy(
                                 "Upload .tdReport file",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "dark",
                                 tooltip =
                                    "<span style='font-size:14px;'>Choose a local file to analyze. The file will be transferred to a temp directory on the server. WARNING: Large .tdReport files can take a long time to upload!</span>"
                              ),
                              accept = c(".tdReport"),
                              multiple = TRUE
                           )
                        ),
                        div(
                           id = "input_local",
                           shinyFilesButton(
                              "tdrep_local",
                              tippy(
                                 "Upload .tdReport file",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "dark",
                                 tooltip =
                                    "<span style='font-size:14px;'>Select a file on the local filesystem to analyze. The file will NOT be uploaded to the server. NOTE: This should only be used if running the Shiny app directly from a local R installation, not when using the web application!</span>"
                              ),
                              "Select one or more .tdReport files",
                              multiple = TRUE
                           ),
                           br(),
                           br(),
                        ),
                        selectInput(
                           "taxon",
                           tippy(
                              "Taxon number",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "dark",
                              tooltip =
                                 "<span style='font-size:14px;'>Choose the NCBI taxon number of the database that should be used to add data from the UniProt/Gene Ontology knowledgebases. This number should match the database used in the TDPortal/ProsightPD analysis.</span>"
                           ),
                           choices =
                              c(
                                 "83333 (E. coli K12)" = "83333",
                                 "4932 (S. cerevisiae)" = "4932",
                                 "3055 (C. reinhardtii)" = "3055",
                                 "2097 (M. genitalium)" = "2097",
                                 "6239 (C. elegans)" = "6239",
                                 "9606 (Homo sapiens)" = "9606",
                                 "10090 (Mus musculus)" = "10090"
                              )
                        ),
                        br(),
                        actionButton(
                           "GUPPIstart",
                           "Analyze tdReport(s)"
                        ),
                        br(), br(),
                        downloadButton("downloadReport", label = "GUPPI Report"),
                        br(),
                        downloadButton("downloadProteinReport", label = "Protein Report"),
                        br(),
                        downloadButton("downloadProteoformReport", label = "Proteoform Report")
                     ),
                     tabPanel(
                        "Advanced",
                        br(),
                        radioGroupButtons(
                           inputId = "tdrep_fileinput",
                           tippy(
                              "File selection",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "dark",
                              tooltip =
                                 "<span style='font-size:14px;'>Chooses whether files are uploaded to the server or read from the local filesystem. Only select local filesystem if running the Shiny app locally from the R console.</span>"
                           ),
                           choices = c(
                              "Upload to server",
                              "Local filesystem"
                           )
                        ),
                        br(),
                        radioGroupButtons(
                           inputId = "tdrep_fracs",
                           label =  
                              tippy(
                                 "Fraction assignment",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "dark",
                                 tooltip =
                                    "<span style='font-size:14px;'>Chooses whether fractions will be assigned to .raw/.mzML input files in the .tdReport based on filename or manually by the user. Automatic assignment will only work if filenames include “f”, “frac”, ”fraction”, “peppi”, or “gf” followed by a one- or two-digit number.</span>"
                              ),
                           choices = c(
                              "Automatic",
                              "Manual"
                           )
                        ),
                        br(),
                        selectInput(
                           "fdr",
                           tippy(
                              "False detection rate cutoff",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "dark",
                              tooltip =
                                 "<span style='font-size:14px;'>Chooses false detection rate cutoff for protein/proteoform results to be considered valid and included in GUPPI reports.</span>"
                           ),
                           choices =
                              c(
                                 "0.1%" = 0.001,
                                 "0.5%" = 0.005,
                                 "1.0%" = 0.01,
                                 "5.0%" = 0.05,
                                 "10%" = 0.1
                              ),
                           selected = 0.01
                        ),
                        radioGroupButtons(
                           inputId = "tdrep_static",
                           label =  
                              tippy(
                                 "Static dashboard",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "dark",
                                 tooltip =
                                    "<span style='font-size:14px;'>Controls whether the GUPPI report generated will be dynamic (i.e. include HTML widgets from Plotly and DataTables) or static (contain only static images and text).</span>"
                              ),
                           choices = c(
                              "TRUE",
                              "FALSE"
                           ),
                           selected = "FALSE"
                        )
                     )
                  )
               ),
               tabPanel(
                  "Visualize",
                  hr(),
                  tabsetPanel(
                     id = "viztoolspanel",
                     type = "tabs",
                     tabPanel(
                        "Make Plot",
                        br(),
                        selectInput(
                           "file1",
                           "Choose a tdReport",
                           choices = "Analyze a tdReport first"
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
                  hr(),
                  includeMarkdown("about.md")
               )
            )
         ),
         
         
         # Main Panel --------------------------------------------------------------
         
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
            
         ),
         
         fluid = FALSE
      )
      
      
   )
)
