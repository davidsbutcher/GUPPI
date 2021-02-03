
library(magrittr)
library(GUPPI)
library(viztools)
library(dplyr)
library(assertthat)
library(tictoc)
library(glue)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
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

options(repos = BiocManager::repositories())

# is_local <- Sys.getenv('SHINY_PORT') == ""

options(shiny.maxRequestSize = 1000*1024^2)

# Functions ---------------------------------------------------------------

find_newest_file <-
   function(
      path
   ) {
      purrr::map_chr(
         path,
         ~fs::dir_info(.x) %>%
            dplyr::filter(modification_time == max(modification_time)) %>%
            dplyr::pull("path")
      )
   }

kickout <-
   function(
      list,
      allowed_ext = c("tdReport", "csv", "xlsx")
   ) {
      
      # This function removes any element from the list of input files
      # which does not have one of the allowed
      # extensions or which has "deprecated" in its name
      
      for (i in rev(seq_along(list))) {
         
         if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {
            
            list[[i]] <- NULL
            
         } else if (any(grep("deprecated", list[[i]], fixed = TRUE)) == TRUE) {
            
            list[[i]] <- NULL
            
         }
      }
      
      return(list)
   }


get_data_path <-
   function(
      filedir,
      filename,
      extension
   ) {
      
      filesindir <-
         fs::dir_ls(
            filedir, recurse = TRUE, type = "file",
            regexp = paste0("[.]", extension, "$")
         ) %>%
         purrr::as_vector()
      
      if (purrr::map(
         filename,
         ~stringr::str_detect(filesindir, .x)
      ) %>%
      purrr::map(any) %>%
      purrr::as_vector() %>%
      all() == FALSE) stop("One or more input files not found")
      
      if (filename %>%
          as.list() %>%
          purrr::map(~stringr::str_subset(filesindir, .x)) %>%
          purrr::map(~any(length(.x) > 1)) %>%
          unlist() %>%
          any() == TRUE) stop("One or more input files found in multiple locations")
      
      filelist <-
         filename %>%
         purrr::map_chr(~stringr::str_subset(filesindir, .x)) %>%
         as.list() %>%
         kickout()
      
      names(filelist) <- seq(1, length(filelist))
      
      return(filelist)
      
   }

read_tdreport_filenames <-
   function(tdreport) {
      
      # Establish database connection. Keep trying until it works!
      
      safe_dbConnect <- safely(dbConnect)
      
      safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
      
      if (is.null(safecon[["result"]]) == TRUE) message("Connection failed, trying again!")
      
      iteration_num <- 1
      
      while (is.null(safecon[["result"]]) == TRUE & iteration_num < 10) {
         
         iteration_num <- iteration_num + 1
         
         message(
            paste0("\nTrying to establish database connection, attempt ", iteration_num)
         )
         
         safecon <-
            safe_dbConnect(
               RSQLite::SQLite(), ":memory:",
               dbname = tdreport,
               synchronous = NULL
            )
         
      }
      
      if (is.null(safecon[["result"]]) == TRUE) {
         
         stop("read_tdreport_filenames could not connect to TDreport")
         
      } else {
         
         message(paste0("\nConnection to ", basename(tdreport), " succeeded"))
         con <- safecon[["result"]]
         
      }
      
      output <-
         dplyr::tbl(con, "DataFile") %>%
         dplyr::select(Name) %>%
         dplyr::collect() %>%
         dplyr::pull()
      
      # Close database connection and return output table
      
      dbDisconnect(con)
      
      message("\nread_tdreport_filenames finished\n")
      
      return(output)
      
   }

add_fraction_shiny <- function(filename) {
   
   # This function attempts to parse a filename to extract information
   # about the fraction that a raw file corresponds to. This is only useful
   # for GELFrEE/PEPPI/other fractionated data
   
   if (stringr::str_detect(
      filename,
      "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}"
   ) == TRUE) {
      return(
         stringr::str_extract(
            filename,
            "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}"
         )
      )
   } else {
      return(NA)
   }
   
}

has_extension_multi <- function(vector, extension) {
   
   res <- 
      purrr::map_chr(
         vector,
         ~assertthat::has_extension(.x, extension)
      )
   
   as.logical(res)
}

# Server ------------------------------------------------------------------

shinyServer(
   function(input, output, session) {
      
      
      # Initial setup -----------------------------------------------------------
      
      
      # Hide panels which are only shown conditionally
      
      hide("input_local")
      hide("input_VT")
      
      # Disable buttons which are selectively enabled later
      
      disable("GUPPIstart")
      disable("VTstart")
      disable("plot_type")
      disable("downloadReport")
      disable("downloadProteinReport")
      disable("downloadProteoformReport")
      disable("downloadPDF")
      disable("downloadSVG")
      disable("downloadPNG")
      
      # Establish params to use for shinyFiles input (local only)
      
      volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
      
      shinyFileChoose(
         input = input,
         "tdrep_local",
         session = session,
         roots = volumes
      )
      
      
      # Modals ------------------------------------------------------------------
      
      modal_extension <- function(failed = FALSE) {
         modalDialog(
            textInput(
               "dataset", "Choose data set",
               placeholder = 'Try "mtcars" or "abc"'
            ),
            span(
               '(Try the name of a valid data object like "mtcars", ',
               'then a name of a non-existent object like "abc")'
            ),
            if (failed)
               div(tags$b("Invalid name of data object", style = "color: red;")),
            
            footer = tagList(
               modalButton("Cancel"),
               actionButton("ok", "OK")
            )
         )
      }      
      
      # Reactives ---------------------------------------------------------------
      
      # Create reactive expression to determine whether file is on local filesystem
      # depends on input$tdrep_fileinput buttons
      
      is_local <-
         reactive(
            {
               
               switch(
                  input$tdrep_fileinput,
                  `Upload to server` = return(FALSE),
                  `Local filesystem` = return(TRUE)
               )
               
            }
         )
      
      is_valid <-
         reactive(
            {
               if (is_local() == FALSE && is.null(tdrep_name()) == FALSE) {
                  all(has_extension_multi(tdrep_name(), "tdReport"))
               } else if (is_local() == TRUE && is.integer(input$tdrep_local) == FALSE) {
                  all(has_extension_multi(tdrep_name(), "tdReport"))
               } else {
                  FALSE
               }
            }
         )
      
      input_filenames <-
         reactive(
            {
               validate(
                  need(
                     is_valid(), "Invalid file"
                  )
               )
               if (is_local() == FALSE) {
                  unlist(
                     map(
                        input$tdrep$datapath,
                        ~read_tdreport_filenames(.x)
                     )
                  )
               } else if (is_local() == TRUE){
                  unlist(
                     map(
                        parseFilePaths(volumes, input$tdrep_local)$datapath,
                        ~read_tdreport_filenames(.x)
                     )
                  )
               }
               
            }
         )
      
      input_fracassign_manual <-
         reactive(
            {
               ui_list <- list()
               
               for (i in seq_along(input_filenames())) {
                  
                  ui_list[[i]] <-
                     selectInput(
                        glue("assfrac{i}"),
                        label = input_filenames()[[i]],
                        choices = c(1:12),
                        width = "150px"
                     )
               }
               
               renderUI(
                  {
                     validate(
                        need(
                           is_valid(), ""
                        )
                     )
                     tagList(
                        h3("Choose fraction for each input file"),
                        ui_list
                     )
                  }
               )
               
            }
         )
      
      input_fracassign_auto <-
         reactive(
            {
               frac_list <- list()
               
               for (i in seq_along(input_filenames())) {
                  
                  frac_list[[i]] <-
                     add_fraction_shiny(input_filenames()[[i]])
                  
                  if (is.na(frac_list[[i]]) == TRUE) {
                     
                     frac_list[[i]] <- ""
                     
                  }
                  
               }
               
               names(frac_list) <- input_filenames()
               
               if (any(frac_list == "") == TRUE) {
                  
                  renderUI(
                     {
                        validate(
                           need(
                              is_valid(), ""
                           )
                        )
                        tagList(
                           h3("Automatically assigned fractions"),
                           tags$p(strong("At least one fraction was not automatically assigned"), style="color:red;font-size:14pt;"),
                           renderTable(
                              tibble::enframe(
                                 frac_list,
                                 name = "filename",
                                 value = "fraction"
                              )
                           )
                        )
                     }
                  )
                  
               } else {
                  
                  renderUI(
                     {
                        validate(
                           need(
                              is_valid(), ""
                           )
                        )
                        tagList(
                           h3("Automatically assigned fractions"),
                           renderTable(
                              tibble::enframe(
                                 frac_list,
                                 name = "filename",
                                 value = "fraction"
                              )
                           )
                        )
                     }
                  )
               }
               
            }
         )
      
      
      # Create reactive expressions for tdreport path, name, GO location type
      
      tdrep_path <-
         reactive(
            {
               if (is_local() == FALSE) {
                  input$tdrep$datapath
               } else if (is_local() == TRUE) {
                  parseFilePaths(volumes, input$tdrep_local)$datapath
               }
            }
         )
      
      tdrep_name <-
         reactive(
            {
               if (is_local() == FALSE) {
                  input$tdrep$name
               } else if (is_local() == TRUE) {
                  parseFilePaths(volumes, input$tdrep_local)$name
               }
            }
         )
      
      
      assignments <-
         reactive(
            {
               if (input$tdrep_fracs == "Manual") {
                  
                  assfrac_list <-
                     as.list(input_filenames())
                  
                  for (i in seq_along(input_filenames())) {
                     
                     names(assfrac_list)[[i]] <-
                        eval(parse(text = glue("input$assfrac{i}")))
                     
                  }
                  
                  return(assfrac_list)
                  
               } else {
                  
                  return(NULL)
                  
               }
            }
         )
      
      GOLocType <-
         reactive(
            {
               switch(
                  isolate(input$taxon),
                  "83333" = "bacteria",
                  "4932" = "eukaryota",
                  "6906" = "eukaryota",
                  "2097" = "bacteria",
                  "6239" = "eukaryota",
                  "9606" = "eukaryota",
                  "10090" = "eukaryota"
               )
            }
         )
      
      
      # Listeners ---------------------------------------------------------------
      
      listener_fracassign <- 
         reactive(
            {
               list(
                  input$tdrep_fracs,
                  input$tdrep,
                  input$tdrep_local
               )
            }
         )     
      
      listener_upload <- 
         reactive(
            {
               list(
                  input$tdrep,
                  input$tdrep_local
               )
            }
         )      
      
      # Observers ---------------------------------------------------------------
      
      observeEvent(
         input$tdrep_fileinput,
         {
            if (is_local() == FALSE) {
               hide("input_local")
               shinyjs::show("input_server")
            } else if (is_local() == TRUE) {
               hide("input_server")
               shinyjs::show("input_local")
            }
         }
      )
      
      observeEvent(
         listener_upload(),
         {
            if(is.null(tdrep_name()) == FALSE && !(all(has_extension_multi(tdrep_name(), "tdReport")))) {
               showModal(
                  modalDialog(
                     title = "Invalid file",
                     "You must upload a tdReport file"
                  )
               )   
            }
         }
      )
      
      # Code to run after a tdRep is uploaded to the SERVER or LOCAL FILESYSTEM
      
      observeEvent(
         listener_upload(),
         {
            
            # Show the names of files that were uploaded
            
            output$ULconfirm <-
               renderUI(
                  {
                     validate(
                        need(
                           is_valid(), "No tdReport uploaded"
                        )
                     )
                     tagList(
                        h3("Uploaded files"),
                        renderTable(
                           data.frame("filename" = tdrep_name())
                        )
                     )
                  }
               )
            
            if (is_valid()) {
               enable("GUPPIstart")
            } else {
               disable("GUPPIstart")
            }
            
            # Get rid of the output plot (for when a new tdRep is analyzed
            # with a plot still shown)
            
            output$outputPlot <-
               NULL
            
         }
      )
      
      # Code to run after a tdRep is chosen on the LOCAL FILESYSTEM
      
      # observeEvent(
      #    input$tdrep_local,
      #    {
      #       
      #       # Show the names of files that were uploaded
      #       
      #       output$ULconfirm <-
      #          renderUI(
      #             {
      #                tagList(
      #                   h3("Uploaded files"),
      #                   renderTable(
      #                      data.frame(
      #                         "filename" = tdrep_name()
      #                      )
      #                   )
      #                )
      #             }
      #          )
      #       
      #       if (all(has_extension_multi(tdrep_name(), "tdReport"))) {
      #          enable("GUPPIstart")
      #       } else {
      #          disable("GUPPIstart")
      #       }
      #       
      #       
      #       # Get rid of the output plot (for when a new tdRep is analyzed
      #       # with a plot still shown)
      #       
      #       output$outputPlot <-
      #          NULL
      #       
      #    }
      # )
      # 
      # observeEvent(
      #    {
      #       listener_fracassign()
      #    },
      #    {
      #       
      #       output$outputPlot <- NULL
      #       
      #       if (input$tdrep_fracs == "Manual") {
      #          
      #          if (is_valid() == TRUE) {
      #             
      #             output$tdrep_fracassign <- input_fracassign_manual()
      #             
      #          }
      #          
      #       } else if (input$tdrep_fracs == "Automatic") {
      #          
      #          if (is_valid() == TRUE) {
      #             
      #             output$tdrep_fracassign <- input_fracassign_auto()
      #             
      #          }
      #          
      #       }
      #    }
      # )
      
      observeEvent(
         {
            listener_fracassign()
         },
         {
            
            if (is_valid() == TRUE) {
               
               if (input$tdrep_fracs == "Manual") {
                  output$tdrep_fracassign <- input_fracassign_manual()
               }
               
               if (input$tdrep_fracs == "Automatic") {
                  output$tdrep_fracassign <- input_fracassign_auto()
               }
               
            } else {
               
               output$tdrep_fracassign <- NULL
               
            }
            
         }
      )
      
      # Push the START button for GUPPI
      
      observeEvent(
         input$GUPPIstart,
         {
            
            # Clear fraction assignment and output plots
            
            output$tdrep_fracassign <-
               NULL
            
            output$outputPlot <-
               NULL
            
            # Hide viztools panel during analysis
            
            hide("input_VT")
            
            # Disable buttons which are selectively enabled later to prevent 
            # crash during analysis
            
            disable("GUPPIstart")
            disable("VTstart")
            disable("plot_type")
            disable("downloadReport")
            disable("downloadProteinReport")
            disable("downloadProteoformReport")
            disable("downloadPDF")
            disable("downloadSVG")
            disable("downloadPNG")
            
            # Update VT/Visualize dropdown menu
            
            updateSelectInput(
               session = session,
               inputId = "file1",
               choices = tdrep_name()
            )
            
            # Initiate progress window
            
            withProgress(
               message = "Analyzing tdReport, please wait...",
               value = 0,
               style = getShinyOption("progress.style", default = "notification"),
               {
                  outputDir <- tempdir()
                  tempReport <- tempfile(fileext = ".html", tmpdir = outputDir)
                  
                  # Rename the file if it's not local, i.e. if it has to be
                  # copied to a temp dir. Shiny gives a random name to
                  # files copied to a temp dir
                  
                  if (is_local() == FALSE) {
                     purrr::map2_chr(
                        tdrep_path(),
                        fs::path(
                           dirname(tdrep_path()),
                           tdrep_name()
                        ),
                        ~file.rename(.x, .y)
                     )
                  }
                  
                  setProgress(value = 0.5)
                  
                  GUPPI::guppi(
                     dirname(tdrep_path())[[1]],
                     tdrep_name(),
                     as.integer(input$taxon),
                     GOLocType = GOLocType(),
                     fractionAssignments = assignments(),
                     outputdir = outputDir,
                     fdr = as.numeric(input$fdr),
                     makeDashboard = T,
                     dashboardPath = tempReport,
                     staticDashboard = as.logical(input$tdrep_static),
                     saveOutput = T
                  )
                  
                  setProgress(value = 0.75)
                  
                  # Download handlers for GUPPI/protein/proteoforms report
                  
                  output$downloadReport <-
                     downloadHandler(
                        filename = 
                           glue::glue(
                              "{format(Sys.time(), '%Y%m%d_%H%M%S')}_GUPPI_report.html"
                           ),
                        content = 
                           function(file) {
                              file.copy(tempReport, file)
                           }
                     )
                  
                  # Download handler for protein report
                  
                  output$downloadProteinReport <-
                     downloadHandler(
                        filename = glue::glue(
                           "{format(Sys.time(), '%Y%m%d_%H%M%S')}_GUPPI_protein_report.xlsx"
                        ),
                        content = function(file) {
                           file.copy(
                              find_newest_file(
                                 fs::path(
                                    outputDir,
                                    "protein_results"
                                 )
                              ),
                              file
                           )
                        }
                     )
                  
                  # Download handler for protein report
                  
                  output$downloadProteoformReport <-
                     downloadHandler(
                        filename = glue::glue(
                           "{format(Sys.time(), '%Y%m%d_%H%M%S')}_GUPPI_proteoform_report.xlsx"
                        ),
                        content = function(file) {
                           file.copy(
                              find_newest_file(
                                 fs::path(
                                    outputDir,
                                    "proteoform_results"
                                 )
                              ),
                              file
                           )
                        }
                     )
                  
                  output$confirm <-
                     renderText(
                        {
                           "tdReport analyzed succesfully"
                        }
                     )
                  
                  shinyjs::show("input_VT")
                  
                  output$tdrep_fracassign <- NULL
                  
                  enable("VTstart")
                  enable("plot_type")
                  enable("downloadReport")
                  enable("downloadProteinReport")
                  enable("downloadProteoformReport")
                  enable("downloadPDF")
                  enable("downloadSVG")
                  enable("downloadPNG")
                  
                  # Save session info (for dev use)
                  
                  systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
                  
                  sessioninfo::session_info() %>%
                     capture.output() %>%
                     writeLines(
                        glue::glue("{outputDir}/session_info/{systime}_shinySessionInfo.txt")
                     )
                  
                  # Finish off the progress bar
                  
                  setProgress(value = 1)
                  
               }
            )
         }
      )
      
      
      observeEvent(
         input$plot_type,
         {
            updateTabsetPanel(session, "params", selected = input$plot_type)
         }
      )
      
      observeEvent(
         input$VTstart,
         {
            
            if (is_local() == FALSE) req(input$tdrep)
            if (is_local() == TRUE) req(input$tdrep_local)
            
            output$confirm <-
               NULL
            
            output$tdrep_fracassign <-
               NULL
            
            output$outputPlot <-
               renderPlot(
                  {
                     switch(
                        isolate(input$plot_type),
                        upset = eval(UpSetPlotExpr),
                        intdeg = eval(IntDegPlotExpr),
                        heatmap = eval(HeatmapPlotExpr),
                        waffle = eval(WafflePlotExpr)
                     )
                  },
                  width = 800,
                  height=800*(input$download_height/input$download_width)
               )
            
         }
      )
      
      # Plot expressions --------------------------------------------------------
      
      UpSetPlotExpr <-
         expr(
            {
               load_UpSet_data <-
                  reactive(
                     {
                        readxl::read_xlsx(
                           path =
                              fs::path(
                                 tempdir(),
                                 "protein_results_allhits",
                                 paste0(
                                    fs::path_ext_remove(input$file1), "_allhits"
                                 ),
                                 ext = "xlsx"
                              )
                        ) %>%
                           dplyr::select(
                              if ((input$upset_name) == "Protein") "AccessionNumber" else "ProteoformRecordNum",
                              fraction
                           ) %>%
                           dplyr::arrange(fraction) %>%
                           tidyr::pivot_wider(
                              names_from = fraction,
                              values_from =
                                 if ((input$upset_name) == "Protein") "AccessionNumber" else "ProteoformRecordNum",
                              values_fn = list,
                              names_prefix = "Frac_"
                           ) %>%
                           purrr::flatten() %>%
                           list() %>%
                           rlang::set_names("1")
                     }
                  )
               
               make_UpSet_plot(
                  isolate(load_UpSet_data()),
                  plotType = isolate(input$upset_name),
                  barColor = isolate(input$upset_barcolor)
               )
            }
         )
      
      IntDegPlotExpr <-
         expr(
            {
               load_intdeg_data <-
                  reactive(
                     {
                        
                        readxl::read_xlsx(
                           path =
                              fs::path(
                                 tempdir(),
                                 "protein_results_allhits",
                                 paste0(
                                    fs::path_ext_remove(input$file1), "_allhits"
                                 ),
                                 ext = "xlsx"
                              )
                        ) %>%
                           dplyr::select(
                              if ((input$intdeg_name) == "Protein") "AccessionNumber" else "ProteoformRecordNum",
                              fraction
                           ) %>%
                           dplyr::arrange(fraction) %>%
                           tidyr::pivot_wider(
                              names_from = fraction,
                              values_from =
                                 if ((input$intdeg_name) == "Protein") "AccessionNumber" else "ProteoformRecordNum",
                              values_fn = list,
                              names_prefix = "Frac_"
                           ) %>%
                           purrr::flatten() %>%
                           purrr::map(unique)
                        
                     }
                  )
               
               make_intersection_degree_plot(
                  isolate(load_intdeg_data()),
                  Yrange = c(0, as.integer(isolate(input$intdeg_yrange))),
                  plotType = isolate(input$intdeg_name),
                  fillColor = isolate(input$intdeg_fillcolor),
                  fontFamily = isolate(input$download_font)
               )
            }
         )
      
      HeatmapPlotExpr <-
         expr(
            {
               load_heatmap_data <-
                  reactive(
                     {
                        
                        allhits_xlsx <-
                           readxl::read_xlsx(
                              path =
                                 fs::path(
                                    tempdir(),
                                    "protein_results_allhits",
                                    paste0(
                                       fs::path_ext_remove(input$file1), "_allhits"
                                    ),
                                    ext = "xlsx"
                                 )
                           )
                        
                        if (input$heatmap_name == "Protein") {
                           allhits_xlsx <-
                              allhits_xlsx %>%
                              dplyr::group_by(
                                 AccessionNumber,
                                 fraction
                              )
                        } else if (input$heatmap_name == "Proteoform"){
                           allhits_xlsx <-
                              allhits_xlsx %>%
                              dplyr::group_by(
                                 ProteoformRecordNum,
                                 fraction
                              )
                        }
                        
                        allhits_xlsx %>%
                           dplyr::filter(
                              `GlobalQvalue` == min(`GlobalQvalue`)
                           ) %>%
                           dplyr::filter(
                              `P-score` == min(`P-score`)
                           ) %>%
                           dplyr::filter(`C-score` == max(`C-score`)) %>%
                           dplyr::ungroup()
                        
                     }
                  )
               
               viztools::make_heatmap(
                  isolate(load_heatmap_data()),
                  plotType = isolate(input$heatmap_name),
                  orientation = isolate(input$heatmap_orientation),
                  binSize = isolate(input$heatmap_binsize),
                  massColname = "ObservedPrecursorMass",
                  fractionColname = "fraction",
                  axisRange = isolate(input$heatmap_axisrange),
                  countRange = isolate(input$heatmap_countrange),
                  fillScale = isolate(input$heatmap_fillScale),
                  fontFamily = isolate(input$download_font)
               )
            }
         )
      
      WafflePlotExpr <-
         expr(
            {
               load_waffle_data <-
                  reactive(
                     {
                        readxl::read_xlsx(
                           path =
                              find_newest_file(
                                 fs::path(
                                    tempdir(),
                                    "protein_results_countsbyfraction"
                                 )
                              )
                        ) %>%
                           dplyr::filter(tdreport_name == input$file1) %>%
                           dplyr::select(-tdreport_name)
                     }
                  )
               
               viztools::waffle_iron(
                  isolate(load_waffle_data()),
                  fraction_colname = "fraction",
                  waffleType = "Protein",
                  fontFamily = isolate(input$download_font)
               )
            }
         )
      
      # Plot download handlers -------------------------------------------------------
      
      output$downloadPDF <-
         downloadHandler(
            filename = glue::glue(
               "{format(Sys.time(), '%Y%m%d_%H%M%S')}_{input$plot_type}_plot.pdf"
            ),
            content = function(file) {
               pdf(
                  file = file,
                  width =
                     switch(
                        input$download_unit,
                        inch = input$download_width,
                        cm = input$download_width/2.54
                     ),
                  height =
                     switch(
                        input$download_unit,
                        inch = input$download_height,
                        cm = input$download_height/2.54
                     ),
                  bg = "transparent",
                  useDingbats = FALSE
               )
               print(
                  eval(
                     switch(
                        input$plot_type,
                        upset = eval(UpSetPlotExpr),
                        intdeg = eval(IntDegPlotExpr),
                        heatmap = eval(HeatmapPlotExpr),
                        waffle = eval(WafflePlotExpr)
                     )
                  )
               )
               dev.off()
            }
         )
      
      output$downloadPNG <-
         downloadHandler(
            filename = glue::glue(
               "{format(Sys.time(), '%Y%m%d_%H%M%S')}_{input$plot_type}_plot.png"
            ),
            content = function(file) {
               png(
                  file = file,
                  width =
                     switch(
                        input$download_unit,
                        inch = input$download_width,
                        cm = input$download_width/2.54
                     ),
                  height =
                     switch(
                        input$download_unit,
                        inch = input$download_height,
                        cm = input$download_height/2.54
                     ),
                  units = "in",
                  bg = "white",
                  res = input$download_dpi
               )
               print(
                  eval(
                     switch(
                        input$plot_type,
                        upset = eval(UpSetPlotExpr),
                        intdeg = eval(IntDegPlotExpr),
                        heatmap = eval(HeatmapPlotExpr),
                        waffle = eval(WafflePlotExpr)
                     )
                  )
               )
               dev.off()
            }
         )
      
      output$downloadSVG <-
         downloadHandler(
            filename = glue::glue(
               "{format(Sys.time(), '%Y%m%d_%H%M%S')}_{input$plot_type}_plot.svg"
            ),
            content = function(file) {
               svg(
                  file = file,
                  width =
                     switch(
                        input$download_unit,
                        inch = input$download_width,
                        cm = input$download_width/2.54
                     ),
                  height =
                     switch(
                        input$download_unit,
                        inch = input$download_height,
                        cm = input$download_height/2.54
                     ),
                  bg = "transparent"
               )
               print(
                  eval(
                     switch(
                        input$plot_type,
                        upset = eval(UpSetPlotExpr),
                        intdeg = eval(IntDegPlotExpr),
                        heatmap = eval(HeatmapPlotExpr),
                        waffle = eval(WafflePlotExpr)
                     )
                  )
               )
               dev.off()
            }
         )      
      
      
   }
)
