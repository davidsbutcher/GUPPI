
library(magrittr)
library(GUPPI)
library(viztools)
library(dplyr)
library(assertthat)
library(tictoc)
library(RPushbullet)
library(glue)
library(feather)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(shinyjs)
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

is_local <- Sys.getenv('SHINY_PORT') == ""

# is_local <- FALSE

options(shiny.maxRequestSize = 1000*1024^2)

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

# Server ------------------------------------------------------------------

shinyServer(
   function(input, output, session) {

      options(shiny.maxRequestSize = 1000*1024^2)

      if (is_local == TRUE) hide("input_server")
      if (is_local == FALSE) hide("input_local")

      disable("GUPPIstart")
      disable("VTstart")
      disable("plot_type")
      disable("downloadReport")
      disable("downloadPDF")
      disable("downloadSVG")
      disable("downloadPNG")

      volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

      shinyFileChoose(
         input = input,
         "tdrep_local",
         session = session,
         roots = volumes
      )

      observeEvent(
         input$tdrep,
         {
            output$ULconfirm <-
               renderTable(
                  data.frame(`Uploaded files` = input$tdrep$name)
               )

            enable("GUPPIstart")

            output$outputPlot <-
               NULL
         }
      )

      observeEvent(
         input$tdrep_local,
         {
            output$ULconfirm <-
               renderTable(
                  data.frame(
                     `Uploaded files` = parseFilePaths(volumes, input$tdrep_local)$name
                  )
               )

            enable("GUPPIstart")

            output$outputPlot <-
               NULL
         }
      )

      if (is_local == FALSE) {

         tdrep_path <-
            reactive(
               {
                  input$tdrep$datapath
               }
            )

         tdrep_name <-
            reactive(
               {
                  input$tdrep$name
               }
            )

      } else if (is_local == TRUE) {

         tdrep_path <-
            reactive(
               parseFilePaths(volumes, input$tdrep_local)$datapath
            )

         tdrep_name <-
            reactive(
               parseFilePaths(volumes, input$tdrep_local)$name
            )

      }

      observeEvent(
         input$GUPPIstart,
         {
            outputDir <- tempdir()
            tempReport <- tempfile(fileext = ".html", tmpdir = outputDir)

            purrr::map2_chr(
               tdrep_path(),
               fs::path(
                  dirname(tdrep_path()),
                  tdrep_name()
               ),
               ~file.rename(.x, .y)
            )

            GUPPI::guppi(
               dirname(tdrep_path())[[1]],
               tdrep_name(),
               as.integer(input$taxon),
               outputdir = outputDir,
               fdr = 0.01,
               makeDashboard = T,
               dashboardPath = tempReport,
               saveOutput = T,
               usePB = F
            )

            output$downloadReport <-
               downloadHandler(
                  filename = glue::glue(
                     "{format(Sys.time(), '%Y%m%d_%H%M%S')}_GUPPI_report.html"
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

            updateSelectInput(
               session = session,
               inputId = "file1",
               choices = tdrep_name()
            )

            enable("VTstart")
            enable("plot_type")
            enable("downloadReport")
            enable("downloadPDF")
            enable("downloadSVG")
            enable("downloadPNG")

         }
      )


      # Isolate input params so plot is not created until VTstart is clicked

      isolate(input$file1)
      isolate(input$download_font)

      # Create expression used to generate plots on-demand

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
                  load_UpSet_data(),
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
                  load_intdeg_data(),
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
                           dplyr::group_by(
                              if ((input$heatmap_name) == "Protein") "AccessionNumber" else "ProteoformRecordNum",
                              fraction
                           ) %>%
                           dplyr::filter(`GlobalQvalue` == min(`GlobalQvalue`)) %>%
                           dplyr::filter(`P-score` == min(`P-score`)) %>%
                           dplyr::filter(`C-score` == max(`C-score`)) %>%
                           dplyr::ungroup()

                     }
                  )

               viztools::make_heatmap(
                  load_heatmap_data(),
                  plotType = isolate(input$heatmap_name),
                  orientation = isolate(input$heatmap_orientation),
                  binSize = isolate(input$heatmap_binsize),
                  massColname = "ObservedPrecursorMass",
                  fractionColname = "fraction",
                  axisRange = isolate(input$heatmap_axisrange),
                  countRange = isolate(input$heatmap_countrange),
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
                           dplyr::select(-tdreport_name, -protein_count)
                     }
                  )

               viztools::waffle_iron(
                  load_waffle_data(),
                  fraction_colname = "fraction",
                  waffleType = "Protein",
                  fontFamily = isolate(input$download_font)
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

            if (is_local == FALSE) req(input$tdrep)
            if (is_local == TRUE) req(input$tdrep_local)

            output$confirm <-
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




   }
)
