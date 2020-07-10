
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

options(repos = BiocManager::repositories())

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

      observeEvent(
         input$tdrep,
         {
            output$ULconfirm <-
               renderTable(
                  data.frame(`Uploaded files` = input$tdrep$name)
               )
         }
      )

      observeEvent(
         input$GUPPIstart,
         {
            validate(
               need(input$tdrep, "Upload a tdReport first")
            )

            outputDir <- tempdir()
            tempReport <- tempfile(fileext = ".html", tmpdir = outputDir)

            purrr::map2_chr(
               input$tdrep$datapath,
               fs::path(
                  dirname(input$tdrep$datapath),
                  input$tdrep$name
               ),
               ~file.rename(.x, .y)
            )

            GUPPI::guppi(
               dirname(input$tdrep$datapath)[[1]],
               input$tdrep$name,
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
               choices = input$tdrep$name
            )

         }
      )


      # Isolate input params so plot is not created until startButton is clicked

      isolate(input$file1)
      isolate(input$download_font)

      # Create expression used to generate plots on-demand

      UpSetPlotExpr <-
         expr(
            if (input$upset_name == "Protein"){
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
                  dplyr::select(AccessionNumber, fraction) %>%
                  dplyr::arrange(fraction) %>%
                  tidyr::pivot_wider(
                     names_from = fraction,
                     values_from = AccessionNumber,
                     values_fn = list,
                     names_prefix = "Frac_"
                  ) %>%
                  purrr::flatten() %>%
                  list() %>%
                  rlang::set_names("1") %>%
                  make_UpSet_plot(
                     plotType = input$upset_name,
                     barColor = input$upset_barcolor
                  )
            } else if (input$upset_name == "Proteoform"){
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
                  dplyr::select(ProteoformRecordNum, fraction) %>%
                  dplyr::arrange(fraction) %>%
                  tidyr::pivot_wider(
                     names_from = fraction,
                     values_from = ProteoformRecordNum,
                     values_fn = list,
                     names_prefix = "Frac_"
                  ) %>%
                  purrr::flatten() %>%
                  list() %>%
                  rlang::set_names("1") %>%
                  make_UpSet_plot(
                     plotType = input$upset_name,
                     barColor = input$upset_barcolor
                  )
            }
         )

      IntDegPlotExpr <-
         expr(
            if (input$intdeg_name == "Protein"){
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
                  dplyr::select(AccessionNumber, fraction) %>%
                  dplyr::arrange(fraction) %>%
                  tidyr::pivot_wider(
                     names_from = fraction,
                     values_from = AccessionNumber,
                     values_fn = list,
                     names_prefix = "Frac_"
                  ) %>%
                  purrr::flatten() %>%
                  purrr::map(
                     unique
                  ) %>%
                  make_intersection_degree_plot(
                     Yrange = c(0, as.integer(input$intdeg_yrange)),
                     plotType = input$intdeg_name,
                     fillColor = input$intdeg_fillcolor,
                     fontFamily = input$download_font
                  )
            } else if (input$intdeg_name == "Proteoform"){
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
                  dplyr::select(ProteoformRecordNum, fraction) %>%
                  dplyr::arrange(fraction) %>%
                  tidyr::pivot_wider(
                     names_from = fraction,
                     values_from = ProteoformRecordNum,
                     values_fn = list,
                     names_prefix = "Frac_"
                  ) %>%
                  purrr::flatten() %>%
                  purrr::map(
                     unique
                  ) %>%
                  make_intersection_degree_plot(
                     Yrange = c(0, as.integer(input$intdeg_yrange)),
                     plotType = input$intdeg_name,
                     fillColor = input$intdeg_fillcolor,
                     fontFamily = input$download_font
                  )
            }
         )

      HeatmapPlotExpr <-
         expr(
            if (input$heatmap_name == "Protein"){
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
                  dplyr::group_by(AccessionNumber, fraction) %>%
                  dplyr::filter(`GlobalQvalue` == min(`GlobalQvalue`)) %>%
                  dplyr::filter(`P-score` == min(`P-score`)) %>%
                  dplyr::filter(`C-score` == max(`C-score`)) %>%
                  dplyr::ungroup() %>%
                  viztools::make_heatmap(
                     plotType = input$heatmap_name,
                     orientation = input$heatmap_orientation,
                     binSize = input$heatmap_binsize,
                     massColname = "ObservedPrecursorMass",
                     fractionColname = "fraction",
                     axisRange = input$heatmap_axisrange,
                     countRange = input$heatmap_countrange,
                     fontFamily = input$download_font
                  )
            } else if (input$heatmap_name == "Proteoform"){
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
                  dplyr::group_by(ProteoformRecordNum, fraction) %>%
                  dplyr::filter(`GlobalQvalue` == min(`GlobalQvalue`)) %>%
                  dplyr::filter(`P-score` == min(`P-score`)) %>%
                  dplyr::filter(`C-score` == max(`C-score`)) %>%
                  dplyr::ungroup() %>%
                  viztools::make_heatmap(
                     plotType = input$heatmap_name,
                     orientation = input$heatmap_orientation,
                     binSize = input$heatmap_binsize,
                     massColname = "ObservedPrecursorMass",
                     fractionColname = "fraction",
                     axisRange = input$heatmap_axisrange,
                     countRange = input$heatmap_countrange,
                     fontFamily = input$download_font
                  )
            }
         )

      WafflePlotExpr <-
         expr(
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
               dplyr::select(-tdreport_name, -protein_count) %>%
               viztools::waffle_iron(
                  fraction_colname = "fraction",
                  waffleType = "Protein",
                  fontFamily = input$download_font
               )
         )


      observeEvent(
         input$plot_type,
         {
            updateTabsetPanel(session, "params", selected = input$plot_type)
         }
      )

      observeEvent(
         input$startButton,
         {
            req(input$tdrep)

            output$confirm <-
               NULL

            output$outputPlot <-
               renderPlot(
                  {
                     switch(
                        input$plot_type,
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
