
#' guppi
#'
#' @return
#' A variety of information sheets.
#' @export
#' @examples

guppi <-
   function(
      filedir,
      filename,
      taxon_number,
      outputdir,
      fdr = 0.01,
      make_dashboard = FALSE,
      use_PB = FALSE
   ) {

      # Assertions --------------------------------------------------------------

      assertthat::assert_that(
         assertthat::is.dir(filedir),
         msg = "filedir is not a recognized path"
      )

      assertthat::assert_that(
         assertthat::is.string(filename),
         msg = "filename is not a string"
      )

      assertthat::assert_that(
         assertthat::is.number(taxon_number),
         msg = "taxon_number is not a number"
      )

      assertthat::assert_that(
         assertthat::is.dir(outputdir),
         msg = "outputdir is not a recognized path"
      )

      assertthat::assert_that(
         assertthat::is.number(fdr) & fdr > 0 & fdr <= 1,
         msg = "fdr should be a number between 0 and 1"
      )

      assertthat::assert_that(
         assertthat::is.flag(make_dashboard),
         msg = "make_dashboard should be TRUE or FALSE"
      )

      assertthat::assert_that(
         assertthat::is.flag(use_PB),
         msg = "use_PB should be TRUE or FALSE"
      )

      # Get path to data file ---------------------------------------------------

      tictoc::tic()

      filelist <-
         get_data_path(
            filedir,
            filename,
            tools::file_ext(filename)
         )

      # Make future workers -----------------------------------------------------

      # Need to run this command for furrr. If 8 is too many sessions for your
      # system try 5. If running <8 files, change workers to be equal to number
      # of files. workers = 1 is equivalent to not using furrr at all

      # if (length(filelist) >= max_workers) {
      #
      #    # future::plan(future::multisession(workers = as.integer(max_workers)))
      #    future::plan(future::multisession(workers = 8L))
      #
      # } else {
      #
      #    # future::plan(future::multisession(workers = as.integer(length(filelist))))
      #    future::plan(future::multisession(workers = 8L))
      #
      # }


      # Load input files --------------------------------------------------------

      ## check for predownloaded UP database

      if (
         file.exists(
            system.file(
               "extdata",
               paste0(taxon_number, "_full_UniProt_database.feather"),
               package = "GUPPI"
            )
         )
      ) {

         message(
            glue::glue(
               "\nLoading UniProt database for taxon {taxon_number} from package directory"
            )
         )

         UPdatabase <-
            feather::read_feather(
               system.file(
                  "extdata",
                  paste0(taxon_number, "_full_UniProt_database.feather"),
                  package = "GUPPI"
               )
            )

      } else {

         UPdatabase <-
            download_UP_database(taxon_number)

      }

      # Load file containing locations corresponding to
      # GO terms

      go_locs <-
         readRDS(
            system.file(
               "extdata",
               "GO_subcellular_locations.rds",
               package = "GUPPI"
            )
         )

      # Load file containing all the PTM information from TDreport
      # (sourced from UniMod, PSI-MOD, etc.)

      tdreport_mods <-
         feather::read_feather(
            system.file(
               "extdata",
               "tdreport_mods_table.feather",
               package = "GUPPI"
            )
         )

      # Run analysis ---------------------------------------------------------------------------------

      # Save start time to variable for use in output filenames

      systime <- format(Sys.time(), "%Y%m%d_%H%M%S")


      # Read Data Files -----------------------------------------------------------------------------

      extension <-
         filelist %>%
         purrr::map(tools::file_ext)

      if (length(unique(extension)) > 1) {

         stop("More than one kind of input file. Try again.")

      } else if (length(extension) == 0) {

         stop("No acceptable input files. Only tdReport, csv, or xlsx are allowed.")

      } else if (extension[[1]] == "csv") {

         message("Reading csv files...")

         proteinlist  <-
            filelist %>%
            purrr::map(readr::read_csv)

         tdreport_file <- FALSE

      } else if (extension[[1]] == "xlsx") {

         message("Reading xlsx files...")

         proteinlist  <-
            filelist %>%
            purrr::map(readxl::read_xlsx)

         tdreport_file <- FALSE

      } else if (extension[[1]] == "tdReport") {

         message("\nReading protein data from tdReport")

         proteinlist <-
            filelist %>%
            purrr::map(
               read_tdreport_protein,
               fdr_cutoff = fdr
            )


         message("\nReading full protein data from tdReport")

         proteinlistfull <-
            filelist %>%
            purrr::map(
               read_tdreport_protein_full,
               fdr_cutoff = fdr
            )

         # message("\nReading full protein data with spectra from tdReport")
         #
         # proteinlistwithspectra <-
         #    filelist %>%
         #    future_map2(
         #       proteinlistfull,
         #       read_tdreport_withspectra,
         #       fdr_cutoff = fdr
         #    )
         #
         # message("\nReading protein data by file name from tdReport")
         #
         # proteinlistbyfilename <-
         #    filelist %>%
         #    future_map(
         #       read_tdreport_byfilename,
         #       fdr_cutoff = fdr
         #    )
         #
         # message("\nAttempting to read protein data by fraction from tdReport")
         #
         # proteinlistbyfraction <-
         #    filelist %>%
         #    future_map(
         #       read_tdreport_byfraction,
         #       fdr_cutoff = fdr
         #    )

         proteoformlist <-
            filelist %>%
            purrr::map(
               read_tdreport_proteoform,
               fdr_cutoff = fdr
            )

         tdreport_file <- TRUE

      } else {

         stop("No acceptable input files. Only tdReport, csv, or xlsx are allowed.")

      }

      names(proteinlist) <- filelist


      # Add Protein Info --------------------------------------------------------

      # Protein results

      results_protein <-
         proteinlist %>%
         purrr::map(add_uniprot_info,
                    taxon = taxon_number,
                    database = UPdatabase,
                    tdrep = tdreport_file) %>%
         purrr::map2(filelist, get_GO_terms, go_locs) %>%
         purrr::map(add_GRAVY) %>%
         purrr::map(add_masses) %>%
         purrr::map(add_fraction)

      names(results_protein) <-
         unlist(filelist) %>%
         basename()

      results_protein[[length(results_protein)+1]] <-
         get_locations_protein(results_protein)

      names(results_protein)[[length(results_protein)]] <-
         "SUMMARY"

      # Protein results, counts by fraction

      results_protein_countsbyfraction <-
         results_protein[1 : length(results_protein) - 1] %>%
         get_locations_byfraction()

      # Protein results, all hits


      results_protein_allhits <-
         proteinlistfull %>%
         purrr::map(add_GRAVY_allhits) %>%
         purrr::map(add_fraction)

      names(results_protein_allhits) <-
         unlist(filelist) %>%
         basename()


      # Proteform results

      results_proteoform <-
         proteoformlist %>%
         purrr::map(
            add_uniprot_info,
            taxon = 83333,
            database = UPdatabase,
            tdrep = TRUE
         ) %>%
         purrr::map(add_fraction) %>%
         purrr::map(
            parse_mods,
            modification = tdreport_mods
         ) %>%
         purrr::map(
            ~dplyr::mutate(
               .x,
               GRAVY = Peptides::hydrophobicity(ProteoformSequence)
            )
         ) %>%
         purrr::map(
            ~dplyr::select(
               .x,
               -ExternalId,
               -IsEndogenousCleavage,
               -IsoformId,
               -ChemicalProteoformId,
               -AggregationLevel,
               -HitId,
               -IsSubsequence,
               -PriorWeight,
               -EntryId,
               -ScoreForDecoy,
               -ObservedPrecursorMassType,
               -ResultSetId,
               -DataFileId,
               -IsActive,
               -Creator,
               -CreationDate,
               -SEQUENCE
            )
         ) %>%
         purrr::map(
            ~dplyr::select(
               .x,
               UNIPROTKB,
               ProteoformRecordNum,
               IntactSequence,
               ProteoformSequence,
               dplyr::everything()
            )
         ) %>%
         purrr::map2(filelist, get_GO_terms, go_locs)


      names(results_proteoform) <-
         unlist(filelist) %>%
         basename()

      results_proteoform[[length(results_proteoform)+1]] <-
         get_locations_proteoform(results_proteoform)

      names(results_proteoform)[[length(results_proteoform)]] <-
         "SUMMARY"

      # Save results ------------------------------------------------------------

      # Protein results

      for (i in seq_along(names(results_protein))) {

         names(results_protein)[i] <-
            stringr::str_replace_all(names(results_protein[i]), "[:punct:]", "")

         names(results_protein)[i] %<>%
            stringr::str_trunc(28, "left") %>% paste(i, "_", ., sep = "")

      }

      if (dir.exists(glue::glue("{outputdir}/protein_results")) == FALSE) {
         dir.create(glue::glue("{outputdir}/protein_results"))
      }

      resultsname <-
         glue::glue("{outputdir}/protein_results/{systime}_protein_results.xlsx")

      message(
         glue::glue("\nSaving protein results to {resultsname}")
      )

      results_protein %>%
         writexl::write_xlsx(path = resultsname)


      # Protein results, counts by fraction

      if (dir.exists(glue::glue("{outputdir}/protein_results_countsbyfraction")) == FALSE) {
         dir.create(glue::glue("{outputdir}/protein_results_countsbyfraction"))
      }

      resultsname <-
         glue::glue(
            "{outputdir}/protein_results_countsbyfraction/{systime}_protein_countsbyfrac.xlsx"
         )


      results_protein_countsbyfraction %>%
         writexl::write_xlsx(path = resultsname)


      # Protein results, all hits


      if (dir.exists(glue::glue("{outputdir}/protein_results_allhits")) == FALSE) {
         dir.create(glue::glue("{outputdir}/protein_results_allhits"))
      }

      filelist %>%
         purrr::map(basename) %>%
         purrr::map(tools::file_path_sans_ext) %>%
         glue::glue_data("{outputdir}/protein_results_allhits/{.}_allhits.xlsx") %>%
         as.list() %>%
         purrr::walk2(
            results_protein_allhits,
            ~writexl::write_xlsx(.y, path = .x)
         )

      # Proteoform results

      for (i in seq_along(names(results_proteoform))) {

         names(results_proteoform)[i] <-
            stringr::str_replace_all(names(results_proteoform[i]), "[:punct:]", "")

         names(results_proteoform)[i] %<>%
            stringr::str_trunc(28, "left") %>% paste(i, "_", ., sep = "")

      }

      if (dir.exists(glue::glue("{outputdir}/proteoform_results")) == FALSE) {
         dir.create(glue::glue("{outputdir}/proteoform_results"))
      }

      resultsname <-
         glue::glue("{outputdir}/proteoform_results/{systime}_proteoform_results.xlsx")

      message(
         glue::glue("\nSaving proteoform results to {resultsname}")
      )

      results_proteoform %>%
         writexl::write_xlsx(path = resultsname)

      # Make Dashboard ----------------------------------------------------------

      if (make_dashboard == TRUE) {


         if (dir.exists(glue::glue("{outputdir}/report")) == FALSE) {
            dir.create(glue::glue("{outputdir}/report"))
         }

         rmarkdown::render(
            system.file(
               "rmd",
               "generate_dashboard_parent.Rmd",
               package = "GUPPI"
            ),
            output_file =
               glue::glue("{outputdir}/report/{systime}_dashboard.html")
         )

      }

      ## Get run time and message using PushBullet

      totaltime <-
         capture.output(tictoc::toc()) %>%
         stringr::str_extract("[0-9]+") %>%
         as.numeric() %>%
         `/`(60) %>%
         round(digits = 2)

      message(paste0("\nElapsed time: ", totaltime, " min"))

      # Optional line used to contact any Pushbullet enabled device.
      # View ?pbSetup for help

      if (use_PB == TRUE) {

         RPushbullet::pbPost(
            "note", "GUPPI Analysis Finished",
            paste0("Elapsed time: ", totaltime, " min \nOutput Dir: ", outputdir)
         )

      }

      # Save Workspace ---------------------------------------------------------

      # Just in case you want to see an image from a particular run of results

      if (dir.exists(glue::glue("{outputdir}/workspace_image")) == FALSE) {
         dir.create(glue::glue("{outputdir}/workspace_image"))
      }

      save.image(
         glue::glue(
            "{outputdir}/workspace_image/{systime}_workspace_image.RData"
         )
      )


      # Save Session Info ------------------------------------------------------

      # Session info for every run is saved to a txt file in the
      # output directory, in case

      if (dir.exists(glue::glue("{outputdir}/session_info")) == FALSE) {
         dir.create(glue::glue("{outputdir}/session_info"))
      }

      sessioninfo::session_info() %>%
         capture.output() %>%
         writeLines(
            glue::glue("{outputdir}/session_info/{systime}_sessionInfo.txt")
         )

   }
