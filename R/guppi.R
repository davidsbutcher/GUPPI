
#' guppi
#'
#' @return
#' A variety of information sheets.
#'
#'
#' @export
#'
#' @examples

guppi <-
   function(
      filedir,
      filename,
      taxonNumber,
      GOLocType = "bacteria",
      fractionAssignments = NULL,
      outputdir,
      fdr = 0.01,
      saveOutput = TRUE,
      makeDashboard = FALSE,
      dashboardPath = glue::glue(
         "{outputdir}/report/{format(Sys.time(), '%Y%m%d_%H%M%S')}_dashboard.html"
      )
   ) {

      # Assertions --------------------------------------------------------------

      assertthat::assert_that(
         assertthat::is.dir(filedir),
         msg = "filedir is not a recognized path"
      )

      purrr::map_chr(
         filename,
         ~assertthat::assert_that(
            assertthat::is.string(.x),
            msg = "A filename is not a string"
         )
      )

      assertthat::assert_that(
         assertthat::is.number(taxonNumber),
         msg = "taxonNumber is not a number"
      )

      assertthat::assert_that(
         GOLocType == "bacteria" | GOLocType == "eukaryota",
         msg = "GOLocType must be 'bacteria' or 'eukaryota'"
      )

      assertthat::assert_that(
         assertthat::is.dir(dirname(outputdir)),
         msg = "outputdir parent is not a recognized path"
      )

      assertthat::assert_that(
         assertthat::is.number(fdr) & fdr > 0 & fdr <= 1,
         msg = "fdr should be a number between 0 and 1"
      )

      assertthat::assert_that(
         assertthat::is.flag(makeDashboard),
         msg = "makeDashboard should be TRUE or FALSE"
      )

      # Get path to data file ---------------------------------------------------

      filelist <-
         get_data_path(
            filedir,
            filename,
            tools::file_ext(filename)
         )

      # Load input files --------------------------------------------------------

      # Start timer

      tictoc::tic()

      # Check for predownloaded UP database

      if (
         file.exists(
            system.file(
               "extdata",
               "UPdatabase",
               paste0(taxonNumber, "_full_UniProt_database.rds"),
               package = "GUPPI"
            )
         )
      ) {

         message(
            glue::glue(
               "\nLoading UniProt database for taxon {taxonNumber} from package directory"
            )
         )

         UPdatabase <-
            readRDS(
               system.file(
                  "extdata",
                  "UPdatabase",
                  paste0(taxonNumber, "_full_UniProt_database.rds"),
                  package = "GUPPI"
               )
            )

      } else {

         UPdatabase <-
            download_UP_database(taxonNumber)

      }

      # Load file containing locations corresponding to
      # GO terms

      # go_locs <-
      #    readRDS(
      #       system.file(
      #          "extdata",
      #          "GO_subcellular_locations.rds",
      #          package = "GUPPI"
      #       )
      #    )

      if (GOLocType == "bacteria") {

         go_locs <-
            read.csv(
               system.file(
                  "extdata",
                  "GO_cellular_component_taxon2_bacteria.csv",
                  package = "GUPPI"
               )
            ) %>%
            dplyr::pull("GO_term")

         go_locs_table <-
            read.csv(
               system.file(
                  "extdata",
                  "GO_cellular_component_taxon2_bacteria.csv",
                  package = "GUPPI"
               )
            )

      } else if (GOLocType == "eukaryota") {

         go_locs <-
            read.csv(
               system.file(
                  "extdata",
                  "GO_cellular_component_taxon2759_eukaryota.csv",
                  package = "GUPPI"
               )
            ) %>%
            dplyr::pull("GO_term")

         go_locs_table <-
            read.csv(
               system.file(
                  "extdata",
                  "GO_cellular_component_taxon2759_eukaryota.csv",
                  package = "GUPPI"
               )
            )

      }

      # Load file containing all the PTM information from TDreport
      # (sourced from UniMod, PSI-MOD, etc.)

      tdreport_mods <-
         readRDS(
            system.file(
               "extdata",
               "tdreport_mods_table.rds",
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

         message("\nReading protein data from tdReport\n")

         proteinlist <-
            filelist %>%
            purrr::map(
               read_tdreport_protein,
               fdr_cutoff = fdr
            )


         message("\nReading full protein data from tdReport\n")

         proteinlistfull <-
            filelist %>%
            purrr::map(
               read_tdreport_protein_full,
               fdr_cutoff = fdr
            )

         message("\nReading proteoform data from tdReport\n")

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

      # Process protein results ------------------------------------------------

      results_protein <-
         proteinlist %>%
         purrr::map(
            add_uniprot_info,
            taxon = taxonNumber,
            database = UPdatabase,
            tdrep = tdreport_file
         ) %>%
         purrr::map2(filelist, get_GO_terms2, go_locs_table) %>%
         purrr::map(add_GRAVY) %>%
         purrr::map(add_masses) %>%
         purrr::map(
            add_fraction,
            assignments = fractionAssignments
         )

      names(results_protein) <-
         unlist(filelist) %>%
         basename()

      results_protein[[length(results_protein)+1]] <-
         get_locations_general(results_protein, go_locs_table)

      names(results_protein)[[length(results_protein)]] <-
         "SUMMARY"

      # Protein results, counts by fraction

      results_protein_countsbyfraction <-
         results_protein[1 : length(results_protein) - 1] %>%
         get_locations_byfraction2(go_locs_table)

      # Protein results, all hits
      # maybe this should be called proteoform allhits?

      results_protein_allhits <-
         proteinlistfull %>%
         purrr::map(add_GRAVY_allhits) %>%
         purrr::map(
            add_fraction,
            assignments = fractionAssignments
         ) %>%
         purrr::map(parse_mods_allhits, modification = tdreport_mods)

      names(results_protein_allhits) <-
         unlist(filelist) %>%
         basename()


      # Proteoform results, best hits per-fraction

      # results_protein_allhits %>%
      #    purrr::map(
      #       ~{
      #          filter(.x, GlobalQvalue <= fdr) %>%
      #             select(-c("HitId")) %>%
      #             select(fraction, everything()) %>%
      #             group_by(fraction, ProteoformRecordNum) %>%
      #             filter(
      #                GlobalQvalue == min(GlobalQvalue) &
      #                   `P-score` == min(`P-score`) &
      #                   `C-score` == max(`C-score`)
      #             ) %>%
      #             distinct() %>%
      #             ungroup()
      #       }
      #    ) %>% .[[1]] %>% View



      # Process proteoform results ----------------------------------------------


      proteoformlist <-
         proteoformlist %>%
         purrr::map(
            ~dplyr::filter(.x, IntactSequence != "DECOY") %>%
               {
                  if (all(.$ProteoformRecordNum) == 0) {
                     dplyr::mutate(
                        .,
                        ProteoformRecordNum =
                           seq_len(length(.$ProteoformRecordNum))
                     )
                  } else {.}
               }
         )

      if (all(proteoformlist$ProteoformRecordNum) == 0) {

         proteoformlist <-
            proteoformlist %>%
            dplyr::mutate(
               ProteoformRecordNum =
                  seq_len(length(proteoformlist$ProteoformRecordNum))
            )

      }

      results_proteoform <-
         proteoformlist %>%
         purrr::map(
            add_uniprot_info,
            taxon = 83333,
            database = UPdatabase,
            tdrep = TRUE
         ) %>%
         purrr::map(
            add_fraction,
            assignments = fractionAssignments
         ) %>%
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
            ~{if (!"IsEndogenousCleavage" %in% names(.x)) dplyr::mutate(.x, IsEndogenousCleavage = NA) else .x} %>%
               {if (!"IsSubsequence" %in% names(.)) dplyr::mutate(., IsSubsequence = NA) else .} %>%
               {if (!"PriorWeight" %in% names(.)) dplyr::mutate(., PriorWeight = NA) else .} %>%
               dplyr::select(
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
         purrr::map2(filelist, get_GO_terms2, go_locs_table)


      names(results_proteoform) <-
         unlist(filelist) %>%
         basename()

      results_proteoform[[length(results_proteoform)+1]] <-
         get_locations_general(results_proteoform, go_locs_table)

      names(results_proteoform)[[length(results_proteoform)]] <-
         "SUMMARY"


      # Save results ------------------------------------------------------------

      if (saveOutput == TRUE) {

         if (dir.exists(outputdir) == FALSE) {
            dir.create(outputdir)
         }

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


      }

      # Make Dashboard ----------------------------------------------------------

      if (makeDashboard == TRUE) {

         if (dir.exists(outputdir) == FALSE) {
            dir.create(outputdir)
         }

         if (dir.exists(glue::glue("{outputdir}/report")) == FALSE) {
            dir.create(glue::glue("{outputdir}/report"))
         }


         # Copy folders to temp output directory and knit there to avoid file
         # permission problems in shiny/shinyapps

         fs::dir_copy(
            system.file(
               "rmd",
               package = "GUPPI"
            ),
            fs::path(
               outputdir,
               "rmd"
            ),
            overwrite = TRUE
         )

         fs::dir_copy(
            system.file(
               "fonts",
               package = "GUPPI"
            ),
            fs::path(
               outputdir,
               "fonts"
            ),
            overwrite = TRUE
         )

         fs::dir_copy(
            system.file(
               "css",
               package = "GUPPI"
            ),
            fs::path(
               outputdir,
               "css"
            ),
            overwrite = TRUE
         )

         rmarkdown::render(
            fs::path(
               outputdir,
               "rmd",
               "generate_dashboard_parent.Rmd"
            ),
            output_file =
               dashboardPath
         )

         fs::dir_delete(
            fs::path(
               outputdir,
               "rmd"
            )
         )

         fs::dir_delete(
            fs::path(
               outputdir,
               "fonts"
            )
         )

         fs::dir_delete(
            fs::path(
               outputdir,
               "css"
            )
         )

      }


      if (saveOutput == TRUE) {

         ## Get run time and message using PushBullet

         totaltime <-
            capture.output(tictoc::toc()) %>%
            stringr::str_extract("[0-9]+") %>%
            as.numeric() %>%
            `/`(60) %>%
            round(digits = 2)

         message(paste0("\nElapsed time: ", totaltime, " min"))

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

   }
