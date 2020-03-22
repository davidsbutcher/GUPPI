
#' download_UP_database
#'
#' @param taxon_number
#'
#' @return
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples

download_UP_database <-
   function(
      taxon_number,
      max_workers = 8L
   ) {

      # Check for a file in the input folder which contains the UniProt.ws
      # taxon data. If it exists, load it. Otherwise, download it SAFELY
      # trying up to 10 times

      if (
         file.exists(
            system.file(
               "extdata",
               glue::glue("UPtaxon{taxon_number}.rds"),
               package = "GUPPI"
            )
         )
      ) {

         message("\nFound corresponding RDS file. Loading...")

         UPtaxon <-
            readRDS(
               system.file(
                  "extdata",
                  glue::glue("UPtaxon{taxon_number}.rds"),
                  package = "GUPPI"
               )
            )

      } else {

         # Establish connection to UniProt WS, safely!

         safe_UniProt.ws <-
            purrr::safely(UniProt.ws::UniProt.ws)

         message("\nTrying to connect to UniProt web service...")

         safeUP <- safe_UniProt.ws(taxId = taxon_number)

         if (is.null(safeUP[["result"]]) == TRUE) message("\nConnection failed, trying again!\n")

         iteration_num <- 1

         while (is.null(safeUP[["result"]]) == TRUE & iteration_num < 11) {

            iteration_num <- iteration_num + 1

            message(
               paste0("Trying to establish database connection, attempt ", iteration_num)
            )

            safeUP <- safe_UniProt.ws(taxId = taxon_number)

         }

         if (is.null(safeUP[["result"]]) == TRUE) stop("\nCould not connect to UniProt")

         if (is.null(safeUP[["result"]]) == FALSE) message("\nCONNECTION SUCCEEDED.")

         UPtaxon <- safeUP[["result"]]

      }


      # Make future workers -----------------------------------------------------

      if (future::nbrOfWorkers() == 1) {

         future::plan(future::multisession(workers = 6))

      }

      # Download database -------------------------------------------------------

      # Use UPtaxon to download the UniProt database
      # Create a safe version of UniProt.ws which will not crash
      # the whole damn program if it fails

      numberofchunks <- ceiling(length(UPtaxon@taxIdUniprots)/100)

      accession_chunks <- chunk2(UPtaxon@taxIdUniprots, numberofchunks)

      safeselect <- purrr::safely(UniProt.ws::select)

      # message(
      #    paste0(
      #       "There are ",
      #       length(accession_chunks),
      #       " chunks"
      #    )
      # )

      colsToQuery <-
         c(
            "ENTRY-NAME", "GENES",
            "PROTEIN-NAMES", "ORGANISM",
            "ORGANISM-ID", "SEQUENCE",
            "FUNCTION", "SUBCELLULAR-LOCATIONS",
            "GO-ID"
         )

      paste0("\nGetting info from UniProt for taxon ", taxon_number) %>%
         message()

      results_safe <-
         furrr::future_map(
            accession_chunks,
            ~safeselect(
               UPtaxon,
               keys = .x,
               columns = colsToQuery,
               keytype = "UNIPROTKB"
            ),
            .progress = TRUE
         )

      results_safe %>%
         purrr::walk(
            ~(
               if (is.null(.x[["result"]]) == TRUE) {
                  stop("NULL result found while downloading UP database")
               }
            )
         )

      UPdatabase <-
         results_safe %>%
         purrr::map(~(.x[["result"]])) %>%
         purrr::reduce(dplyr::union_all) %>%
         tibble::as_tibble()

      # Output and Cleanup ------------------------------------------------------

      message(
         paste0(
            system.file("extdata", package = "GUPPI"),
            "/UPtaxon",
            taxon_number,
            ".rds"
         )
      )

      saveRDS(
         UPtaxon,
         file =
            paste0(
               system.file("extdata", package = "GUPPI"),
               "/UPtaxon",
               taxon_number,
               ".rds"
            )
      )

      feather::write_feather(
         UPdatabase,
         paste0(
            system.file("extdata", package = "GUPPI"),
            "/",
            taxon_number,
            "_full_UniProt_database.feather"
         )
      )

      return(UPdatabase)

   }

