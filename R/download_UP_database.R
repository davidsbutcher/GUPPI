
#' download_UP_database
#'
#' @param taxon_number NCBI taxon number of UniProt database to download.
#' @param saveUPdata Logical, controls whether the downloaded taxon info and 
#' database are saved to the package directory.
#'
#' @return
#'
#' @import magrittr
#' 
#' @noRd
#' 


download_UP_database <-
   function(
      taxon_number,
      saveUPdata = TRUE
   ) {
      
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
            paste0("\nTrying to establish database connection, attempt ", iteration_num)
         )
         
         safeUP <- safe_UniProt.ws(taxId = taxon_number)
         
      }
      
      if (is.null(safeUP[["result"]]) == TRUE) stop("\nCould not connect to UniProt")
      
      if (is.null(safeUP[["result"]]) == FALSE) message("\nCONNECTION SUCCEEDED.")
      
      UPtaxon <- safeUP[["result"]]
      
      
      # Make future workers -----------------------------------------------------
      
      if (future::nbrOfWorkers() < 10) {
         
         future::plan(future::multisession(workers = 10))
         
      }
      
      # Download database -------------------------------------------------------
      
      # Use UPtaxon to download the UniProt database
      # Create a safe version of UniProt.ws which will not crash
      # the whole damn program if it fails
      
      numberofchunks <- ceiling(length(UPtaxon@taxIdUniprots)/50)
      
      glue::glue("\nCutting taxon into {numberofchunks} chunks") %>%
         message()
      
      accession_chunks <- chunk2(UPtaxon@taxIdUniprots, numberofchunks)
      
      safeselect <- purrr::safely(UniProt.ws::select)
      
      colsToQuery <-
         c(
            "ENTRY-NAME", "GENES",
            "PROTEIN-NAMES", "ORGANISM",
            "ORGANISM-ID", "SEQUENCE",
            "FUNCTION", "SUBCELLULAR-LOCATIONS",
            "GO-ID"
         )
      
      paste0("\nGetting info from UniProt for taxon ", taxon_number, "\n") %>%
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
         dplyr::bind_rows() %>%
         tibble::as_tibble()
      
      # Output and Cleanup ------------------------------------------------------
      
      if (saveUPdata == TRUE) {
         
         # Save UPtaxon as an RDS
         
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
            fs::path(
               system.file("extdata", package = "GUPPI"),
               "UPdatabase",
               paste0("UPtaxon", taxon_number, ".rds")
            )
         )
         
         # Save UniProt database as an RDS
         
         message("\nWriting UniProt database to rds file\n")
         
         saveRDS(
            UPdatabase,
            fs::path(
               system.file("extdata", package = "GUPPI"),
               "UPdatabase",
               paste0(taxon_number, "_full_UniProt_database.rds")
            )
         )
         
      }
      
      future::plan(future::multisession(workers = 1))
      
      return(UPdatabase)
      
   }

