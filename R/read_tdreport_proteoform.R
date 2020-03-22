#' read_tdreport_proteoform
#'
#' @param tdreport Full path to tdReport file
#' @param fdr_cutoff Cutoff for false detection rate, default 0.01
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples

read_tdreport_proteoform <-
   function(
      tdreport,
      fdr_cutoff = 0.01
   ) {

      message(
         glue::glue("\nEstablishing connection to {basename(tdreport)}...")
      )

      #Establish database connection. Keep trying until it works!

      safe_dbConnect <- purrr::safely(DBI::dbConnect)

      safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)

      if (is.null(safecon[["result"]]) == TRUE) {

         message("\nConnection failed, trying again!")

      }

      iteration_num <- 1

      while (is.null(safecon[["result"]]) == TRUE & iteration_num < 100) {

         iteration_num <- iteration_num + 1

         message(glue("Trying to establish database connection, attempt {iteration_num}"))
         safecon <- safe_dbConnect(
            RSQLite::SQLite(), ":memory:",
            dbname = tdreport,
            synchronous = NULL
         )

      }

      if (is.null(safecon[["result"]]) == TRUE) stop("Failed to connect using SQLite!")

      con <- safecon[["result"]]

      message("\nConnection to tdReport succeeded")

      # Generate SQL query using dbplyr

      protein_accession <-
         dplyr::tbl(con, "Isoform") %>%
         dplyr::left_join(
            dplyr::tbl(con, "GlobalQualitativeConfidence"),
            by = c("Id" = "ExternalId")
         ) %>%
         dplyr::filter(GlobalQvalue <= fdr_cutoff) %>%
         dplyr::collect() %>%
         .$AccessionNumber

      proteoform_results <-
         dplyr::tbl(con, "BiologicalProteoform") %>%
         dplyr::left_join(
            dplyr::tbl(con, "GlobalQualitativeConfidence"),
            by = c("Id" = "ExternalId")
         ) %>%
         dplyr::rename("ExternalId" = Id.x) %>%
         dplyr::left_join(
            dplyr::tbl(con, "ChemicalProteoform"),
            by = c("ChemicalProteoformId" = "Id")
         ) %>%
         dplyr::rename("ProteoformSequence" = Sequence) %>%
         dplyr::filter(GlobalQvalue <= fdr_cutoff) %>%
         dplyr::left_join(
            dplyr::tbl(con, "Isoform"),
            by = c("IsoformId" = "Id")
         ) %>%
         dplyr::rename("IntactSequence" = Sequence) %>%
         dplyr::left_join(
            dplyr::tbl(con, "Hit"),
            by = c("HitId" = "Id")
         ) %>%
         dplyr::left_join(
            dplyr::tbl(con, "ResultSet"),
            by = c("ResultSetId" = "Id")
         ) %>%
         dplyr::left_join(
            dplyr::tbl(con, "DataFile"),
            by = c("DataFileId" = "Id")
         ) %>%
         dplyr::collect() %>%
         dplyr::filter(AccessionNumber %in% protein_accession) %>%
         dplyr::rename(
            "ResultSet" = Name.x, "filename" = Name.y,
            "ChemicalProteoformId" = ChemicalProteoformId.x
         ) %>%
         dplyr::select(
            -Description.x, -Id.y, -Description.y, -ChemicalProteoformId.y,
            -Description, -FilePath
         )

      DBI::dbDisconnect(con)

      return(proteoform_results)

   }
