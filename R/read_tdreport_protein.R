#' read_tdreport_protein
#'
#' @param tdreport Full path to tdReport file
#' @param fdr_cutoff Cutoff for false detection rate, default 0.01
#'
#' @import magrittr
#'
#' @return
#' @export
#'
#' @examples

read_tdreport_protein <-
   function(
      tdreport,
      fdr_cutoff = 0.01
   ) {

      # Uses dbplyr instead of dplyr. Not much faster than read_tdreport2, possibly
      # has less memory usage

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

      output <-
         dplyr::tbl(con, "Isoform") %>%
         dplyr::left_join(
            dplyr::tbl(con, "GlobalQualitativeConfidence"),
            by = c("Id" = "ExternalId")
         ) %>%
         dplyr::filter(GlobalQvalue <= fdr_cutoff) %>%
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
         dplyr::select(
            "AccessionNumber",
            "GlobalQvalue",
            "ObservedPrecursorMass",
            "Name.y"
         ) %>%
         dplyr::rename("filename" = Name.y)

      # Close database connection and return output table

      DBI::dbDisconnect(con)

      message("read_tdreport_protein Finished!")

      return(output)

   }
