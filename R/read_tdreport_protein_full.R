#' read_tdreport_protein_full
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

read_tdreport_protein_full <-
   function(
      tdreport,
      fdr_cutoff = 0.01
   ) {

      # This function will return ALL hits above Q value threshold, not only the
      # one with the lowest Q value.
      # Output is a tibble with all proteins hits below FDR cutoff
      # Output should match "Hit Report" from TDViewer
      # This version used dbplyr to generate SQL

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

      hitscores <-
         dplyr::tbl(con, "Hit") %>%
         dplyr::rename("HitId" = Id) %>%
         dplyr::select("HitId") %>%
         dplyr::left_join(
            dplyr::tbl(con, "HitScore") %>%
               dplyr::left_join(
                  dplyr::tbl(con, "ScoreType"),
                  by = c("ScoreTypeId" = "Id")
               ) %>%
               dplyr::filter(Name == "P-score") %>%
               dplyr::select("HitId", "Value")
         ) %>%
         dplyr::rename("P-score" = Value) %>%
         dplyr::left_join(
            dplyr::tbl(con, "HitScore") %>%
               dplyr::left_join(
                  dplyr::tbl(con, "ScoreType"),
                  by = c("ScoreTypeId" = "Id")
               ) %>%
               dplyr::filter(Name == "C-score") %>%
               dplyr::select("HitId", "Value")
         ) %>%
         dplyr::rename("C-score" = Value) %>%
         dplyr::collect()

      allproteinhits <-
         dplyr::tbl(con, "Hit") %>%
         dplyr::rename("HitId" = Id) %>%
         dplyr::left_join(dplyr::tbl(con, "GlobalQualitativeConfidence")) %>%
         dplyr::select(
            "HitId", "ObservedPrecursorMass", "ResultSetId",
            "DataFileId", "GlobalQvalue", "ChemicalProteoformId"
         ) %>%
         dplyr::filter(GlobalQvalue <= fdr_cutoff) %>%
         dplyr::left_join(
            dplyr::tbl(con, "DataFile"),
            by = c("DataFileId" = "Id")
         ) %>%
         dplyr::rename("filename" = Name) %>%
         dplyr::select(-c("Description", "FilePath")) %>%
         dplyr::left_join(
            dplyr::tbl(con, "ResultSet"),
            by = c("ResultSetId" = "Id")
         ) %>%
         dplyr::rename("ResultSet" = Name) %>%
         dplyr::left_join(
            dplyr::tbl(con, "ChemicalProteoform"),
            by = c("ChemicalProteoformId" = "Id")
         ) %>%
         dplyr::rename("ProteoformSequence" = Sequence) %>%
         dplyr::left_join(
            dplyr::tbl(con, "BiologicalProteoform"),
            by = c("ChemicalProteoformId")
         ) %>%
         dplyr::rename("BiologicalProteoformId" = Id) %>%
         dplyr::left_join(
            dplyr::tbl(con, "Isoform"),
            by = c("IsoformId" = "Id")
         ) %>%
         dplyr::group_by(HitId) %>%
         dplyr::filter(GlobalQvalue == min(GlobalQvalue)) %>%
         dplyr::collect() %>%
         dplyr::left_join(hitscores) %>%
         dplyr::select(
            "HitId", "ProteoformRecordNum", "AccessionNumber", "GlobalQvalue",
            "P-score", "C-score", "filename", "ObservedPrecursorMass",
            "MonoisotopicMass", "AverageMass", "ProteoformSequence",
            IntactSequence = Sequence, "IsSubsequence", "ResultSet",
            "ModificationHash", "NTerminalModificationSetId",
            "NTerminalModificationId", "CTerminalModificationSetId",
            "CTerminalModificationId"
         ) %>%
         dplyr::ungroup()

      # Close database connection and return output table

      DBI::dbDisconnect(con)

      message("read_tdreport_protein_full Finished!")

      return(allproteinhits)

   }
