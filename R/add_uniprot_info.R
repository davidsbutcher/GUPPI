#' add_uniprot_info
#'
#' @param listofproteins Output from read_tdreport function family.
#' @param database UniProt database to use for adding data.
#' @param tdrep Boolean value, indicates whether the input file is a tdReport.
#'
#' @noRd
#'

add_uniprot_info <-
   function(
      rawProteinList,
      database = NULL,
      tdrep = TRUE
   ) {

      message("Adding UniProt info to proteins")

      # Find column in the input tibble which has "accession"
      # in it and use it to get info from UniProt
      
      if (all(is.na(rawProteinList) == TRUE)) return(NA)
      
      accession_name <-
         grep(
            "accession",
            names(rawProteinList),
            ignore.case = TRUE,
            value = TRUE
         )

      if (tdrep == TRUE) {

         rawProteinList <-
            rawProteinList %>%
            dplyr::rename(
               "UNIPROTKB" := !!accession_name
            )

      } else {

         rawProteinList <-
            tibble(
               UNIPROTKB = pull(rawProteinList, accession_name)
            )

      }

      dplyr::left_join(rawProteinList, database)

   }
