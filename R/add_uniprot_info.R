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
      listofproteins,
      database = NULL,
      tdrep = TRUE
   ) {

      message("Adding UniProt info to proteins")

      # Find column in the input tibble which has "accession"
      # in it and use it to get info from UniProt

      accession_name <-
         grep(
            "accession",
            names(listofproteins),
            ignore.case = TRUE,
            value = TRUE
         )

      if (tdrep == TRUE) {

         listofproteins <-
            listofproteins %>%
            dplyr::rename(
               "UNIPROTKB" := !!accession_name
            )

      } else {

         listofproteins <-
            tibble(
               UNIPROTKB = pull(listofproteins, accession_name)
            )

      }

      dplyr::left_join(listofproteins, database)

   }
