#' add_uniprot_info
#'
#' @param listofproteins
#' @param taxon
#' @param database
#' @param tdrep
#'
#' @return
#' @export
#'
#' @examples

add_uniprot_info <-
   function(
      listofproteins,
      taxon = NULL,
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
