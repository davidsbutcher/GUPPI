
#' kickout
#'
#' @param list A list of strings.
#' @param allowed_ext Vector of strings containing allowed file extensions.
#'
#' @return
#'
#' @examples
#'

kickout <-
   function(
      list,
      allowed_ext = c("tdReport", "csv", "xlsx")
   ) {

      # This function removes any element from the list of input files
      # which does not have one of the allowed
      # extensions or which has "deprecated" in its name


      for (i in rev(seq_along(list))) {

         if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {

            list[[i]] <- NULL

         } else if (any(grep("deprecated", list[[i]], fixed = TRUE)) == TRUE) {

            list[[i]] <- NULL

         }
      }

      return(list)
   }


#' get_data_path
#'
#' @param filedir Directory to search for file
#' @param filename Filename to get full path for
#' @param extension Vector of allowed file extensions
#'
#' @importFrom magrittr %>%
#'
#' @return
#' A vector of full file paths
#'
#' @examples
#'

get_data_path <-
   function(
      filedir,
      filename,
      extension
   ) {

      filesindir <-
         fs::dir_ls(
            filedir, recurse = TRUE, type = "file",
            regexp = paste0("[.]", extension, "$")
         ) %>%
         purrr::as_vector()

      if (purrr::map(
         filename,
         ~stringr::str_detect(filesindir, .x)
      ) %>%
      purrr::map(any) %>%
      purrr::as_vector() %>%
      all() == FALSE) stop("One or more input files not found")

      if (filename %>%
          as.list() %>%
          purrr::map(~stringr::str_subset(filesindir, .x)) %>%
          purrr::map(~any(length(.x) > 1)) %>%
          unlist() %>%
          any() == TRUE) stop("One or more input files found in multiple locations")

      filelist <-
         filename %>%
         purrr::map_chr(~stringr::str_subset(filesindir, .x)) %>%
         as.list() %>%
         kickout()

      names(filelist) <- seq(1, length(filelist))

      return(filelist)

   }

#' chunk2
#'
#' @param x Vector
#' @param n Number of chunks to split into
#'
#' @return
#'
#' @examples

chunk2 <- function(x,n) {
   split(x, cut(seq_along(x), n, labels = FALSE))
}

#' get_GO_terms
#'
#' @param tbl
#' @param filelist
#'
#' @return
#'
#' @examples

get_GO_terms <-
   function(tbl, filelist, go_locs) {

      message(
         paste0("Getting GO subcellular locations for ", basename(filelist))
      )

      temptbl <-
         tbl %>%
         tibble::add_column(GO_term = NA) %>%
         tibble::add_column(GO_subcell_loc = NA)

      for (i in seq_along(tbl$`GO-ID`)) {

         temptbl$GO_term[i] <-
            unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
            trimws() %>%
            AnnotationDbi::Term() %>%
            paste(collapse = "; ")

         temptbl$GO_subcell_loc[i] <-
            unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
            trimws() %>%
            AnnotationDbi::Term() %>%
            .[. %in% go_locs] %>%
            paste(collapse = "; ")

      }

      return(temptbl)
   }

#' add_GRAVY
#'
#' @param tbl
#'
#' @return
#'
#' @examples

add_GRAVY <- function(tbl) {

   message("\nAdding hot GRAVY to list of proteins")

   tbl %>%
      dplyr::mutate(
         GRAVY = Peptides::hydrophobicity(SEQUENCE)
      )


}

#' add_GRAVY_allhits
#'
#' @param tbl
#'
#' @return
#'
#' @examples

add_GRAVY_allhits <- function(tbl) {

   message("\nAdding hot GRAVY to list of all protein hits")

   tbl %>%
      dplyr::mutate(
         GRAVY = Peptides::hydrophobicity(ProteoformSequence)
      )


}


#' add_masses
#'
#' @param tbl
#'
#' @return
#'
#' @examples

add_masses <- function(tbl) {

   # This function uses the Peptides package to determine
   # average and monoisotopic masses based on protein sequence.
   # These values diverge from those in the TDreport by <0.00005 Da.

   dplyr::mutate(
      tbl, MonoisotopicMass = Peptides::mw(tbl$SEQUENCE, monoisotopic = TRUE),
      AverageMass = Peptides::mw(tbl$SEQUENCE, monoisotopic = FALSE)
   )

}

#' add_fraction
#'
#' @param tbl
#'
#' @return
#'
#' @examples

add_fraction <- function(tbl) {

   # This function attempts to parse the filenames to extract information
   # about the fraction that a raw file corresponds to. This is only useful
   # for GELFrEE/PEPPI/other fractionated data

   message("\nAdding fraction numbers by parsing filenames")

   tbl %>%
      dplyr::mutate(
         fraction = dplyr::case_when(
            stringr::str_detect(filename,
                                "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}") == TRUE ~
               stringr::str_extract(filename,
                                    "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}"),
            TRUE ~ "NA"
         )
      )

}

#' get_locations_protein
#'
#' @param resultslist
#'
#' @return
#'
#' @examples

get_locations_protein <- function(resultslist) {

   # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
   # GO terms pulled from UniProt for each unique accession number.

   counts <- tibble::tibble(
      filename = basename(names(resultslist)),
      protein_count = NA,
      cytosol_count = NA,
      membrane_count = NA,
      periplasm_count = NA,
      NOTA_count = NA
   )

   for (i in seq_along(resultslist)) {

      tempresults <- tibble::tibble()

      # For every proteoform in each output, get the count of proteoforms
      # whose GO terms include any of
      # "cytosol|cytoplasm|ribosome",
      # "membrane & !membrane-bounded periplasmic space",
      # "periplasm"
      # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!

      tempresults <-
         resultslist[[i]] %>%
         dplyr::mutate(
            cytosol =
               stringr::str_detect(
                  resultslist[[i]]$GO_subcell_loc,
                  c("cytosol|cytoplasm|ribosome")
               )
         ) %>%
         dplyr::mutate(
            membrane =
               stringr::str_detect(resultslist[[i]]$GO_subcell_loc, c("membrane")) &
               !stringr::str_detect(resultslist[[i]]$GO_subcell_loc, c("membrane-bounded periplasmic space"))
         ) %>%
         dplyr::mutate(
            periplasm =
               stringr::str_detect(
                  resultslist[[i]]$GO_subcell_loc,
                  c("periplasm")
               )
         )

      # For cytosol_count and membrane_count, ONLY count the accessions which are NOT
      # found in the list of accessions including both "cytosol|cytoplasm" and "membrane".
      # This prevents double-counting of proteoforms by localization

      counts$protein_count[i] <- length(tempresults$UNIPROTKB)

      counts$cytosol_count[i] <- sum(tempresults$cytosol)

      counts$membrane_count[i] <- sum(tempresults$membrane)

      counts$periplasm_count[i] <- sum(tempresults$periplasm)

      counts$NOTA_count[i] <-
         dplyr::tally(
            tempresults %>%
               dplyr::filter(cytosol == FALSE) %>%
               dplyr::filter(membrane == FALSE) %>%
               dplyr::filter(periplasm == FALSE)
         ) %>%
         .$n

   }

   return(counts)
}

#' get_locations_proteoform
#'
#' @param resultslist
#'
#' @return
#' @export
#'
#' @examples

get_locations_proteoform <- function(resultslist) {

   # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
   # GO terms pulled from UniProt for each unique accession number.

   counts <-
      tibble::tibble(
         filename = basename(names(resultslist)),
         proteoform_count = NA,
         cytosol_count = NA,
         membrane_count = NA,
         periplasm_count = NA,
         NOTA_count = NA
      )

   for (i in seq_along(resultslist)) {

      # For every proteoform in each output, get the count of proteoforms whose GO terms
      # include "cytosol" OR "cytoplasm", "membrane", or BOTH.
      # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!

      tempresults <-
         resultslist[[i]] %>%
         dplyr::mutate(
            cytosol = stringr::str_detect(
               resultslist[[i]]$GO_subcell_loc,
               c("cytosol|cytoplasm|ribosome")
            )
         ) %>%
         dplyr::mutate(
            membrane = stringr::str_detect(
               resultslist[[i]]$GO_subcell_loc,
               c("membrane")
            ) &
               !stringr::str_detect(
                  resultslist[[i]]$GO_subcell_loc,
                  c("membrane-bounded periplasmic space")
               )
         ) %>%
         dplyr::mutate(
            periplasm = stringr::str_detect(
               resultslist[[i]]$GO_subcell_loc,
               c("periplasm")
            )
         )

      # For cytosol_count and membrane_count, ONLY count the accessions which are NOT
      # found in the list of accessions including both "cytosol|cytoplasm" and "membrane".
      # This prevents double-counting of proteoforms by localization

      counts$proteoform_count[i] <- length(tempresults$UNIPROTKB)

      counts$cytosol_count[i] <- sum(tempresults$cytosol)

      counts$membrane_count[i] <- sum(tempresults$membrane)

      counts$periplasm_count[i] <- sum(tempresults$periplasm)

      counts$NOTA_count[i] <-
         dplyr::tally(
            tempresults %>%
               dplyr::filter(cytosol == FALSE) %>%
               dplyr::filter(membrane == FALSE) %>%
               dplyr::filter(periplasm == FALSE)
         ) %>%
         .$n

   }

   return(counts)
}

#' get_locations_byfraction
#'
#' @param resultslist
#'
#' @return
#' @export
#'
#' @examples
#'
get_locations_byfraction <-
   function(resultslist) {

      # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
      # GO terms pulled from UniProt for each unique accession number.

      counts <- tibble::tibble()

      for (i in seq_along(resultslist)) {

         tempresults <- tibble::tibble()

         # For every proteoform in each output, get the count of proteoforms whose GO terms
         # include "cytosol" OR "cytoplasm", "membrane", or BOTH.
         # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!

         tempresults <-
            resultslist[[i]] %>%
            dplyr::mutate(tdreport_name = names(resultslist)[[i]]) %>%
            dplyr::mutate(
               cytosol =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("cytosol|cytoplasm|ribosome")
                  )
            ) %>%
            dplyr::mutate(
               membrane =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane-bounded periplasmic space")
                  )
            ) %>%
            dplyr::mutate(
               periplasm =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("periplasm")
                  )
            ) %>%
            dplyr::mutate(
               NOTA =
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("cytosol|cytoplasm|ribosome")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane-bounded periplasmic space")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("periplasm")
                  )
            )

         tempresultssummary <-
            tempresults %>%
            dplyr::group_by(tdreport_name, fraction) %>%
            dplyr::summarize(
               protein_count = dplyr::n(),
               cytosol_count = sum(cytosol),
               membrane_count = sum(membrane),
               periplasm_count = sum(periplasm),
               NOTA_count = sum(NOTA)
            )

         counts <-
            dplyr::union_all(tempresultssummary, counts)

      }

      return(counts)
   }

#' get_locations_byfraction2
#'
#' @param resultslist
#'
#' @return
#' @export
#'
#' @examples
#'
get_locations_byfraction2 <-
   function(resultslist) {

      # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
      # GO terms pulled from UniProt for each unique accession number.

      # THIS IS AN EXPERIMENTAL VERSION which gets locations based on GLOBAL BEST HITS,
      # NOT PER FRACTION

      counts <- tibble::tibble()

      for (i in seq_along(resultslist)) {

         tempresults <- tibble::tibble()

         # For every proteoform in each output, get the count of proteoforms whose GO terms
         # include "cytosol" OR "cytoplasm", "membrane", or BOTH.
         # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!

         tempresults <-
            resultslist[[i]] %>%
            dplyr::mutate(tdreport_name = names(resultslist)[[i]]) %>%
            dplyr::mutate(
               cytosol =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("cytosol|cytoplasm|ribosome")
                  )
            ) %>%
            dplyr::mutate(
               membrane =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane-bounded periplasmic space")
                  )
            ) %>%
            dplyr::mutate(
               periplasm =
                  stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("periplasm")
                  )
            ) %>%
            dplyr::mutate(
               NOTA =
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("cytosol|cytoplasm|ribosome")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("membrane-bounded periplasmic space")
                  ) &
                  !stringr::str_detect(
                     resultslist[[i]]$GO_subcell_loc,
                     c("periplasm")
                  )
            )

         tempresultssummary <-
            tempresults %>%
            dplyr::group_by(tdreport_name, fraction) %>%
            dplyr::summarize(
               protein_count = dplyr::n(),
               cytosol_count = sum(cytosol),
               membrane_count = sum(membrane),
               periplasm_count = sum(periplasm),
               NOTA_count = sum(NOTA)
            )

         counts <-
            dplyr::union_all(tempresultssummary, counts)

      }

      return(counts)
   }



#' coalesce_by_column
#'
#' @param df
#'
#' @return
#'
#' @examples

coalesce_by_column <- function(df) {

   return(dplyr::coalesce(!!! as.list(df)))

}
