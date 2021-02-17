
#' Connect to a TDReport
#'
#' @param tdreport_path Full path to TDReport file.
#'
#' @return An object of class SQLiteConnection. See ?RSQLite.
#'
#' @examples

connect_tdreport <- function(
   tdreport_path,
   max_attempt = 100L
) {
   
   message(
      glue::glue("\nEstablishing connection to {basename(tdreport_path)}...")
   )
   
   #Establish database connection. Keep trying until it works!
   
   safe_dbConnect <- purrr::safely(DBI::dbConnect)
   
   safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport_path)
   
   if (is.null(safecon[["result"]]) == TRUE) {
      
      message("\nConnection failed, trying again!")
      
   }
   
   iteration_num <- 1
   
   while (is.null(safecon[["result"]]) == TRUE & iteration_num < max_attempt) {
      
      iteration_num <- iteration_num + 1
      
      message(glue("Trying to establish database connection, attempt {iteration_num}"))
      safecon <- safe_dbConnect(
         RSQLite::SQLite(), ":memory:",
         dbname = tdreport_path,
         synchronous = NULL
      )
      
   }
   
   if (is.null(safecon[["result"]]) == TRUE) stop("Failed to connect using SQLite!")
   
   con <- safecon[["result"]]
   
   message("\nConnection to tdReport succeeded")
   
   return(con)
   
}


#' kickout
#'
#' @param list A list of strings.
#' @param allowed_ext Vector of strings containing allowed file extensions.
#'
#' @noRd
#'
#' @return
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
#' @param filename Filenames to get full paths for
#' @param extension Vector of allowed file extensions
#'
#' @noRd
#'
#' @import magrittr
#'
#' @return
#' A vector of full file paths
#'
#'

get_data_path <-
   function(
      filedir,
      filename,
      extension
   ) {
      
      filesindir <-
         fs::dir_ls(
            filedir, recurse = TRUE, type = "file"
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
         as.list()
      
      names(filelist) <- seq(1, length(filelist))
      
      return(filelist)
      
   }

#' chunk2
#'
#' @param x Vector
#' @param n Number of chunks to split into
#'
#' @noRd
#'
#'

chunk2 <- function(x,n) {
   split(x, cut(seq_along(x), n, labels = FALSE))
}

#' get_GO_terms2
#'
#' @param tbl A data frame of protein/proteoform IDs with GO IDs.
#' @param filelist Name of the input file.
#'
#' @noRd
#'
#' @return
#'

get_GO_terms2 <-
   function(tbl, filelist, GO_locs_table) {
      
      if (all(is.na(tbl)) == TRUE) return(NA)
      
      message(
         paste0("Getting GO subcellular locations for ", basename(filelist))
      )
      
      temptbl <-
         tbl %>%
         tibble::add_column(
            GO_term = NA,
            GO_subcell_loc = NA,
            GUPPI_loc = NA
         )
      
      GO_loc_vec <-
         GO_locs_table$GO_ID %>%
         rlang::set_names(GO_locs_table$GUPPI_loc)
      
      for (i in seq_along(tbl$`GO-ID`)) {
         
         temptbl$GO_term[i] <-
            unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
            stringr::str_trim() %>%
            AnnotationDbi::Term() %>%
            paste(collapse = "; ")
         
         temptbl$GO_subcell_loc[i] <-
            unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
            stringr::str_trim() %>%
            AnnotationDbi::Term() %>%
            .[. %in% GO_locs_table$GO_term] %>%
            paste(collapse = "; ")
         
         temptbl$GUPPI_loc[i] <-
            unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
            stringr::str_trim() %>%
            {names(GO_loc_vec)[match(., stringr::str_trim(GO_loc_vec))]} %>%
            {.[is.na(.) == FALSE & nchar(.) != 0]} %>%
            strsplit(";") %>%
            unlist() %>%
            stringr::str_trim() %>%
            unique() %>%
            paste(collapse = "; ")
         
      }
      
      return(temptbl)
   }

#' add_GRAVY
#'
#' @param tbl
#'
#' @noRd
#'
#' @return
#'

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
#' @noRd
#'
#' @return
#'

add_GRAVY_allhits <- function(tbl) {
   
   if (all(is.na(tbl)) == TRUE) return(NA)
   
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
#' @noRd
#'
#' @return
#'

add_masses <- function(sequence, monoisotopic = TRUE) {
   
   # This function uses the Peptides package to determine
   # average and monoisotopic masses based on protein sequence.
   # These values diverge from those in the TDreport by <0.00005 Da.
   
   Peptides::mw(sequence, monoisotopic = monoisotopic)
   
}

#' add_fraction
#'
#' @param tbl
#'
#' @noRd
#'
#' @return
#'

add_fraction <- function(tbl, assignments = NULL) {
   
   # This function attempts to parse the filenames to extract information
   # about the fraction that a raw file corresponds to. This is only useful
   # for GELFrEE/PEPPI/other fractionated data
   
   if (all(is.na(tbl) == TRUE)) return(NA)
   
   if (is.null(assignments) == TRUE) {
      
      message("\nAdding fraction numbers by parsing filenames")
      
      tbl %>%
         {if (!"filename" %in% names(.)) dplyr::mutate(., filename = NA) else .} %>% 
         dplyr::mutate(
            fraction = dplyr::case_when(
               stringr::str_detect(
                  filename,
                  "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}"
               ) == TRUE ~
                  stringr::str_extract(
                     filename,
                     "(?i)(?<=gf|gf_|peppi|peppi_|frac|fraction|f|f_)[0-9]{1,2}"
                  ),
               TRUE ~ "NA"
            )
         )
      
   } else {
      
      message("\nAdding fraction numbers from filename assignments")
      
      assign_tbl <-
         assignments %>%
         tibble::enframe(
            name = "fraction",
            value = "filename"
         ) %>%
         tidyr::unnest(cols = c(filename))
      
      dplyr::left_join(tbl, assign_tbl, by = "filename")
      
   }
   
}

#' get_locations_protein
#'
#' @param resultslist
#'
#' @noRd
#'
#' @return
#'

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
               !stringr::str_detect(
                  resultslist[[i]]$GO_subcell_loc, c("membrane-bounded periplasmic space")
               )
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

#' get_locations_general
#'
#' @param resultslist
#'
#' @noRd
#'
#' @return
#'

get_locations_general <-
   function(resultslist, GO_locs_table) {
      
      GO_loc_vec <-
         unique(GO_locs_table$GUPPI_loc) %>%
         .[. != ""]
      
      location_list <- 
         vector(
            mode = "list",
            length = length(resultslist)
         )
      
      for (i in seq_along(resultslist)) {
         
         if(all(is.na(resultslist[[i]]))) {
            
            location_list[[i]] <- rep(NA, length(GO_loc_vec))
            next()
            
         }
         
         tempresults <-
            purrr::map(
               as.list(resultslist[[i]]$GUPPI_loc),
               ~stringr::str_count(.x, GO_loc_vec)
            )
         
         location_vec <- vector()
         
         for (j in seq_along(GO_loc_vec)) {
            location_vec[[j]] <-
               purrr::map(
                  tempresults,
                  ~.[[j]]
               ) %>%
               unlist() %>%
               sum()
         }
         
         location_list[[i]] <- location_vec
         
      }
      
      names(location_list) <- names(resultslist)
      
      tibble::enframe(location_list) %>%
         tidyr::unnest(cols = c(value)) %>%
         tibble::add_column(loc = rep(GO_loc_vec, length(resultslist))) %>%
         tidyr::pivot_wider(names_from = loc) %>%
         dplyr::rename("results_file_name" = name)
      
   }

#' get_locations_proteoform
#'
#' @param resultslist
#'
#' @noRd
#'
#' @return
#'

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
#' @noRd
#'
#' @return
#'
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
#' @noRd
#'
#' @return
#'

get_locations_byfraction2 <-
   function(resultslist, GO_locs_table) {
      
      GO_loc_vec <-
         unique(GO_locs_table$GUPPI_loc) %>%
         .[. != ""] %>%
         stringr::str_split(";") %>%
         unlist() %>%
         stringr::str_trim()
      
      get_count <-
         function(vec, pattern) {
            stringr::str_count(vec, pattern) %>%
               sum()
         }
      
      resultslist_TEMP <-
         resultslist %>%
         purrr::imap(
            ~dplyr::mutate(.x, results_file_name = .y)
         ) %>%
         purrr::reduce(dplyr::union_all)
      
      resultslist_TEMP3 <-
         resultslist_TEMP %>%
         dplyr::group_by(results_file_name, fraction) %>%
         dplyr::summarize()
      
      for (i in seq_along(GO_loc_vec)) {
         
         resultslist_TEMP2 <-
            resultslist_TEMP %>%
            dplyr::group_by(results_file_name, fraction) %>%
            dplyr::summarize(get_count(GUPPI_loc, GO_loc_vec[[i]])) %>%
            dplyr::rename(!!rlang::sym(GO_loc_vec[[i]]) := dplyr::last_col())
         
         resultslist_TEMP3 <-
            dplyr::left_join(resultslist_TEMP3, resultslist_TEMP2)
         
      }
      
      return(resultslist_TEMP3)
      
   }

#' get_locations_byfraction_exp
#'
#' @param resultslist
#'
#' @noRd
#'
#' @return
#'

get_locations_byfraction_exp <-
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
#' @noRd
#'
#' @return
#'

coalesce_by_column <- function(df) {
   
   return(dplyr::coalesce(!!! as.list(df)))
   
}



get_result_parameters <-
   function(
      tdreport
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
      
      output <-
         dplyr::tbl(con, "ResultParameter") %>%
         dplyr::left_join(
            dplyr::tbl(con, "ResultSet"),
            by = c("ResultSetId" = "Id")
         ) %>%
         dplyr::collect() %>%
         dplyr::select(
            -tidyselect::any_of(
               c(
                  "Id", "ResultSetId"
               )
            ),
            tidyselect::any_of(
               c(
                  "GroupName", "Name.x", "Value", "Name.y"
               )
            )
         ) %>%
         dplyr::rename(
            "Name" = tidyselect::any_of("Name.x"),
            "ResultSet" = tidyselect::any_of("Name.y")
         )
      
      
      # Close database connection and return output table
      
      DBI::dbDisconnect(con)
      
      message("get_result_parameters Finished!")
      
      return(output)
      
      
   }


#' Determine top-down software suite
#'
#' @param resultsFile A results file taken from top-down analysis by TopPIC,
#' MSPathfinder, etc.
#'
#' @return String indicating the software suite the results file came from.
#'
#'

determine_software <- 
   function(
      resultsFile
   ){
      
      # Assertions --------------------------------------------------------------
      
      assertthat::assert_that(
         assertthat::is.readable(resultsFile),
         msg = "resultsFile is not a readable file"
      ) 
      
      # Try to determine software -----------------------------------------------
      
      # Check if the extension if tdReport
      
      if (tolower(fs::path_ext(resultsFile)) == "tdreport") return("tdreport")
      
      # Read the top 50 lines of the results. This should be more than enough
      # to determine identity, and saves time
      
      lines <- 
         readr::read_lines(
            resultsFile,
            n_max = 50
         )
      
      # Check for column names found in MSPathfinder output
      
      mspath_colnames <- 
         c(
            "Scan", "Pre", "Sequence", "Post", "Modifications", "Composition",
            "ProteinName", "ProteinDesc", "ProteinLength", "Start", "End",
            "Charge", "MostAbundantIsotopeMz", "Mass", "Ms1Features",
            "#MatchedFragments", "Probability", "SpecEValue", "EValue",
            "QValue", "PepQValue"
         )
      
      results_colnames <- 
         stringr::str_split(lines[1], "\t") %>% 
         purrr::as_vector()
      
      if (all(mspath_colnames %in% results_colnames)) return("mspathfinder")
      
      
      # Check for presence of "** Parameters **" line found in TopPIC output
      
      toppic_format <- 
         stringr::str_detect(lines, stringr::fixed("** Parameters **")) %>% 
         any()
      
      if (toppic_format == TRUE) return("toppic")
      
      # Check to see if its just a list of accessions
      
      resultsFileSpread <- 
         switch(
            fs::path_ext(resultsFile),
            "xlsx" = readxl::read_xlsx(resultsFile),
            "tsv" = readr::read_tsv(resultsFile),
            "csv" = readr::read_csv(resultsFile)
         )
      
      accessionlist_format <- 
         all(
            c(
               length(resultsFile) == 1,
               dplyr::pull(resultsFileSpread, 1) %>% 
                  typeof() == "character"
            )
         )
      
      if (length(resultsFile) == 1) return("accessionlist")
      
      # If all tests above fail, return "unknown"
      
      return("unknown")   
      
   }

toppic_lines_to_skip <- 
   function(
      resultsFile
   ){
      
      # Assertions --------------------------------------------------------------
      
      assertthat::assert_that(
         assertthat::is.readable(resultsFile),
         msg = "resultsFile is not a readable file"
      ) 
      
      # Find last line where "** Parameters **" occurs
      
      readr::read_lines(
         resultsFile,
         n_max = 50
      ) %>% 
         stringr::str_detect(stringr::fixed("** Parameters **")) %>% 
         which() %>% 
         dplyr::last()
      
      
   }


reshape_toppic_protein <- 
   function(
      rawProteinList
   ) {
      
      rawProteinList %>% 
         tidyr::separate(
            `Protein name`,
            sep = "\\|",
            into = c("Protein name.1", "UNIPROTKB", "Protein name.3")
         ) %>% 
         dplyr::select(
            -c(
               'Protein name.1',
               'Protein name.3'
            )
         ) %>% 
         dplyr::select(
            UNIPROTKB,
            dplyr::everything()
         ) %>% 
         dplyr::rename(
            'filename' = dplyr::any_of('Data file name'),
            'ObservedPrecursorMass' = dplyr::any_of('Precursor mass'),
            'ProteoformSequencewithMods' = dplyr::any_of('Proteoform'),
            'Pvalue' = dplyr::any_of('`P-value`'),
            'Evalue' = dplyr::any_of('`E-value`'),
            'Qvalue' = dplyr::any_of('`Q-value (spectral FDR)`')
         ) %>% 
         dplyr::group_by(UNIPROTKB) %>% 
         dplyr::arrange(
            `Proteoform FDR`,
            `E-value`,
            `P-value`,
            .by_group = TRUE
         ) %>% 
         dplyr::slice_head()
      
   }

reshape_mspathfinder_protein <- 
   function(
      rawProteinList
   ) {
      
      rawProteinList %>% 
         tidyr::separate(
            ProteinName,
            sep = "\\|",
            into = c("ProteinName.1", "UNIPROTKB", "ProteinName.3")
         ) %>% 
         dplyr::select(
            -c(
               'ProteinName.1',
               'ProteinName.3',
               Sequence
            )
         ) %>% 
         dplyr::select(UNIPROTKB, dplyr::everything()) %>% 
         dplyr::rename(
            'Evalue' = dplyr::any_of('EValue'),
            'Qvalue' = dplyr::any_of('QValue')
         ) %>% 
         dplyr::mutate(filename = "NA") %>% 
         dplyr::group_by(UNIPROTKB) %>% 
         dplyr::arrange(
            PepQValue,
            Qvalue,
            Evalue,
            SpecEValue,
            .by_group = TRUE
         ) %>% 
         dplyr::slice_head()
      
   }

reshape_toppic_proteoform <- 
   function(
      rawProteoformList
   ) {
      
      rawProteoformList %>% 
         tidyr::separate(
            `Protein name`,
            sep = "\\|",
            into = c("Protein name.1", "UNIPROTKB", "Protein name.3")
         ) %>% 
         dplyr::select(
            -c(
               'Protein name.1',
               'Protein name.3'
            )
         ) %>% 
         dplyr::select(
            UNIPROTKB,
            dplyr::everything()
         ) %>% 
         dplyr::rename(
            'filename' = dplyr::any_of('Data file name'),
            'ObservedPrecursorMass' = dplyr::any_of('Precursor mass'),
            'ProteoformSequencewithMods' = dplyr::any_of('Proteoform'),
            'Pvalue' = dplyr::any_of('`P-value`'),
            'Evalue' = dplyr::any_of('`E-value`'),
            'Qvalue' = dplyr::any_of('`Q-value (spectral FDR)`')
         )
      
   }

reshape_mspathfinder_proteoform <- 
   function(
      rawProteoformList
   ) {
      
      rawProteoformList %>% 
         tidyr::separate(
            ProteinName,
            sep = "\\|",
            into = c("ProteinName.1", "UNIPROTKB", "ProteinName.3")
         ) %>% 
         dplyr::select(
            -c(
               'ProteinName.1',
               'ProteinName.3'
            )
         ) %>% 
         dplyr::select(UNIPROTKB, dplyr::everything()) %>% 
         dplyr::rename(
            'Evalue' = dplyr::any_of('EValue'),
            'Qvalue' = dplyr::any_of('QValue')
         ) %>% 
         dplyr::mutate(filename = "NA") 
      
   }

replace_missing_PFR <- 
   function(
      tdrep_proteoforms
   ) {
      
      tdrep_proteoforms %>% 
         dplyr::filter(IntactSequence != "DECOY") %>% 
         {
            if (all(.$ProteoformRecordNum) == 0) {
               dplyr::mutate(
                  .,
                  ProteoformRecordNum =
                     seq_len(length(.$ProteoformRecordNum))
               )
            } else {.}
         }
      
   }
