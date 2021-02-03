#' parse_mods_special
#'
#' @noRd


parse_mods_special <-
   function(
      tbl,
      modification = NULL
   ) {

      if (is.null(modification) == TRUE) {

         stop("No modifications table specified for parse_mods")

      }

      modification_trunc <-
         modification %>%
         dplyr::select(ModificationSetId, Id, Name, Formula) %>%
         dplyr::rename(
            "ModSet" = ModificationSetId,
            "ModID" = Id,
            "ModificationName" = Name,
            "ModFormula" = Formula
         ) %>%
         dplyr::mutate(ModFormula = stringr::str_remove_all(ModFormula, " ")) %>%
         dplyr::mutate(
            ModFormula = stringr::str_remove_all(ModFormula, stringr::fixed("+"))
         )

      # Get info on modified residues, organize and format them properly

      result1a <-
         tbl %>%
         tidyr::separate_rows(
            "ModificationHash",
            sep = "\\|"
         ) %>%
         tidyr::separate(
            "ModificationHash",
            c("ModSet", "ModID", "ModLoc"),
            sep = "([\\:\\@])",
            extra = "drop",
            remove = FALSE
         ) %>%
         dplyr::mutate(ModID = as.double(ModID)) %>%
         dplyr::mutate(ModLoc = as.integer(ModLoc)) %>%
         dplyr::mutate(ModLoc = ModLoc + 1) %>%
         dplyr::left_join(
            modification_trunc,
            by = c("ModSet", "ModID")
         )

      result1a %>%
         dplyr::filter(is.na(ModSet) == FALSE) %>%
         dplyr::group_by(ModSet, ModID, ModificationName) %>%
         dplyr::summarize(Count = dplyr::n()) %>%
         dplyr::arrange(dplyr::desc(Count))

   }
