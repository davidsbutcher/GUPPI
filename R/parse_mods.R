#' parse_mods
#'
#' @param tbl
#' @param modification
#'
#' @return
#' @export
#'
#' @examples

parse_mods <-
   function(
      tbl,
      modification = NULL
   ) {

      if (is.null(modification) == TRUE) {

         stop("No modifications table specified for parse_mods")

      }

      modification_trunc <-
         modification %>%
         dplyr::select(ModificationSetId, Id, Name) %>%
         dplyr::rename(
            "ModSet" = ModificationSetId,
            "ModID" = Id,
            "ModificationName" = Name
         )

      # Get info on modified residues, organize and format them properly

      result1 <-
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
         ) %>%
         dplyr::mutate(
            ModificationName =
               dplyr::case_when(
                  is.na(ModLoc) == FALSE ~
                     paste0(ModificationName, "@", ModLoc)
               )
         ) %>%
         dplyr::group_by(ProteoformRecordNum) %>%
         dplyr::mutate(ModNumber = seq_along(ModificationName)) %>%
         dplyr::mutate(ModNumber = paste("_Mod", ModNumber, sep = "")) %>%
         tidyr::pivot_wider(names_from = ModNumber, values_from = ModificationName) %>%
         dplyr::select(-(ModificationHash:ModLoc)) %>%
         dplyr::summarize_all(coalesce_by_column) %>%
         tidyr::unite("Modification", tidyr::contains("_Mod"), sep = "; ", na.rm = TRUE)

      # Organize and format info on N-terminal mods

      result2 <-
         result1 %>%
         dplyr::left_join(
            modification_trunc,
            by = c(
               "NTerminalModificationSetId" = "ModSet",
               "NTerminalModificationId" = "ModID"
            )
         ) %>%
         dplyr::rename("NTerminalModification" = ModificationName) %>%
         dplyr::left_join(
            modification_trunc,
            by = c(
               "CTerminalModificationSetId" = "ModSet",
               "CTerminalModificationId" = "ModID"
            )
         ) %>%
         dplyr::rename("CTerminalModification" = ModificationName) %>%
         dplyr::select(-(NTerminalModificationSetId:CTerminalModificationId))

      return(result2)
   }
