#' parse_mods
#'
#' @noRd

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

      modification_trunc2 <-
         modification %>%
         dplyr::select(ModificationSetId, Id, Name, DiffFormula) %>%
         dplyr::rename(
            "ModSet" = ModificationSetId,
            "ModID" = Id,
            "ModificationName" = Name
         ) %>%
         dplyr::mutate(DiffFormula = stringr::str_remove_all(DiffFormula, " ")) %>%
         dplyr::mutate(
            DiffFormula = stringr::str_remove_all(DiffFormula, stringr::fixed("+"))
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
         ) %>%
         dplyr::mutate(
            ModificationName =
               dplyr::case_when(
                  is.na(ModLoc) == FALSE ~
                     paste0(ModificationName, "@", ModLoc)
               )
         )

      result1b <-
         result1a %>%
         dplyr::group_by(ProteoformRecordNum) %>%
         dplyr::mutate(ModNumber = seq_along(ModificationName)) %>%
         dplyr::mutate(ModNumber = paste("_Mod", ModNumber, sep = "")) %>%
         tidyr::pivot_wider(
            names_from = ModNumber, values_from = c(ModificationName, ModFormula)
         ) %>%
         dplyr::select(-(ModificationHash:ModLoc)) %>%
         dplyr::summarize_all(coalesce_by_column) %>%
         dplyr::left_join(
            modification_trunc2,
            by = c(
               "NTerminalModificationSetId" = "ModSet",
               "NTerminalModificationId" = "ModID"
            )
         ) %>%
         dplyr::rename("ModFormula__ModNterm" = DiffFormula) %>%
         dplyr::select(-ModificationName) %>%
         dplyr::left_join(
            modification_trunc2,
            by = c(
               "CTerminalModificationSetId" = "ModSet",
               "CTerminalModificationId" = "ModID"
            )
         ) %>%
         dplyr::rename("ModFormula__ModCterm" = DiffFormula) %>%
         dplyr::select(-ModificationName) %>%
         tidyr::unite(
            "Modification",
            tidyr::contains("ModificationName__Mod"),
            sep = "; ",
            na.rm = TRUE
         ) %>%
         tidyr::unite(
            "ModFormula",
            tidyr::contains("ModFormula__Mod"),
            sep = "; ",
            na.rm = TRUE
         )

      result1c <-
         result1b %>%
         dplyr::mutate(
            ModFormula_C =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=C)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_H =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=H)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_N =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=N)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_O =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=O)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_P =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=P)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_S =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=S)[0-9]{1,4}"), ~sum(as.integer(.x))
               ),
            ModFormula_Se =
               purrr::map(stringr::str_extract_all(
                  ModFormula, "(?i)(?<=Se)[0-9]{1,4}"), ~sum(as.integer(.x))
               )
         ) %>%
         dplyr::mutate(
            ModFormula =
               glue::glue(
                  "C{ModFormula_C}H{ModFormula_H}N{ModFormula_N}O{ModFormula_O}P{ModFormula_P}S{ModFormula_S}Se{ModFormula_Se}"
               )
         ) %>%
         dplyr::mutate(
            ModFormula =
               stringr::str_remove_all(ModFormula, c("C0|H0|N0|O0|P0|S0|Se0"))
         ) %>%
         dplyr::select(
            -tidyselect::contains("ModFormula_")
         ) %>%
         dplyr::rename("ModificationFormula" = ModFormula)

      # dplyr::mutate(
      #    ModFormula_C =
      #       dplyr::case_when(
      #          stringr::str_detect(ModFormula, "(?i)(?<=C)[0-9]{1,3}") == TRUE ~
      #             stringr::str_extract_all(ModFormula, "(?i)(?<=C)[0-9]{1,3}") %>% unlist(),
      #          TRUE ~ "0"
      #       )
      # ) %>% View

      # Organize and format info on N-terminal & C-terminal mods

      result2 <-
         result1c %>%
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
         dplyr::select(
            -(NTerminalModificationSetId:CTerminalModificationId),
            -ModFormula.x,
            -ModFormula.y
         )

      return(result2)
   }
