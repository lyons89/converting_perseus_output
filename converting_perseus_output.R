
library(tidyverse)
library(openxlsx)


unfiltered = "data/PC1175_VEEV_nonImputed_preFiltered_persesus_output.txt"

unFiltered = read_delim(unfiltered, delim = "\t", skip = 4, show_col_types = FALSE, 
                        col_names = names(read_delim(unfiltered,show_col_types = FALSE)))


imputed = "data/PC1175_VEEV_ttest_persesus_output.txt"

imputedDf = read_delim(imputed, delim = "\t", skip = 4, show_col_types = FALSE, 
                        col_names = names(read_delim(imputed,show_col_types = FALSE)))


Protein = unFiltered %>%
  mutate(Gene = str_trim(
        str_remove_all(
        str_extract_all(Description, "GN=(.*?) ", simplify = TRUE),
        "GN="))) %>%
  dplyr::select(., any_of(c("Accession", "Description", "Gene", "Coverage [%]", "# Peptides",
                            "# Unique Peptides", "# PSMs", "MW [kDa]")), 
                starts_with("Scaled"), starts_with("Norm"))

protein_norm = Protein %>%
  dplyr::select(., Accession, starts_with("Norm"))

Imputed = imputedDf %>%
  mutate(Gene = str_trim(
        str_remove_all(
        str_extract_all(Description, "GN=(.*?) ", simplify = TRUE),
        "GN="))) %>%
  dplyr::select(., any_of(c("Accession", "Description", "Gene", "Coverage [%]", "# Peptides",
                            "# Unique Peptides", "# PSMs", "MW [kDa]")), 
                contains("Difference"), contains("p-value"), contains("q-value"), contains("Significant"),
                starts_with("Norm")) %>%
  mutate(., across(contains("p-value"), ~10^-.x)) %>%
  mutate(., across(contains("Significant"), ~replace(.x, is.na(.x), ""))) %>%
  rename_with(., ~gsub("Norm", "imputed", .x)) %>%
  left_join(., protein_norm, by = "Accession")
  
  
hs = createStyle(textDecoration = "Bold", wrapText = TRUE)


write.xlsx(list(Protein, Imputed), "PC1175_VEEV.xlsx", headerStyle = hs)
