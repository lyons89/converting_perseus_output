
library(tidyverse)
library(openxlsx)


# read in non-Imputed data
pers_NI_loc = "data/PC1159-Braunstein_no1xrev_log2_ttests.txt"

# read in imputed data
pers_I_loc = "data/PC1159-Braunstein_no1xrev_log2_2validvalues1group_impute_ttests.txt"

# read in perseus exported data
# ignore an warning messages

import_pers = function(file_loc){
  
  firstRead = read.delim(file_loc, sep = "\t", na.strings = c("", "NaN", "NA"), check.names = FALSE)
  
  finalRead = read.delim(file_loc, sep = "\t", skip = sum(grepl("#!", firstRead[,1])), 
                         col.names = names(firstRead), na.strings = c("", "NaN", "NA"), check.names = FALSE)
  
}


unImputed = import_pers(file_loc = pers_NI_loc)
imputed = import_pers(file_loc = pers_I_loc)



## Two branching points.
## 1. project type. APMS, Protein-ID, PTM-ID, LFQ
## 2. software. MQ, FP, PD



APMS_MQ = function(df){
  

  untilt = df %>%
    dplyr::select(., any_of(c("Majority protein IDs", "Protein names", "Fasta headers", "Gene names", "Gene name", "Potential contaminant", "Peptides",  
                  "Razor + unique peptides", "Unique peptides")), 
                  contains("Difference"), contains("p-value"), contains("Significant"), starts_with("LFQ intensity"), starts_with("MS/MS count")) %>%
    mutate("Summed LFQ Intensity" = rowSums(2^across(.cols = starts_with("LFQ intensity")), na.rm=TRUE)) %>%
    dplyr::select(., any_of(c("Majority protein IDs", "Protein names", "Fasta headers", "Gene names", "Gene name", "Razor + unique peptides", "Potential contaminant")),
                  contains("Difference"), contains("p-value"), contains("Significant"), starts_with("LFQ intensity"), "Peptides",
                  "Unique peptides", starts_with("Sequence coverage"), "Summed LFQ Intensity", starts_with("MS/MS count"),
                  -contains("significant", ignore.case=FALSE)) %>%
    arrange(.,desc(`Summed LFQ Intensity`))
  
  
}


cleaned_unimputed = APMS_MQ(df = unImputed)
cleaned_imputed = APMS_MQ(df = imputed)


# filter and create signif tabs for each comparison

# creates the additional tabs for each comparison
# df = use imputed data
# filter_type. 2 options, "pval" or "fc_pval"
# "pval" will give same results if you only filtered for Perseus' Student's T-test Significant column.
signif_tabs_MQ = function(df, filter_type){
  
  log_names = names(dplyr::select(df, contains("Difference")))
  pvalue_names = names(dplyr::select(df, contains("p-value")))
  
  if(filter_type == "pval"){
    
    pval_filt = map2(log_names, pvalue_names, function(x,y){
      int = df %>%
        dplyr::filter(., !!as.symbol(y) < 0.05) %>%
        dplyr::arrange(., desc((!!as.symbol(x))))
    })
    return(pval_filt)
  }
  if(filter_type == "fc_pval"){
    
    fc_pval_filt = map2(log_names, pvalue_names, function(x,y){
      
      int = df %>%
        dplyr::filter(., !!as.symbol(x) > 1 | !!as.symbol(y) < -1) %>%
        dplyr::filter(., !!as.symbol(y) < 0.05) %>% # this method to evaluate "column name" as a variable. first turn into symbol, then !! inject it into expression
        dplyr::arrange(., desc(!!as.symbol(x)))
    })
    return(fc_pval_filt)
  }
  
}

sig_tabs = signif_tabs_MQ(df = cleaned_imputed, filter_type = "pval")

# export to excel
hs = createStyle(textDecoration = "Bold", wrapText = TRUE)

tabNames = str_replace_all(str_remove_all(names(dplyr::select(cleaned_imputed, contains("p-value"))), "Student's T-test p-value "), 
                           pattern = "_", replacement = " v ")

write.xlsx(c(list(cleaned_unimputed, cleaned_imputed), sig_tabs), 
           file = paste0(gsub(".*/|.txt", "", pers_I_loc), "_PROCESSED", ".xlsx"), 
           headerStyle = hs,
           sheetName = c("Proteins", "Imputed", tabNames))




##### Don't look below here ######


Protein = unFiltered %>%
  mutate(Gene = str_trim(
        str_remove_all(
        str_extract_all(Description, "GN=(.*?) ", simplify = TRUE),
        "GN="))) %>%
  dplyr::select(., any_of(c("Accession", "Description", "Gene", "Coverage [%]", "# Peptides",
                            "# Unique Peptides", "# PSMs", "MW [kDa]")), 
                starts_with("Scaled"), starts_with("Norm")) %>%
  mutate(., across(starts_with("scaled"), ~replace(., is.nan(.x), 0))) %>%
  mutate(., across(starts_with("Norm"), ~replace(., is.nan(.x), NA))) 

protein_norm = Protein %>%
  dplyr::select(., Accession, starts_with("Norm")) %>%
  mutate(., across(starts_with("Norm"), ~2^.x))

Imputed = imputedDf %>%
  mutate(Gene = str_trim(
        str_remove_all(
        str_extract_all(Description, "GN=(.*?) ", simplify = TRUE),
        "GN="))) %>%
  dplyr::select(., any_of(c("Accession", "Description", "Gene", "Coverage [%]", "# Peptides",
                            "# Unique Peptides", "# PSMs", "MW [kDa]")), 
                contains("Difference"), contains("p-value"), contains("q-value"), contains("Significant"),
                starts_with("Norm")) %>%
  #mutate(., across(contains("p-value"), ~10^-.x)) %>%
  mutate(., across(contains("Significant"), ~replace(.x, is.na(.x), ""))) %>%
  rename_with(., ~gsub("Norm", "imputed", .x)) %>%
  rename_with(., ~gsub("Student's T-test Difference", "log2 fold change:", .x)) %>%
  rename_with(., ~gsub("Student's T-test p-value", "p-value:", .x)) %>%
  rename_with(., ~gsub("Student's T-test q-value ", "q-value:", .x)) %>%
  left_join(., protein_norm, by = "Accession")
  
log_names = names(dplyr::select(Imputed, starts_with("log2")))
pvalue_names = names(dplyr::select(Imputed, starts_with("p-value")))

imputed2 = map2(log_names, pvalue_names, function(y,z){
  
  int = Imputed %>%
    dplyr::filter(., !!as.symbol(y) > 1 | !!as.symbol(y) < -1) %>%
    dplyr::filter(., !!as.symbol(z) < 0.05) %>% # this method to evaluate "column name" as a variable. first turn into symbol, then !! inject it into expression
    dplyr::arrange(., desc(!!as.symbol(y)))
})


hs = createStyle(textDecoration = "Bold", wrapText = TRUE)

tabNames = str_replace_all(str_remove_all(names(dplyr::select(Imputed, starts_with("p-value"))), "p-value: "), pattern = "_", replacement = " v ")

write.xlsx(c(list(Protein, Imputed), imputed2), paste(unfiltered, headerStyle = hs,
           sheetName = c("Proteins", "Imputed", tabNames)))
