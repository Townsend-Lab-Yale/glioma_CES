
library(tidyverse)
library(philentropy)

# load("~/Downloads/from_ssh/UCEC__selection_output.RData")
# load("~/Downloads/from_ssh/UCEC_population_scaled_effect_per_tumor_data.RData")


# trinuc_weights <- analysis@trinucleotide_mutation_weights$signatures_output_list
# NRSI_per_tumor <- population_scaled_effect_per_tumor_data$NRSI_prop

divergence_calculation <- function(trinuc_weights, NRSI_per_tumor){
  
  # remove if there are any NA
  
  if(length(which(is.na(NRSI_per_tumor$pct)))>0){
  message("Some NRSI_per_tumor are NA, removing them.")
  NRSI_per_tumor <- NRSI_per_tumor[-which(is.na(NRSI_per_tumor$pct)),]
  }
  
  NRSI_per_tumor_signature_lvl <- NRSI_per_tumor %>%
    group_by(Tumor, Signature) %>%
    summarize(pct_sum_over_sig = sum(pct))
  
  # signature_weights <- lapply(trinuc_weights, function(x) x$weights)
  # signature_weights <- do.call(what = rbind, args = signature_weights)
  signature_weights <- trinuc_weights %>%
    mutate(Tumor = Unique_Patient_Identifier) %>%
    select(-total_snvs,-sig_extraction_snvs, -group_avg_blended,-Unique_Patient_Identifier) 
    
  # signature_weights$Tumor <- rownames(signature_weights)
  signature_weights <- tidyr::pivot_longer(signature_weights,-Tumor,names_to="Signature",values_to="Weight")
  
  # need weights to be a proportion of total and sum to 1 for each tumor
  
  signature_weights <- signature_weights %>%
    group_by(Tumor) %>%
    mutate(Weight_prop = Weight / (sum(Weight))) %>%
    ungroup()
  

  # NRSI_per_tumor_signature_lvl has less unique tumors than signature_weights
  # because NRSI_per_tumor_signature_lvl only has tumors with recurrent variants

  tumor_JSD_df <- data.frame(Tumor = as.character(unique(NRSI_per_tumor_signature_lvl$Tumor)),JSD=NA,stringsAsFactors = F)
    
  for(tumor_ind in 1:length(unique(NRSI_per_tumor_signature_lvl$Tumor))){
    
    NRSIs <- NRSI_per_tumor_signature_lvl[as.character(NRSI_per_tumor_signature_lvl$Tumor) == tumor_JSD_df$Tumor[tumor_ind],]
    sig_weights <- signature_weights[signature_weights$Tumor == tumor_JSD_df$Tumor[tumor_ind],]
    
    merged_NRSIs_sig_weights <- merge(NRSIs, sig_weights)
    
    tumor_JSD_df$JSD[tumor_ind] <- philentropy::JSD(x = rbind(merged_NRSIs_sig_weights$pct_sum_over_sig,
                               merged_NRSIs_sig_weights$Weight_prop))
    
  }
    
  
  return(list(tumor_JSD_df=tumor_JSD_df,
              tumor_signature_weights=signature_weights))
  
}


