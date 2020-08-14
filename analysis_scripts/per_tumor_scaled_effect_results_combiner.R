

# script to combine per-tumor scaled effect data 

cores_files <- dir(full.names = T,recursive = T)[grep(pattern = "per_tumor_data.RData",x = dir(recursive = T))][1:4]

tumor_names <- unlist(strsplit(cores_files,split = "/"))[c(F,F,T)]
tumor_names <- unlist(strsplit(tumor_names,split = "_population"))[c(T,F)]


per_tumor_data_main <- NULL
recurrent_var_per_tumor_main <- NULL

library(tidyverse)


for(j in 1:length(cores_files)){
  
  load(cores_files[j])
  
  per_tumor_data <- population_scaled_effect_per_tumor_data$NRSI_prop
  
  per_tumor_data$Variant <- unlist(per_tumor_data$Variant)
  
  
  
  per_tumor_data_var_split <- strsplit(x = per_tumor_data$Variant,split="_")
  recurrent_data_fill <- rep(NA,nrow(per_tumor_data))
  
  coding_ind <- which(sapply(X = per_tumor_data_var_split,FUN = length) == 3)
  
  recurrent_data_fill[coding_ind] <-
    sapply(per_tumor_data_var_split[coding_ind], function(x) paste(x[1:2],collapse = " "))
  
  per_tumor_data$Name <- per_tumor_data$Variant
  per_tumor_data$Name[coding_ind] <- recurrent_data_fill[coding_ind]
  
  per_tumor_data$Gene <- NA
  
  per_tumor_data$Gene[coding_ind] <- sapply(per_tumor_data_var_split[coding_ind], function(x) paste(x[1]))
  
  per_tumor_data$tumor_type <- tumor_names[j]

  per_tumor_data_main <- rbind(per_tumor_data_main,per_tumor_data)
  
  
  recurrent_var_per_tumor <- population_scaled_effect_per_tumor_data$NRSI_prop %>%
    group_by(Tumor) %>%
    summarize(count=n_distinct(Variant))
  
  recurrent_var_per_tumor <- rbind(recurrent_var_per_tumor,
                                   data.frame(
                                     Tumor = population_scaled_effect_per_tumor_data$tumors_without_recurrent_variants, count = 0))
  
  recurrent_var_per_tumor$tumor_type <- tumor_names[j]
  
  recurrent_var_per_tumor_main <- rbind(recurrent_var_per_tumor_main,recurrent_var_per_tumor)
  
  print(tumor_names[j])
}

save(per_tumor_data_main,file = "output_data/combined_per_tumor_proportion_scaled_selection.RData")

save(recurrent_var_per_tumor_main,file = "output_data/recurrent_var_per_tumor.RData")


