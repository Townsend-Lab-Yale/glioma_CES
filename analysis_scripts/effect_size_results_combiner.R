

# script to combine effect size data for glioma analyses

library(cancereffectsizeR)


cores_files <- dir(full.names = T,recursive = T)[grep(pattern = "__selection_output.RData",x = dir(recursive = T))]

tumor_names <- unlist(strsplit(cores_files,split = "/"))[c(F,F,T)]
tumor_names <- unlist(strsplit(tumor_names,split = "__"))[c(T,F)]

selection_data_main <- NULL
total_substitutions <- vector(mode = "list",length = length(tumor_names))
names(total_substitutions) <- tumor_names


library(tidyverse)


for(j in 1:length(cores_files)){
  
  
  load(cores_files[j])
  
  # selection_data <- cancereffectsizeR::selection_results_converter(analysis)
  maf <- analysis@maf[Variant_Type=="SNV"]
  total_substitutions[[j]] <- table(maf$Unique_Patient_Identifier)
  
  selection_data <- as.data.frame(analysis@selection_results[tumors_with_variant > 1])
  # getting rid of double annotations
  
  annot_number <- sapply(analysis@maf$genes,length)
  annot_index <- which(annot_number > 1)
  if(length(annot_index)>0){
    for(annot in 1:length(annot_index)){
      
      
      # just deduplicate recurrent variants (that is all we have effect size information from) 
      if(length(which(analysis@maf$snv_id==analysis@maf$snv_id[annot_index[annot]]))>1){
        # ZNF321P_N42H_ENSP00000375656
        
        # if the variant is non-coding... 
        if(all(is.na(unlist(analysis@maf$assoc_aa_mut[annot_index[annot]])))){
          
          # pick first gene for the annotation
          
          analysis@maf$genes[annot_index[annot]] <- unlist(analysis@maf$genes[annot_index[annot]])[1]
          
          
        }else{ # if the variant has an associated amino acid change 
          these_variants <- unlist(analysis@maf$assoc_aa_mut[annot_index[annot]])
          # if("19:53472796_T>A" %in% analysis@maf$snv_id[annot_index[annot]]){break}
          # if(unlist(analysis@maf$genes[annot_index[annot]])  == c("ZNF321P", "ZNF816")){break}
          
          selection_ind <- which(selection_data$variant %in% these_variants)
          selection_pick <- which.max(selection_data$selection_intensity[selection_ind])
          winner <- selection_ind[selection_pick]
          
          winner_variant <- selection_data$variant[winner]
          
          winner_maf_ind <- which(analysis@maf$assoc_aa_mut[annot_index[annot]][[1]] == winner_variant)
          
          analysis@maf$assoc_aa_mut[annot_index[annot]][[1]] <- analysis@maf$assoc_aa_mut[annot_index[annot]][[1]][winner_maf_ind]
          analysis@maf$genes[annot_index[annot]][[1]] <- analysis@maf$genes[annot_index[annot]][[1]][winner_maf_ind]
          
          if(length(selection_ind[-selection_pick])>0){
            selection_data <- selection_data[-selection_ind[-selection_pick],]
          }
          
        }
        # if(nrow(selection_data) == 0){break}
      }
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  selection_data_var_split <- strsplit(x = selection_data$variant,split="_")
  recurrent_data_fill <- rep(NA,nrow(selection_data))
  
  coding_ind <- which(sapply(X = selection_data_var_split,FUN = length) == 3)
  
  recurrent_data_fill[coding_ind] <- 
    sapply(selection_data_var_split[coding_ind], function(x) paste(x[1:2],collapse = " "))
  
  selection_data$Name <- selection_data$variant
  selection_data$Name[coding_ind] <- recurrent_data_fill[coding_ind]
  
  selection_data$Gene <- NA
  
  selection_data$Gene[coding_ind] <- sapply(selection_data_var_split[coding_ind], function(x) paste(x[1]))
  # assign gene names to noncoding 
  
  
  # [1] "ZNF321P" "ZNF816"  problem
  
  for(this_row in 1:nrow(selection_data)){
    if(is.na(selection_data$Gene[this_row])){
      # selection_data$Gene <- analysis$mutations$snv[selection_data$variant, genes]
      
      # if(length(unlist(analysis@maf[which(analysis@maf$snv_id == selection_data$variant[this_row])[1],genes]))>1){
      # break
      # }
      # should always be one because we got rid of duplicate annotations above
      selection_data$Gene[this_row] <- unlist(analysis@maf[which(analysis@maf$snv_id == selection_data$variant[this_row])[1],genes])
      
    }
  }
  
  
  
  selection_data$dndscv_q <- NA
  
  for(unique_ind in 1:length(unique(selection_data$Gene))){
    selection_data$dndscv_q[which(selection_data$Gene == unique(selection_data$Gene)[unique_ind])] <- 
      analysis@dndscv_out_list$stageless$sel_cv$qallsubs_cv[which(analysis@dndscv_out_list$stageless$sel_cv$gene_name == unique(selection_data$Gene)[unique_ind])]
  }
  
  
  
  selection_data$tumor_type <- tumor_names[j]
  # dndscv_results <- selection_output$dndscvout$sel_cv
  # rownames(dndscv_results) <- dndscv_results$gene_name
  # 
  # selection_data_df$dndscv_q  <- dndscv_results[selection_data_df$V5,"qallsubs_cv"]
  
  selection_data_main <- rbind(selection_data_main,selection_data)
  
  
  # View(selection_data_df)
  print(tumor_names[j])
}

save(selection_data_main,file = "output_data/combined_selection_results.RData")

save(total_substitutions, file = "output_data/combined_substitution_data.RData")
