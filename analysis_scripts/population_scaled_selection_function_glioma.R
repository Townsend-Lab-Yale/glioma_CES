

net_realized_selection_calculation <- function(selection_output, 
                                               selection_output_df=NULL,
                                               signatures="signatures_cosmic_v3", 
                                               subset_freq_greater_than=1, 
                                               substitution_tumor_cutoff=-1, 
                                               median_calc=F,
                                               cores_for_contribution_mat=NULL){
  # signatures.cosmic found in the deconstructSigs package
  
  if(signatures == "signatures_cosmic_v3"){
    data("signatures_cosmic_May2019",package = "cancereffectsizeR")
    signatures <- signatures_cosmic_May2019
  }
  
  # convert selection output to usable format
  if(is.null(selection_output_df)){
    selection_data <- as.data.frame(selection_output@selection_results[tumors_with_variant > 0])
  }else{
    selection_data <- selection_output_df
  }
  
  contribution_matrix_list <- vector(mode = "list",length = nrow(selection_data))
  names(contribution_matrix_list) <- selection_data$variant
  
  
  # Giving the MAF a consistent label with selection output
  # library(data.table)
  selection_output_df <- as.data.frame(selection_output@maf)
  
  annot_number <- sapply(selection_output_df$assoc_aa_mut,length)
  annot_index <- which(annot_number > 1)
  if(length(annot_index)>0){
    for(annot in 1:length(annot_index)){
      # if(annot == 88){break}
      
      # if(!"ABHD14A-ACY1_L354P_ENSP00000420487" %in% selection_data$variant){break}
      
      these_variants <- unlist(selection_output_df$assoc_aa_mut[annot_index[annot]])
      
      # if("ABHD14A-ACY1_L354P_ENSP00000420487" %in% these_variants){break}
      
      selection_ind <- which(selection_data$variant %in% these_variants)
      selection_pick <- which.max(selection_data$selection_intensity[selection_ind])
      winner <- selection_ind[selection_pick]
      
      winner_variant <- selection_data$variant[winner]
      
      winner_maf_ind <- which(selection_output_df$assoc_aa_mut[annot_index[annot]][[1]] == winner_variant)
      
      selection_output_df$assoc_aa_mut[annot_index[annot]][[1]] <- selection_output_df$assoc_aa_mut[annot_index[annot]][[1]][winner_maf_ind]
      selection_output_df$genes[annot_index[annot]][[1]] <- selection_output_df$genes[annot_index[annot]][[1]][winner_maf_ind]
      
      if(length(selection_ind[-selection_pick])>0){
        selection_data <- selection_data[-selection_ind[-selection_pick],]
      }
      # if(nrow(selection_data) == 0){break}
      
      
    }
  }
  
  
  # adding variant matcher
  
  selection_output_df$variant_matcher <- selection_output_df$assoc_aa_mut
  selection_output_df$variant_matcher[is.na(selection_output_df$variant_matcher)] <- 
    selection_output_df$snv_id[is.na(selection_output_df$variant_matcher)]
  
  # selection_output_df$label <- selection_output_df$unique_variant_ID
  # selection_output_df$label[which(sapply(strsplit(selection_output_df$unique_variant_ID_AA,split=" "), function(x) length(x)==1))] <- paste(selection_output_df$Gene_name[which(sapply(strsplit(selection_output_df$unique_variant_ID_AA,split=" "), function(x) length(x)==1))], selection_output_df$unique_variant_ID_AA[which(sapply(strsplit(selection_output_df$unique_variant_ID_AA,split=" "), function(x) length(x)==1))])
  # 
  
  NRSI_matrix <- matrix(data = NA,nrow = nrow(selection_data),ncol = nrow(signatures))
  rownames(NRSI_matrix) <- selection_data$variant
  colnames(NRSI_matrix) <- rownames(signatures)
  
  
  # contribution matrix function
  
  #' Contribution matrix stores: for each substitution, for each tumor, for each signature
  #' the contribution of the total flux of that variant's trinucleotide context from each signature
  #' in that tumor.
  
  mutational_weights <- as.data.frame(selection_output@trinucleotide_mutation_weights$signature_weights_table)
  
  contribution_mat_function <- function(sub_index,median_calc){
    
    
    this_annot_maf <- selection_output_df[which(selection_output_df$variant_matcher==selection_data$variant[sub_index]),]
    
    contribution_matrix <- matrix(data = 0,
                                  nrow = nrow(this_annot_maf),
                                  ncol = nrow(signatures))
    rownames(contribution_matrix) <- this_annot_maf$Unique_Patient_Identifier
    colnames(contribution_matrix) <- rownames(signatures)
    
    
    for(this_tumor_index in 1:nrow(contribution_matrix)){
      
      # if there is no substitution count data structure, 
      # then there was zero substitutions.
      # if(is.null(selection_output@trinucleotide_mutation_weights$signatures_output_list[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count)){
      #   selection_output@trinucleotide_mutation_weights$signatures_output_list[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count <- 0
      # }
      # 
      
      # selection_output@trinucleotide_mutation_weights$signature_weights_table[Unique_Patient_Identifier==this_annot_maf$Unique_Patient_Identifier[this_tumor_index],]
      
      #TODO: remove after JM fixes zero substitution weight
      # if(!any(mutational_weights$Unique_Patient_Identifier== rownames(contribution_matrix)[this_tumor_index])){
      #   next
      # }

      # if(selection_output@trinucleotide_mutation_weights$signatures_output_list[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count > substitution_tumor_cutoff) {
      
      for(this_signature_index in 1:ncol(contribution_matrix)){
        
        contribution_matrix[this_tumor_index,this_signature_index] <- as.numeric(
          (mutational_weights[
            mutational_weights$Unique_Patient_Identifier== rownames(contribution_matrix)[this_tumor_index]
            ,rownames(signatures)[this_signature_index]] *
             signatures[this_signature_index,this_annot_maf[this_tumor_index,"trinuc_mut"]]) /selection_output@trinucleotide_mutation_weights$trinuc_proportion_matrix[rownames(contribution_matrix)[this_tumor_index],this_annot_maf[this_tumor_index,"trinuc_mut"]]
        )
        
        
      }
      # }else{
      #   contribution_matrix[this_tumor_index,] <- 0
      #   
      # }
      
    }
    
    
    if(median_calc){
      median_vec <- apply(contribution_matrix,2,median)
      
      NRSI_vec <- median_vec * selection_data[sub_index,"selection_intensity"] * (selection_data[sub_index,"variant_freq"]/ length(unique(selection_output$MAF$Unique_patient_identifier)))
    }else{
      
      NRSI_vec <- selection_data[sub_index,"selection_intensity"] * (apply(contribution_matrix,2,sum)/length(unique(selection_output_df$Unique_Patient_Identifier)))
      
      
    }
    
    
    
    
    
    
    return(list(contribution_matrix = contribution_matrix,
                NRSI_vec = NRSI_vec))
  }
  
  
  
  if(is.null(cores_for_contribution_mat)){
    cores_for_contribution_mat <- parallel::detectCores()
  }
  
  contribution_mat_function_output <- parallel::mclapply(X = 1:nrow(selection_data),
                                                         FUN = contribution_mat_function,
                                                         mc.cores = cores_for_contribution_mat,
                                                         median_calc=median_calc)
  
  # contribution_mat_function_output <- pbapply::pblapply(X = 1:nrow(selection_data),
  #                                                        FUN = contribution_mat_function,
  #                                                        cl = cores_for_contribution_mat,
  #                                                        median_calc=median_calc)
  
  # out_list <- list()
  # for (sub_index in 1:nrow(selection_data)){
  #   out_list[[sub_index]] <- contribution_mat_function(sub_index = sub_index,median_calc = F)
  # }
  # 
  # 79 breaks
  
  # 
  # # calculate net realized selection intensity for each substitution
  # for(sub_index in 1:nrow(selection_data)){
  #   contribution_matrix <- matrix(data = 0,
  #                                 nrow = selection_data$variant_freq[sub_index],
  #                                 ncol = nrow(signatures))
  #   rownames(contribution_matrix) <- selection_output$MAF[which(selection_output$MAF$label==selection_data$variant[sub_index]),"Unique_patient_identifier"]
  #   colnames(contribution_matrix) <- rownames(signatures)
  #   
  #   
  #   for(this_tumor_index in 1:nrow(contribution_matrix)){
  #     
  #     # if there is no substitution count data structure, 
  #     # then there was zero substitutions.
  #     if(is.null(selection_output$trinuc_data$signatures_output[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count)){
  #       selection_output$trinuc_data$signatures_output[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count <- 0
  #     }
  #     
  #     if(selection_output$trinuc_data$signatures_output[[rownames(contribution_matrix)[this_tumor_index]]]$substitution_count > substitution_tumor_cutoff) {
  #       
  #       for(this_signature_index in 1:ncol(contribution_matrix)){
  #         
  #         contribution_matrix[this_tumor_index,this_signature_index] <- as.numeric(
  #           (selection_output$trinuc_data$signatures_output[[
  #             rownames(contribution_matrix)[this_tumor_index]
  #             ]]$signatures_output$weights[this_signature_index] *
  #              signatures[this_signature_index,selection_output$MAF[which(selection_output$MAF$label==selection_data$variant[sub_index])[this_tumor_index],"trinuc_dcS"]]) /selection_output$trinuc_data$trinuc_proportion_matrix[rownames(contribution_matrix)[this_tumor_index],selection_output$MAF[which(selection_output$MAF$label==selection_data$variant[sub_index])[this_tumor_index],"trinuc_dcS"]]
  #         )
  #         
  #         
  #       }
  #     }else{
  #       contribution_matrix[this_tumor_index,] <- 0
  #       
  #     }
  #     
  #   }
  #   
  #   # store for later
  #   contribution_matrix_list[[sub_index]] <- contribution_matrix
  #   
  #   if(median_calc){
  #     median_vec <- apply(contribution_matrix,2,median)
  #     
  #     NRSI_vec <- median_vec * selection_data[sub_index,"selection_intensity"] * (selection_data[sub_index,"variant_freq"]/ length(unique(selection_output$MAF$Unique_patient_identifier)))
  #   }else{
  #     
  #     NRSI_vec <- selection_data[sub_index,"selection_intensity"] * (apply(contribution_matrix,2,sum)/length(unique(selection_output$MAF$Unique_patient_identifier)))
  #     
  #     
  #   }
  #   NRSI_matrix[sub_index,] <- as.numeric(NRSI_vec)
  #   
  #   # print(sub_index /nrow(selection_data))
  #   
  # }
  # 
  
  
  
  contribution_matrix_list <- sapply(contribution_mat_function_output, function(x) x$contribution_matrix)
  names(contribution_matrix_list) <- selection_data$variant
  
  NRSI_matrix <-  t(sapply(contribution_mat_function_output, function(x) x$NRSI_vec))
  
  rownames(NRSI_matrix) <- selection_data$variant
  colnames(NRSI_matrix) <- rownames(signatures)
  
  return(list(NRSI_matrix = NRSI_matrix,
              contribution_matrix= contribution_matrix_list))
}



