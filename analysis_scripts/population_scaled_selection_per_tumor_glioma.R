# load("../../cancereffectsizeR/cluster_work/results/trinuc_all/LUAD_FALSE_NRSI_output_ML_cores_trinuc.RData")
# load("../../cancereffectsizeR/cluster_work/results/trinuc_all/LUAD_FALSE_selection_output_ML_cores_trinuc.RData")

# load("../cluster_work/results/trinuc_all/BLCA_FALSE_NRSI_output_ML_cores_trinuc.RData")
# load("../cluster_work/results/trinuc_all/BLCA_FALSE_selection_output_ML_cores_trinuc.RData")
# names(NRSI_output$contribution_matrix)
# head(NRSI_output$NRSI_matrix)
# NRSI_output$contribution_matrix[["KRAS G12V"]]

library(tidyverse)
library(cowplot)


# head(selection_output@maf)


# 
# selection_data <- selection_output
# NRSI_data <- NRSI_output

population_scaled_effect_per_tumor <- function(selection_output, 
                                               NRSI_output,
                                               counter=F,
                                               orig_30_sigs=F,
                                               cores_for_contribution_mat=NULL,
                                               remove_silent_away_from_splice=F){
  
  # selection_data <- cancereffectsizeR::selection_results_converter(selection_output,min_recurrence = 2)
  
  selection_data <- as.data.frame(selection_output@selection_results[tumors_with_variant > 1])
  
  
  selection_output_df <- as.data.frame(selection_output@maf)
  
  annot_number <- sapply(selection_output_df$assoc_aa_mut,length)
  annot_index <- which(annot_number > 1)
  if(length(annot_index)>0){
    for(annot in 1:length(annot_index)){

        these_variants <- unlist(selection_output_df$assoc_aa_mut[annot_index[annot]])
        if(length(which(selection_data$variant %in% these_variants)>0)){
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
        
      }
      
    }
  }
  # 
  
  # adding variant matcher
  
  selection_output_df$variant_matcher <- selection_output_df$assoc_aa_mut
  selection_output_df$variant_matcher[is.na(selection_output_df$variant_matcher)] <- 
    selection_output_df$snv_id[is.na(selection_output_df$variant_matcher)]
  
  
  # removing the silent variants that are not by a splice site
  # if(remove_silent_away_from_splice){
  #   
  #   # getting AA change from variant
  #   selection_data$variant_minus_gene <- NA
  #   for(row_ind in 1:nrow(selection_data)){
  #     selection_data$variant_minus_gene[row_ind] <- 
  #       trimws(gsub(pattern = selection_data$gene[row_ind],replacement = "",x = selection_data$variant[row_ind]))
  #   }
  #   
  #   selection_data$to_remove <- F
  #   
  #   # telling if the variant is silent AA 
  #   for(row_ind in 1:nrow(selection_data)){
  #     if(length(strsplit(x = selection_data$variant_minus_gene[row_ind]," ")[[1]])==1){
  #       var_1 <- substr(x = selection_data$variant_minus_gene[row_ind],1,1)
  #       var_2 <- substr(x = selection_data$variant_minus_gene[row_ind],
  #                       nchar(selection_data$variant_minus_gene[row_ind]),
  #                       nchar(selection_data$variant_minus_gene[row_ind]))
  #       
  #       if(var_1 == var_2){
  #         
  #         if(!any(selection_output@maf$next_to_splice[
  #           which(selection_output@maf$Gene_name==selection_data$gene[row_ind] &
  #                 selection_output@maf$unique_variant_ID_AA == selection_data$variant_minus_gene[row_ind])
  #           ])){
  #          
  #           selection_data$to_remove[row_ind] <- T
  #            
  #         }
  #         
  #         
  #       }
  #     }
  #     
  #      
  #   }
  #   
  #   
  #   # remove them
  #   
  #   selection_data <- selection_data[-which(selection_data$to_remove==T),]
  #   
  #   
  # }
  
  
  rownames(selection_data) <- selection_data$variant
  NRSI_data <- NRSI_output
  
  
  tumors_without_recurrent_variants <- NULL
  
  # NRSI_tumor_list <- vector(mode = "list",length = length(unique(selection_output@maf$U)))
  # names(NRSI_tumor_list) <- unique(selection_output@maf$Unique_Patient_Identifier)
  
  
  
  # mutational_flux_tumor_list <- vector(mode = "list",length = length(unique(selection_output@maf$Unique_Patient_Identifier)))
  # names(mutational_flux_tumor_list) <- unique(selection_output@maf$Unique_Patient_Identifier)
  
  
  
  recurrent_variants <- selection_data$variant
  
  # MAF_from_analysis <- selection_output@maf
  # 
  # MAF_from_analysis %>%
  #   filter(Unique_Patient_Identifier == "TCGA-73-4670" & Start_Position == 19378375)
  
  MAF_from_analysis <- selection_output_df
  
  # MAF_from_analysis %>%
  #   filter(Unique_Patient_Identifier == "TCGA-73-4670" & Start_Position == 19378375)
  
  
  # MAF_from_analysis$variant_ID_to_match <- MAF_from_analysis$unique_variant_ID_AA
  # MAF_from_analysis[which(sapply(strsplit(MAF_from_analysis$unique_variant_ID_AA,split = " "), function(x) length(x))==1),"variant_ID_to_match"] <- paste(MAF_from_analysis$Gene_name[which(sapply(strsplit(MAF_from_analysis$unique_variant_ID_AA,split = " "), function(x) length(x))==1)],MAF_from_analysis$unique_variant_ID_AA[which(sapply(strsplit(MAF_from_analysis$unique_variant_ID_AA,split = " "), function(x) length(x))==1)],sep = " ")
  # 
  
  
  
  NRSI_tumor_list_function <- function(tumor_index,counter){
    
    if(length(
      which(
        MAF_from_analysis$Unique_Patient_Identifier==unique(MAF_from_analysis$Unique_Patient_Identifier)[tumor_index] & 
        MAF_from_analysis$variant_matcher %in% recurrent_variants)
    ) > 0){
      
      
      NRSI_tumor_list <- cbind(
        unique(MAF_from_analysis$Unique_Patient_Identifier)[tumor_index],
        expand.grid(
          MAF_from_analysis$variant_matcher[
            which(MAF_from_analysis$Unique_Patient_Identifier == unique(MAF_from_analysis$Unique_Patient_Identifier)[tumor_index] & 
                    MAF_from_analysis$variant_matcher %in% recurrent_variants)],
          colnames(NRSI_data$NRSI_matrix),stringsAsFactors = F))
      NRSI_tumor_list$NRSI <- NA
      colnames(NRSI_tumor_list) <- c("Tumor","Variant","Signature","NRSI")
      
      
      mutational_flux_tumor_list <- NRSI_tumor_list
      colnames(mutational_flux_tumor_list)[4] <- "Mutational_flux"
      
      
      unique_variant_length <- length(unique(NRSI_tumor_list$Variant))
      unique_variants <- unique(NRSI_tumor_list$Variant)
      this_tumor <- unique(MAF_from_analysis$Unique_Patient_Identifier)[tumor_index]
      for(unique_variant in 1:unique_variant_length){
        
        NRSI_tumor_list$NRSI[NRSI_tumor_list$Variant==unlist(unique_variants)[unique_variant]] <- 
          NRSI_data$contribution_matrix[[unlist(unique_variants[unique_variant])]][this_tumor,] * 
          selection_data[unlist(unique_variants)[unique_variant],"selection_intensity"]
        
        
        mutational_flux_tumor_list$Mutational_flux[NRSI_tumor_list$Variant==unlist(unique_variants)[unique_variant]] <- NRSI_data$contribution_matrix[[unlist(unique_variants)[unique_variant]]][
          this_tumor,]
        
      }
      
      # for(variant_index in 1:nrow(NRSI_tumor_list)){
      #   
      #   NRSI_tumor_list$NRSI[variant_index] <-
      #     selection_data[NRSI_tumor_list$Variant[variant_index],"selection_intensity"] *
      #     NRSI_data$contribution_matrix[[NRSI_tumor_list$Variant[variant_index]]][as.character(NRSI_tumor_list[variant_index,"Tumor"]),NRSI_tumor_list$Signature[variant_index]]
      #   
      #   
      #   mutational_flux_tumor_list$Mutational_flux[variant_index] <-
      #     NRSI_data$contribution_matrix[[NRSI_tumor_list$Variant[variant_index]]][as.character(NRSI_tumor_list[variant_index,"Tumor"]),NRSI_tumor_list$Signature[variant_index]]
      #   
      #   
      # }
      
      
      return(list(NRSI_tumor_list=NRSI_tumor_list,mutational_flux_tumor_list=mutational_flux_tumor_list))
    }else{
      # message(paste("No recurrent variants within ", names(NRSI_tumor_list)[tumor_index]))
      
      return(list(tumors_without_recurrent_variants=unique(MAF_from_analysis$Unique_Patient_Identifier)[tumor_index]))
      
      # tumors_without_recurrent_variants <- c(tumors_without_recurrent_variants,names(NRSI_tumor_list)[tumor_index])
    }
    if(counter){
      print(tumor_index/length(unique(selection_output@maf$Unique_Patient_Identifier)))
    }
    
    # return()
  }
  
  if(is.null(cores_for_contribution_mat)){
    cores_for_contribution_mat <- parallel::detectCores()
  }
  
  NRSI_per_tumor_list <- parallel::mclapply(X = 1:length(unique(selection_output@maf$Unique_Patient_Identifier)),
                                            FUN = NRSI_tumor_list_function,
                                            mc.cores = cores_for_contribution_mat,
                                            counter=counter)
  
  # NRSI_per_tumor_list <- pbapply::pblapply(X = 1:length(unique(selection_output@maf$Unique_Patient_Identifier)),
  #                                           FUN = NRSI_tumor_list_function,
  #                                           cl = cores_for_contribution_mat,
  #                                           counter=F)
  
  # 
  # 
  # for(tumor_index in 1:length(NRSI_tumor_list)){
  #   # build data structure
  #   
  #   if(length(
  #     which(
  #       MAF_from_analysis$Unique_Patient_Identifier==names(NRSI_tumor_list)[tumor_index] & 
  #       MAF_from_analysis$variant_ID_to_match %in% recurrent_variants)
  #   ) > 0){
  #     NRSI_tumor_list[[tumor_index]] <- cbind(
  #       names(NRSI_tumor_list)[tumor_index],
  #       expand.grid(
  #         MAF_from_analysis$variant_ID_to_match[
  #           which(MAF_from_analysis$Unique_Patient_Identifier == names(NRSI_tumor_list)[tumor_index] & 
  #                   MAF_from_analysis$variant_ID_to_match %in% recurrent_variants)],
  #         colnames(NRSI_data$NRSI_matrix),stringsAsFactors = F))
  #     NRSI_tumor_list[[tumor_index]]$NRSI <- NA
  #     colnames(NRSI_tumor_list[[tumor_index]]) <- c("Tumor","Variant","Signature","NRSI")
  #     
  #     
  #     mutational_flux_tumor_list[[tumor_index]] <- NRSI_tumor_list[[tumor_index]]
  #     colnames(mutational_flux_tumor_list[[tumor_index]])[4] <- "Mutational_flux"
  #     
  #     for(variant_index in 1:nrow(NRSI_tumor_list[[tumor_index]])){
  #       
  #       NRSI_tumor_list[[tumor_index]]$NRSI[variant_index] <-
  #         selection_data[NRSI_tumor_list[[tumor_index]]$Variant[variant_index],"selection_intensity"] *
  #         NRSI_data$contribution_matrix[[NRSI_tumor_list[[tumor_index]]$Variant[variant_index]]][names(NRSI_tumor_list)[tumor_index],NRSI_tumor_list[[tumor_index]]$Signature[variant_index]]
  #       
  #       
  #       mutational_flux_tumor_list[[tumor_index]]$Mutational_flux[variant_index] <-
  #         NRSI_data$contribution_matrix[[NRSI_tumor_list[[tumor_index]]$Variant[variant_index]]][names(NRSI_tumor_list)[tumor_index],NRSI_tumor_list[[tumor_index]]$Signature[variant_index]]
  #       
  #       
  #     }
  #     
  #   }else{
  #     # message(paste("No recurrent variants within ", names(NRSI_tumor_list)[tumor_index]))
  #     tumors_without_recurrent_variants <- c(tumors_without_recurrent_variants,names(NRSI_tumor_list)[tumor_index])
  #   }
  #   if(counter){
  #     print(tumor_index/length(NRSI_tumor_list))
  #   }
  # }
  
  tumors_without_recurrent_variants <- unlist(sapply(NRSI_per_tumor_list, function(x) x$tumors_without_recurrent_variants))
  
  NRSI_tumor_list <- sapply(NRSI_per_tumor_list, function(x) x$NRSI_tumor_list) 
  NRSI_per_tumor_expanded <- do.call(what = rbind,args = NRSI_tumor_list)
  
  mutational_flux_tumor_list <- sapply(NRSI_per_tumor_list, function(x) x$mutational_flux_tumor_list) 
  mutational_flux_tumor_expanded <- do.call(what = rbind,args = mutational_flux_tumor_list)
  
  
  message(paste("Percent of tumors without any recurrent variants:",round(length(tumors_without_recurrent_variants)/length(unique(MAF_from_analysis$Unique_Patient_Identifier)),2)))
  
  # NRSI_per_tumor_expanded <- do.call(what = rbind,args = NRSI_tumor_list)
  
  # mutational_flux_tumor_expanded <- do.call(what = rbind, args = mutational_flux_tumor_list)
  
  if(orig_30_sigs){
    signature_mat <- matrix(data=c(paste("Signature.",1:30,sep=""),c("Spontaneous deamination, aging (1)",
                                                                     "APOBEC (2)",
                                                                     "BRCA1 and BRCA2 (3)",
                                                                     "Tobacco (4)",
                                                                     "Unknown (5)",
                                                                     "Defective DNA mismatch repair (6)",
                                                                     "Ultraviolet light exposure (7)",
                                                                     "Unknown (8)",
                                                                     "Polymerase eta (9)",
                                                                     "Polymerase POLE (10)",
                                                                     "Alkylating agents (11)",
                                                                     "Unknown (12)",
                                                                     "APOBEC (13)",
                                                                     "Unknown (14)",
                                                                     "Defective DNA mismatch repair (15)",
                                                                     "Unknown (16)",
                                                                     "Unknown (17)",
                                                                     "Unknown (18)",
                                                                     "Unknown (19)",
                                                                     "Defective DNA mismatch repair (20)",
                                                                     "Unknown (21)",
                                                                     "Aristolochic acid (22)",
                                                                     "Unknown (23)",
                                                                     "Aflatoxin (24)",
                                                                     "Unknown (25)",
                                                                     "Defective DNA mismatch repair (26)",
                                                                     "Unknown (27)",
                                                                     "Unknown (28)",
                                                                     "Tobacco chewing (29)",
                                                                     "Unknown (30)")),nrow=30)
    rownames(signature_mat) <- signature_mat[,1]
    
    
    # NRSI_per_tumor_expanded <- do.call(what = rbind,args = NRSI_per_tumor$)
    # head(NRSI_per_tumor$NRSI_per_tumor_expanded)
    
    # NRSI_per_tumor_expanded <- NRSI_per_tumor$NRSI_per_tumor_expanded
    
    NRSI_per_tumor_expanded$Signature <- as.character(signature_mat[NRSI_per_tumor_expanded$Signature,2])
    
    mutational_flux_tumor_expanded$Signature <- as.character(signature_mat[mutational_flux_tumor_expanded$Signature,2])
  }
  
  
  NRSI_plot_all_prop <- NRSI_per_tumor_expanded %>%
    group_by(Tumor) %>%
    mutate(pct = NRSI / sum(NRSI))
  
  mut_flux_plot_all_prop <- mutational_flux_tumor_expanded %>%
    group_by(Tumor) %>%
    mutate(pct = Mutational_flux / sum(Mutational_flux))
  
  NRSI_plot_all_prob_max <- NRSI_plot_all_prop %>%
    group_by(Signature) %>%
    summarize(total_signature = sum(pct)) %>%
    arrange(desc(total_signature))
  
  NRSI_descending <- NRSI_plot_all_prop %>%
    filter(Signature == NRSI_plot_all_prob_max$Signature[1]) %>%
    group_by(Tumor) %>%
    summarize(total_sig = sum(pct)) %>%
    arrange(desc(total_sig))
  
  
  
  
  NRSI_plot_all_prop$Tumor <- factor(x = NRSI_plot_all_prop$Tumor,levels = NRSI_descending$Tumor)
  
  NRSI_per_tumor_expanded$Tumor <- factor(x = NRSI_per_tumor_expanded$Tumor,levels = as.character(NRSI_descending$Tumor))
  
  NRSI_plot <- ggplot(data = NRSI_per_tumor_expanded) +
    geom_bar(aes(x=Tumor,y=NRSI,fill=Signature),color="white",stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
    labs(title="Tumor-specific NRSSC",x="Tumor",y="NRSCC")
  
  NRSI_plot_pct <- ggplot(data = NRSI_plot_all_prop) +
    geom_bar(aes(x=Tumor,y=pct,fill=Signature),stat="identity", color="white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
    labs(title="Tumor-specific NRSSC, proportion of total",x="Tumor",y="NRSCC")
  
  
  subs_per_tumor <- NRSI_per_tumor_expanded %>%
    group_by(Tumor) %>%
    summarize(variant_count = length(unique(Variant)))
  
  subs_per_tumor$Tumor <- factor(subs_per_tumor$Tumor, levels = NRSI_descending$Tumor)
  
  
  NRSI_subs_per_tumor <- ggplot(data = subs_per_tumor) +
    geom_bar(aes(x=Tumor,y=variant_count),stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
    labs(title="Recurrent substitutions within each tumor",x="Tumor",y="Substitutions per tumor")
  
  
  plot_combined <- cowplot::plot_grid(NRSI_plot_pct + theme(legend.position = "none"),
                                      NRSI_plot + theme(legend.position = "none"),
                                      NRSI_subs_per_tumor + theme(legend.position = "none"),
                                      align = 'hv',ncol = 1,labels = "AUTO")
  
  legend <- cowplot::get_legend(plot = NRSI_plot_pct)
  plot_combined_legend <- cowplot::plot_grid(plot_combined, legend, rel_widths = c(20,1))
  
  return(list(NRSI_tumor_list=NRSI_tumor_list,
              NRSI_per_tumor_expanded=NRSI_per_tumor_expanded,
              NRSI_plot=NRSI_plot,
              NRSI_plot_pct=NRSI_plot_pct,
              NRSI_subs_per_tumor=NRSI_subs_per_tumor,
              plot_combined_legend=plot_combined_legend,
              NRSI_prop = NRSI_plot_all_prop,
              mutational_flux_tumor_list=mutational_flux_tumor_list,
              mut_flux_plot_all_prop=mut_flux_plot_all_prop,
              tumors_without_recurrent_variants=tumors_without_recurrent_variants))
  
}


# NRSI_per_tumor <- NRSI_per_tumor(selection_output = selection_output,NRSI_output  = NRSI_output)
# # 
# 
# 
# signature_mat <- matrix(data=c(paste("Signature.",1:30,sep=""),c("Spontaneous deamination, aging (1)",
#                                                                  "APOBEC (2)",
#                                                                  "BRCA1 and BRCA2 (3)",
#                                                                  "Tobacco (4)",
#                                                                  "Unknown (5)",
#                                                                  "Defective DNA mismatch repair (6)",
#                                                                  "Ultraviolet light exposure (7)",
#                                                                  "Unknown (8)",
#                                                                  "Polymerase eta (9)",
#                                                                  "Polymerase POLE (10)",
#                                                                  "Alkylating agents (11)",
#                                                                  "Unknown (12)",
#                                                                  "APOBEC (13)",
#                                                                  "Unknown (14)",
#                                                                  "Defective DNA mismatch repair (15)",
#                                                                  "Unknown (16)",
#                                                                  "Unknown (17)",
#                                                                  "Unknown (18)",
#                                                                  "Unknown (19)",
#                                                                  "Defective DNA mismatch repair (20)",
#                                                                  "Unknown (21)",
#                                                                  "Aristolochic acid (22)",
#                                                                  "Unknown (23)",
#                                                                  "Aflatoxin (24)",
#                                                                  "Unknown (25)",
#                                                                  "Defective DNA mismatch repair (26)",
#                                                                  "Unknown (27)",
#                                                                  "Unknown (28)",
#                                                                  "Tobacco chewing (29)",
#                                                                  "Unknown (30)")),nrow=30)
# rownames(signature_mat) <- signature_mat[,1]
# 
# 
# # NRSI_per_tumor_expanded <- do.call(what = rbind,args = NRSI_per_tumor$)
# # head(NRSI_per_tumor$NRSI_per_tumor_expanded)
# 
# NRSI_per_tumor_expanded <- NRSI_per_tumor$NRSI_per_tumor_expanded
# 
# NRSI_per_tumor_expanded$Signature <- as.character(signature_mat[NRSI_per_tumor_expanded$Signature,2])
# 
# NRSI_plot_all_prop <- NRSI_per_tumor_expanded %>%
#   group_by(Tumor) %>%
#   mutate(pct = NRSI / sum(NRSI))
# 
# NRSI_descending_sig4 <- NRSI_plot_all_prop %>%
#   filter(Signature == "Tobacco (4)") %>%
#   group_by(Tumor) %>%
#   summarize(total_sig4 = sum(pct)) %>%
#   arrange(desc(total_sig4))
# 
# NRSI_plot_all_prop$Tumor <- factor(x = NRSI_plot_all_prop$Tumor,levels = NRSI_descending_sig4$Tumor)
# 
# NRSI_per_tumor_expanded$Tumor <- factor(x = NRSI_per_tumor_expanded$Tumor,levels = as.character(NRSI_descending_sig4$Tumor))
# 
# 
# 
# # library(tidyverse)
# 
# 
# NRSI_per_tumor_expanded_nozero <- NRSI_per_tumor_expanded[-which(NRSI_per_tumor_expanded$NRSI==0),]
# 
# NRSI_plot <- ggplot(data = NRSI_per_tumor_expanded) +
#   geom_bar(aes(x=Tumor,y=NRSI,fill=Signature),color="white",stat="identity") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
#   labs(title="LUAD tumor-specific NRSSC",x="Tumor",y="NRSCC")
# 
# ggsave(filename = "../local_work/NRSI_plot.pdf",width = 80,limitsize = F)
# 
# NRSI_plot_pct <- ggplot(data = NRSI_plot_all_prop) +
#   geom_bar(aes(x=Tumor,y=pct,fill=Signature),stat="identity", color="white") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
#   labs(title="LUAD tumor-specific NRSSC",x="Tumor",y="NRSCC")
# 
# ggsave(filename = "../local_work/NRSI_plot_pct.pdf",width = 80,limitsize = F,plot = NRSI_plot_pct)
# 
# subs_per_tumor <- NRSI_per_tumor$NRSI_per_tumor_expanded %>%
#   group_by(Tumor) %>%
#   summarize(variant_count = length(unique(Variant)))
# 
# subs_per_tumor$Tumor <- factor(subs_per_tumor$Tumor, levels = NRSI_descending_sig4$Tumor)
# 
# NRSI_subs_per_tumor <- ggplot(data = subs_per_tumor) +
#   geom_bar(aes(x=Tumor,y=variant_count),stat="identity") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
#   labs(title="Substitutions within LUAD tumors",x="Tumor",y="Substitutions per tumor")
# 
# ggsave(filename = "../local_work/NRSI_subs_per_tumor.pdf",width = 80,limitsize = F,plot = NRSI_subs_per_tumor)
# 
# library(cowplot)
# 
# plot_combined <- plot_grid(NRSI_plot_pct + theme(legend.position = "none"),
#                            NRSI_plot + theme(legend.position = "none"),
#                            NRSI_subs_per_tumor + theme(legend.position = "none"),
#                            align = 'hv',ncol = 1,labels = "AUTO")
# 
# legend <- get_legend(plot = NRSI_plot_pct)
# plot_combined_legend <- plot_grid(plot_combined, legend, rel_widths = c(20,1))
# save_plot(filename = "../local_work/NRSI_plot_combined.pdf",plot = plot_combined_legend,base_height = 3.35*3,base_width = 90,limitsize=F)
# 
