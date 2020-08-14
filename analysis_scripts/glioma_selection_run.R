#!/usr/bin/env Rscript

# set up data ----- 

inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
rownames(inputs) <- args[seq(1,length(args),2)]
inputs[,1] <- trimws(inputs[,1])


library(cancereffectsizeR)
analysis <- cancereffectsizeR::CESAnalysis(genome = "hg19")


# load in MAFs ---- 

if(!("-ready_MAF" %in% rownames(inputs))){
  MAF_for_analysis <- read.delim(file = inputs["-NCI_MAF",], header = T, stringsAsFactors = F,skip = inputs["-NCI_skip",])
  MAF_for_analysis$Tumor_Sample_Barcode = cancereffectsizeR::consolidate_tcga_tumors_by_patient(MAF_for_analysis$Tumor_Sample_Barcode)
  
  analysis = load_maf(analysis, maf = MAF_for_analysis, 
                      chain_file = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain")
  
  analysis@maf <- dplyr::distinct(analysis@maf) # take out duplications imposed by liftOver
  # MAF_for_analysis <- cancereffectsizeR::hg_converter(chain = "/ysm-gpfs/pi/townsend/general_genome_info/hg38ToHg19.over.chain",maf_to_convert = MAF_for_analysis)
  
  if("-Local_MAF" %in% rownames(inputs)){
    local_maf <- read.delim(file = inputs["-Local_MAF",],header = T,stringsAsFactors = F)
    
    local_maf <- dplyr::distinct(local_maf)
    if("Patient_ID" %in% colnames(local_maf)){
      colnames(local_maf)[which(colnames(local_maf) == "Patient_ID") ] <- "Tumor_Sample_Barcode"
    }
    
    
    if("Chrom" %in% colnames(local_maf)){
      colnames(local_maf)[which(colnames(local_maf) == "Chrom") ] <- "Chromosome"
    }
    
    # MAF_for_analysis <- cancereffectsizeR::merging_TCGA_and_local_MAFdata_function(NCI_data = MAF_for_analysis,Local_data = local_maf,check_for_same_tumor = T)
    analysis <- load_maf(analysis, maf = local_maf)
  }
  
}else{
  MAF_for_analysis <- get(load(inputs["-ready_MAF",]))
  
  MAF_for_analysis <- dplyr::distinct(MAF_for_analysis[,c("Chromosome","Start_Position","Tumor_Sample_Barcode","Tumor_Seq_Allele2","Tumor_Allele","Reference_Allele")])
  
  analysis <- load_maf(analysis, maf = MAF_for_analysis)
  
  
}


MAF_for_analysis <- analysis@maf

save(MAF_for_analysis, file = paste("output_data/",inputs["-tumor_name",],"_","MAF.RData",sep=""))


signatures_to_remove_lgg = cancereffectsizeR::suggest_cosmic_v3_signatures_to_remove(cancer_type = "CNS-LGG", treatment_naive = T)
signatures_to_remove_gbm = cancereffectsizeR::suggest_cosmic_v3_signatures_to_remove(cancer_type = "CNS-GBM", treatment_naive = T)

# extract that artifact signatures


load("output_data/tumor_type_designation.RData")

# just LGG
lgg_tumors <- tumor_type_designation$tumor_ID[tumor_type_designation$tumor_type=="LGG"]
lgg_maf <- MAF_for_analysis[MAF_for_analysis$Unique_Patient_Identifier %in% lgg_tumors,]


lgg_ces = load_maf(cesa = CESAnalysis("hg19"), maf = lgg_maf)
lgg_ces = trinuc_mutation_rates(lgg_ces,cores = as.numeric(inputs["-cores",]),
                                signatures_to_remove = signatures_to_remove_lgg,use_dS_exome2genome = T)

# just GBM
gbm_tumors <- tumor_type_designation$tumor_ID[tumor_type_designation$tumor_type=="GBM"]
gbm_maf <- MAF_for_analysis[MAF_for_analysis$Unique_Patient_Identifier %in% gbm_tumors,]

gbm_ces = load_maf(cesa = CESAnalysis("hg19"), maf = gbm_maf)
gbm_ces = trinuc_mutation_rates(gbm_ces,cores = as.numeric(inputs["-cores",]),
                                signatures_to_remove = signatures_to_remove_gbm,use_dS_exome2genome = T)



analysis <- CESAnalysis("hg19")
analysis <- load_maf(cesa = analysis, maf = lgg_maf)
analysis <- load_maf(cesa = analysis, maf=gbm_maf)

combined_rates <- rbind(lgg_ces$trinuc_rates, gbm_ces$trinuc_rates)
analysis = set_trinuc_rates(cesa = analysis, trinuc_rates = combined_rates)

analysis@trinucleotide_mutation_weights$signature_weights_table <- rbind(get_signature_weights(cesa = lgg_ces, include_tumors_without_data = T),get_signature_weights(cesa = gbm_ces, include_tumors_without_data = T))

# Calculate trinucleotide mutation weightings using deconstructSigs ----

# Calculate gene-level mutation rates using dNdScv

analysis <- cancereffectsizeR::gene_mutation_rates(analysis, 
                                                   covariate_file = if("-cov" %in% rownames(inputs)) inputs["-cov",] else NULL)

# Assign genes to MAF, keeping assignments consistent with dndscv when possible
analysis <- cancereffectsizeR::annotate_variants(analysis)



analysis <-  ces_snv(cesa = analysis, 
                     cores=as.numeric(inputs["-cores",]),
                     include_nonrecurrent_variants = F) 

save(analysis,file = paste("output_data/",inputs["-tumor_name",],"_","_selection_output.RData",sep=""))

results <- as.data.frame(analysis@selection_results[tumors_with_variant > 0])



save(results,file = paste("output_data/",inputs["-tumor_name",],"_","_selection_output_df.RData",sep=""))



source("analysis_scripts/population_scaled_selection_function_glioma.R")





# #
NRSI_output <- net_realized_selection_calculation(selection_output = analysis)


save(NRSI_output,file = paste("output_data/",inputs["-tumor_name",],"_","NRSI_output.RData",sep=""))


source("analysis_scripts/population_scaled_selection_per_tumor_glioma.R")

population_scaled_effect_per_tumor_data <- population_scaled_effect_per_tumor(selection_output = analysis,
                                                                              NRSI_output = NRSI_output,
                                                                              counter = F,
                                                                              remove_silent_away_from_splice = F)

save(population_scaled_effect_per_tumor_data,file = paste("output_data/",inputs["-tumor_name",],"_population_scaled_effect_per_tumor_data.RData",sep=""))

source("analysis_scripts/divergence_calculation_glioma.R")


tumor_JSD_df <- divergence_calculation(trinuc_weights = 
                                         rbind(cancereffectsizeR::get_signature_weights(cesa = gbm_ces,include_tumors_without_data = T),
                                               cancereffectsizeR::get_signature_weights(cesa = lgg_ces,include_tumors_without_data = T)),
                                       NRSI_per_tumor = population_scaled_effect_per_tumor_data$NRSI_prop)

tumor_JSD_df$tumor_type <- inputs["-tumor_name",]

save(tumor_JSD_df, file = paste("output_data/",inputs["-tumor_name",],"_","tumor_JSD_df.RData",sep=""))


signature_weights_all_tumors <- rbind(get_signature_weights(cesa = lgg_ces, include_tumors_without_data = T),
                                      get_signature_weights(cesa = gbm_ces, include_tumors_without_data = T))
save(signature_weights_all_tumors, file = paste("output_data/",inputs["-tumor_name",],"_","signature_weights_all_tumors.RData",sep=""))

