
# run all effect size calculations


# IDH_mutant_M_glioma ----- 

args <- c("-ready_MAF", 
          "output_data/IDH_mutant_M_glioma.RData", 
          "-tumor_name", 
          "IDH_mutant_M",
          "-cores","10","-cov_name", "lgg_pca")

source("analysis_scripts/glioma_selection_run.R")

rm(list=ls())

# IDH_mutant_F_glioma ------ 

args <- c("-ready_MAF", 
          "output_data/IDH_mutant_F_glioma.RData", 
          "-tumor_name", 
          "IDH_mutant_F",
          "-cores","10","-cov_name", "lgg_pca")

source("analysis_scripts/glioma_selection_run.R")

rm(list=ls())

# IDH_WT_F_glioma -----

args <- c("-ready_MAF", 
          "output_data/IDH_WT_F_glioma.RData", 
          "-tumor_name", 
          "IDH_WT_F",
          "-cores","10","-cov_name", "lgg_pca")

source("analysis_scripts/glioma_selection_run.R")

rm(list=ls())

# IDH_WT_M_glioma ---- 


args <- c("-ready_MAF", 
          "output_data/IDH_WT_M_glioma.RData", 
          "-tumor_name", 
          "IDH_WT_M",
          "-cores","10","-cov_name", "lgg_pca")

source("analysis_scripts/glioma_selection_run.R")

rm(list=ls())
