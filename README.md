

These scripts reproduce the analyses and figures implemented within the manuscript `Environmental and sex-specific molecular signatures of glioma causation` by Elizabeth B. Claus MD PhD\* and Vincent Cannataro PhD\*, Stephen G. Gaffney PhD, and Jeffrey P. Townsend PhD

\*: co-first-authors

Steps in the analysis:

1. In the `input_data` directory, you need to download the freely available raw data from GLASS and TCGA. These files are called within `data_preprocess.Rmd`, so check the comments within that script for more information

    1.1. From GLASS https://www.synapse.org/#!Synapse:syn17038081/tables/ : `variants_anno.csv`, `variants_passgeno.csv`, `clinical_cases.csv`, `clinical_surgeries.csv` 
    
    1.2. From TCGA: `TCGA.LGG.mutect.1e0694ca-fcde-41d3-9ae3-47cfaf527f25.DR-10.0.somatic.maf` and `TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf`. These are the two open MAF files annotated via mutect for the LGG and GBM TCGA projects, respectively.
    
2. Run `data_preprocess.Rmd` to gather all the glioma data, and split the data into the various sex and IDH status designations. 
 
3. Run `analysis_scrips/glioma_selection_run_wrap.R`

4. Run `analysis_scripts/effect_size_results_combiner.R`

5. Run `analysis_scripts/per_tumor_scaled_effect_results_combiner.R`

6. Run 


