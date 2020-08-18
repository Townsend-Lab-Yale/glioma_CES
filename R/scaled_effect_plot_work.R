# load("../cluster_results/NRSI_from_cluster/LUAD_FALSE_NRSI_output_ML_cores_trinuc.RData")

library(tidyverse)


# this_tumor_type <- "IDH_mutant_F"
# 
# NRSI_input <- per_tumor_data_main[per_tumor_data_main$tumor_type==this_tumor_type,]
# top_variant_number=15
# make_plot=T
# tumor_type=""
# sqrt_scale=F
# plot_text_size=8
# signature_version="v3"
# legend_rows=4
# top_sigs_to_plot=5
# sig_text_vjust=-2
# sig_text_hjust=0
# pct_text_vjust=2
# pct_text_hjust=-.5
# sig_text_angle=90
# pct_text_angle=90
# given_new_df=NULL
# y_boost=10
# 
# tumortype_ind <- 1
# this_tumor_type <- unique(per_tumor_data_main$tumor_type)[tumortype_ind]
# 
# NRSI_input = per_tumor_data_main[per_tumor_data_main$tumor_type==this_tumor_type,]
# top_variant_number = 10
# tumor_type = this_tumor_type
# plot_text_size = 11
# legend_rows = 6
# top_sigs_to_plot = 67
# sig_text_vjust = 0.5
# sig_text_hjust = 0
# pct_text_vjust = 0
# pct_text_hjust = 0
# pct_text_angle = 0
# sig_text_angle = 90
# # given_new_df = test_plot$new_NRSI_df,
# y_boost = 15
# 
# NRSI_input$Variant <- NRSI_input$Name

population_scaled_selection_melter_and_plotter <- function(NRSI_input, 
                                                           top_variant_number=10, 
                                                           # make_plot=T,
                                                           tumor_name_for_title="",
                                                           # sqrt_scale=F,
                                                           plot_text_size=8
                                                           # legend_rows=4,
                                                           # top_sigs_to_plot=5,
                                                           # sig_text_vjust=-2,
                                                           # sig_text_hjust=0,
                                                           # pct_text_vjust=2,
                                                           # pct_text_hjust=-.5,
                                                           # sig_text_angle=90,
                                                           # pct_text_angle=90,
                                                           # given_new_df=NULL,
                                                           # y_boost=10
                                                           ){
  
  # NRSI_input <- NRSI_output
  # top_variant_number <- 15
  # make_plot <- T 
  # 
  
  
  NRSI_input <- NRSI_input %>%
    mutate(Variant = Name)
  # top_variant_number <- 10
  # plot_text_size <- 8
  # tumor_name_for_title = this_tumor_type
  
  
  data(signatures_names_matrix,package = "cancereffectsizeR")
  signature_mat <- signatures_names_matrix
  
  
  
  NRSI_input$Signature <- signature_mat[NRSI_input$Signature,2]
  
  NRSI_melted <- NRSI_input %>%
    group_by(Variant,Signature) %>%
    summarize(Population_scaled_selection=sum(pct)) %>%
    arrange(desc(Population_scaled_selection))
  
  
  
  # NRSI_melted <- NRSI_melted %>% 
  # arrange(desc(pct))
  
  top_variants <- as.character(unique(NRSI_melted$Variant)[1:top_variant_number])
  
  
  
  sum_of_signatures <- NRSI_melted %>%
    group_by(Signature) %>%
    summarize(signature_sum=sum(Population_scaled_selection,na.rm = T)) %>%
    mutate(signature_percent=round((signature_sum/sum(signature_sum,na.rm = T))*100,1)) %>%
    arrange(desc(signature_sum))
  
  new_NRSI_df <- expand.grid(Variant=c(top_variants,"Other variants"),Signature=as.character(sum_of_signatures$Signature))
  
  new_NRSI_df$Population_scaled_selection <- NA
  
  
  # populate this data frame
  
  
  
  for(i in 1:nrow(new_NRSI_df)){
    if(as.character(new_NRSI_df[i,"Variant"]) != "Other variants"){
      if(length(which(NRSI_melted[,"Variant"] == as.character(new_NRSI_df[i,"Variant"]) &
                      NRSI_melted[,"Signature"] == as.character(new_NRSI_df[i,"Signature"]))) > 0){
        
        new_NRSI_df[i,"Population_scaled_selection"] <- NRSI_melted[which(NRSI_melted[,"Variant"] == as.character(new_NRSI_df[i,"Variant"]) &
                                                                            NRSI_melted[,"Signature"] == as.character(new_NRSI_df[i,"Signature"])),"Population_scaled_selection"]
        
      }
      
    }else{
      
      new_NRSI_df[i,"Population_scaled_selection"] <- sum_of_signatures[which(sum_of_signatures$Signature==new_NRSI_df[i,"Signature"]),"signature_sum"] - sum(new_NRSI_df[as.character(new_NRSI_df$Signature)==as.character(new_NRSI_df[i,"Signature"]),"Population_scaled_selection"],na.rm = T)
      
    }
    
  }
  
  
  
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  
  
  new_NRSI_df$Signature <- factor( new_NRSI_df$Signature, levels= unique(new_NRSI_df$Signature))
  
  sum_of_signatures <- new_NRSI_df %>%
    group_by(Signature) %>%
    summarize(signature_sum=sum(Population_scaled_selection,na.rm = T)) %>%
    mutate(signature_percent=round((signature_sum/sum(signature_sum,na.rm = T))*100,1)) %>%
    arrange(desc(signature_sum)) %>%
    mutate(index = 1:nrow(.))
  
  index_mat <- sum_of_signatures$index
  names(index_mat) <- sum_of_signatures$Signature
  
  sum_of_signatures_sqrt <- new_NRSI_df %>%
    mutate(sqrt_nrsi = sqrt(Population_scaled_selection)) %>%
    group_by(Signature) %>%
    summarize(signature_sum=sum(sqrt_nrsi,na.rm = T)) %>%
    mutate(signature_percent=round((signature_sum/sum(signature_sum,na.rm = T))*100,1)) %>%
    arrange(desc(signature_sum))
  
  sum_of_signatures_sqrt$index <- index_mat[as.character(sum_of_signatures_sqrt$Signature)]
  
  sum_of_signatures_sqrt <- sum_of_signatures_sqrt %>%
    arrange((index))
  
  sum_of_signatures$new_y <- sum_of_signatures_sqrt$signature_sum
  
  library(scales)
  # This code was found at 
  # https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4 - Thanks Brian Diggs! 
  # And discussed and edited here: https://stackoverflow.com/a/24241954/8376488 - Thanks Jack Aidley! 
  
  
  
  
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    l <- gsub("0e\\+00","0",l)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
  }
  
  
  mysqrt_trans <- function() {
    domain <- c(0, Inf)
    transform <- function(x) x^(1/1.5)
    range <- transform(domain)
    trans_new("mysqrt",
              transform = transform,
              inverse = function(x) squish(x, range=range)^1.5,
              domain = domain)
  }
  
  
  new_NRSI_df <- new_NRSI_df %>%
    mutate(Population_scaled_selection_pct = Population_scaled_selection/sum(Population_scaled_selection))
  
  
  
  new_NRSI_df$Variant <- factor(new_NRSI_df$Variant,levels = unique(new_NRSI_df$Variant))
  
  
  
  geom_text_size <- (5/14)*plot_text_size
  # if(sqrt_scale){
  #   NRSI_bars <- ggplot(data = new_NRSI_df) +
  #     geom_bar(aes(x=Signature,y=sqrt(Population_scaled_selection),fill=Variant),stat="identity",color="black") +
  #     # geom_text(data=sum_of_signatures,aes(x=Signature,y=new_y,label=signature_percent,vjust=-.25)) +
  #     theme_classic() +
  #     scale_fill_manual(values = c(gg_color_hue(n = top_variant_number),"black"),name="Variant") +
  #     # scale_y_log10()+
  #     # scale_y_continuous(trans="mysqrt") +
  #     theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust=1,size=plot_text_size),
  #           axis.text.y = element_text(size=plot_text_size),
  #           axis.title  = element_text(size=plot_text_size))  +
  #     labs(title=tumor_type,x="Mutational process (signature)",y="Population scaled cancer effect size\n(square root scale)")
  # }else{
  #   
  # new_NRSI_plot_sub <- new_NRSI_df[new_NRSI_df$Signature %in% names(index_mat)[1:top_sigs_to_plot],]
  
  new_NRSI_plot_sub <- new_NRSI_df %>% 
    filter(Population_scaled_selection_pct>0)
  
  sum_of_signatures_toplot <- sum_of_signatures %>%
    filter(signature_sum>0)
  
  # new_NRSI_plot_sub$Variant <- forcats::fct_relabel(new_NRSI_plot_sub$Variant,
  #                                                   gsub(pattern = "_",replacement = " ")
  
  new_NRSI_plot_sub <- new_NRSI_plot_sub %>%
    filter(Population_scaled_selection_pct>0) %>%
    mutate(Signature = forcats::fct_rev(forcats::fct_drop(Signature))) 
  
  
  NRSI_bars <- ggplot(data = new_NRSI_plot_sub)+
    geom_bar(aes(x=Signature,
                 y=Population_scaled_selection_pct*100,
                 fill=Variant),
             stat="identity",
             color="white") +
    coord_flip() + 
    geom_text(data=sum_of_signatures_toplot,
              aes(x=Signature,y=0,label=Signature),
              size=geom_text_size,
              # angle=sig_text_angle,
              # vjust=1,
              hjust=0,
              fontface="bold")  +
    theme_classic() +
    scale_fill_manual(values = c(gg_color_hue(n = top_variant_number),"gray80"),name="Variant") +
    theme(axis.text.x = element_text(angle=0,vjust = 1,hjust=0.5),
          axis.text.y = element_blank(),
          text = element_text(size=plot_text_size),
          legend.text = element_text(size=plot_text_size*(8/10)),
          legend.title = element_blank())  +
    labs(x="Mutational process (signature)",y="Population scaled cancer effect size, percent of total",title=tumor_name_for_title) + 
    scale_y_continuous(expand = c(0,0)) + 
    # coord_cartesian(ylim=c(0,max(sum_of_signatures$signature_percent)+y_boost)) + 
    theme(legend.justification = c(1, 0), legend.position = c(1, .01)) + 
    guides(fill = guide_legend(ncol = 1)) + 
    theme(
      # axis.text.x=element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5)
    ) + 
    theme(legend.key.size = unit(0.3, "cm"),
          legend.box.margin=margin(0,0,0,0),
          legend.margin = margin(0,0,0,0))
  
  
  #TODO: return more
  return(list(NRSI_plot=NRSI_bars,
              sum_of_signatures=sum_of_signatures,
              new_NRSI_df=new_NRSI_df
              
  ))
  
}