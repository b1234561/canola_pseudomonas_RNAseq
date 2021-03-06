# Figures of top GO enrichments based on (1) up-regulated and (2) down-regulated genes.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("scripts/canola_pseudomonas_R_code.R")

library(cowplot)
library(ggplot2)
library(Hmisc)
library(stringr)


go_enrichment_barplot <- function(go_output, title, num_to_keep=15, go2ignore=c(), x_max=13) {
  
  # Remove and GO ids that are specified (which would be because they are 100% redundant with another GO in the plot already)
  row2remove <- which(go_output$filt_results$GO.ID %in% go2ignore)
  if(length(row2remove) > 0) {
    go_output$filt_results <- go_output$filt_results[-row2remove, ]
  }
  
  if(nrow(go_output$filt_results) > num_to_keep) {
    go_output$filt_results <- go_output$filt_results[1:num_to_keep, ]
  }
  
  go_output$filt_results$capitalized <- capitalize(go_output$filt_results$Term)
  
  go_output$filt_results$capitalized <- str_wrap(go_output$filt_results$capitalized, 30)
  
  go_output$filt_results$capitalized <- factor(go_output$filt_results$capitalized,
                                               levels=rev(go_output$filt_results$capitalized))
  return(
    ggplot(go_output$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
      geom_bar(stat="identity") +
      geom_hline(yintercept=1, linetype="dotted") +
      ylim(0, x_max) +
      coord_flip() +
      xlab("GO annotation") +
      ylab("Fold enrichment") +
      labs(fill=expression('log'[10]*'(q)')) +
      scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, log10(0.05))) +
      ggtitle(title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.y = element_text(size=10))
  )
}

go_enrichment_results <- readRDS("At_GO/go_enrichment_results.rds")

# A lot of terms are cut-off in the dataframe so this list is a dictionary to correct them.
terms2replace <- list()
terms2replace[["generation of precursor metabolites and ..."]] <- "Generation of precursor metabolites and energy"
terms2replace[["regulation of jasmonic acid mediated sig..."]] <- "Regulation of jasmonic acid mediated signaling pathway"
terms2replace[["benzene-containing compound metabolic pr..."]] <- "Benzene-containing compound metabolic process"
terms2replace[["defense response, incompatible interacti..."]] <- "Defense response, incompatible interaction"
terms2replace[["photosynthesis, light harvesting in phot..."]] <- "Photosynthesis, light harvesting in photosystem I"
terms2replace[["cellular polysaccharide metabolic proces..."]] <- "Cellular polysaccharide metabolic process"
terms2replace[["plant-type cell wall organization or bio..."]] <- "Plant-type cell wall organization or biogenesis"
terms2replace[["external encapsulating structure organiz..."]] <- "External encapsulating structure organization"
terms2replace[["branched-chain amino acid catabolic proc..."]] <- "Branched-chain amino acid catabolic process"
terms2replace[["glutamine family amino acid catabolic pr..."]] <- "Glutamine family amino acid catabolic process"
terms2replace[["regulation of response to water deprivat..."]] <- "Regulation of response to water deprivation"
terms2replace[["glutamine family amino acid catabolic pr..."]] <- "Glutamine family amino acid catabolic process"
terms2replace[["glucosamine-containing compound cataboli..."]] <- "Glucosamine-containing compound catabolic process"
terms2replace[["aromatic amino acid family catabolic pro..."]] <- "Aromatic amino acid family catabolic process"
terms2replace[["phenol-containing compound metabolic pro..."]] <- "Phenol-containing compound metabolic process"
terms2replace[["energy derivation by oxidation of organi..."]] <- "Energy derivation by oxidation of organic compounds"
terms2replace[["glucosamine-containing compound metaboli..."]] <- "Glucosamine-containing compound metabolic process"
terms2replace[["cell wall macromolecule catabolic proces..."]] <- "Cell wall macromolecule catabolic process"
terms2replace[["defense response by callose deposition i..."]] <- "Defense response by callose deposition in cell wall"
terms2replace[["photosynthetic electron transport in pho..."]] <- "Photosynthetic electron transport in photosystem I"
terms2replace[["regulation of photosynthesis, light reac..."]] <- "Regulation of photosynthesis, light reaction"
terms2replace[["unsaturated fatty acid biosynthetic proc..."]] <- "Unsaturated fatty acid biosynthetic process"
terms2replace[["regulation of generation of precursor me..."]] <- "Regulation of generation of precursor metabolites and energy"


for(category in names(go_enrichment_results)) {
  for(str2replace in names(terms2replace)) {
    if(str2replace %in% go_enrichment_results[[category]]$filt_results$Term) {
      matching_row <- which(go_enrichment_results[[category]]$filt_results$Term == str2replace)
      go_enrichment_results[[category]]$filt_results[matching_row, "Term"] <- terms2replace[[str2replace]]
    }
  }
}
  



root_day1_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day1,
                                      title="Root Day 1 Up-regulated",
                                      go2ignore = c("GO:0036294", "GO:0046217", "GO:0052314", "GO:0052315"))

root_day3_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day3,
                                      title="Root Day 3 Up-regulated",
                                      go2ignore = c("GO:0036294", "GO:1900056", "GO:1905622", "GO:0006030", "GO:0006032", "GO:0046348", "GO:1901072"))                                      

root_day5_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day5,
                                      title="Root Day 5 Up-regulated",
                                      go2ignore = c())                                      

shoot_day1_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day1,
                                      title="Shoot Day 1 Up-regulated",
                                      go2ignore = c("GO:0006030", "GO:0006032", "GO:0046348", "GO:1901072", "GO:0042343"))

shoot_day3_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day3,
                                      title="Shoot Day 3 Up-regulated",
                                      go2ignore = c("GO:0015985"))

shoot_day5_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day5,
                                      title="Shoot Day 5 Up-regulated",
                                      go2ignore = c("GO:0052317", "GO:0006030", "GO:0006032", "GO:0009700",
                                                    "GO:0046217", "GO:0046348", "GO:0052314", "GO:0052315",
                                                    "GO:1901072", "GO:0052482", "GO:1902074", "GO:0019759",
                                                    "GO:0019762", "GO:0009074"))




# Downregulated plots
root_day1_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day1,
                                      title="Root Day 1 Down-regulated",
                                      go2ignore = c("GO:0019757", "GO:0019760", "GO:1990066", "GO:0018119", "GO:0042044"),
                                      x_max = 21)

root_day3_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day3,
                                      title="Root Day 3 Down-regulated",
                                      go2ignore = c(), x_max = 21)                                     

root_day5_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day5,
                                      title="Root Day 5 Down-regulated",
                                      go2ignore = c("GO:0045490", "GO:0048765"),
                                      x_max = 21)                                    

shoot_day1_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day1,
                                       title="Shoot Day 1 Down-regulated",
                                       go2ignore = c(), x_max = 21)

shoot_day3_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day3,
                                       title="Shoot Day 3 Down-regulated",
                                       go2ignore = c(), x_max = 21)                                     

shoot_day5_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day5,
                                       title="Shoot Day 5 Down-regulated",
                                       go2ignore = c("GO:0019758", "GO:0019761"), x_max = 21)  

# Write out plots
pdf(file = "plots/main/Figure3_upregulated_GO.pdf", width=18, height=10, onefile=FALSE)
plot_grid(root_day1_up, root_day3_up, root_day5_up, shoot_day1_up, shoot_day3_up, shoot_day5_up,
          nrow=2,
          labels=c('A', 'B', 'C', 'D', 'E', 'F'))
dev.off()


pdf(file = "plots/main/Figure4_downregulated_GO.pdf", width=18, height=6, onefile=FALSE)
plot_grid(root_day1_down, root_day5_down, shoot_day5_down,
          nrow=1,
          labels=c('A', 'B', 'C'))
dev.off()

