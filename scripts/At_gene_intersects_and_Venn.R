### Commands to determine how significant genes intersect between samples and also to generate Venn diagrams.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")


dir.create("sig_gene_sets", showWarnings = FALSE)
dir.create("plots/venn", showWarnings = FALSE)


DE_categories <- c("shoot_up", "shoot_down", "shoot_de",
                   "root_up", "root_down", "root_de",
                   "day1_up", "day1_down", "day1_de",
                   "day3_up", "day3_down", "day3_de",
                   "day5_up", "day5_down", "day5_de")

for(DE_category in DE_categories) {
  dir.create(paste("sig_gene_sets", DE_category, sep="/"), showWarnings = FALSE)
}

# Load in DESeq2 output files.
day1_results_SC_SI_shrink <- read.table(file = "At_deseq2_outfiles/day1_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day1_results_RC_RI_shrink <- read.table(file = "At_deseq2_outfiles/day1_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

day3_results_SC_SI_shrink <- read.table(file = "At_deseq2_outfiles/day3_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day3_results_RC_RI_shrink <- read.table(file = "At_deseq2_outfiles/day3_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

day5_results_SC_SI_shrink <- read.table(file = "At_deseq2_outfiles/day5_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day5_results_RC_RI_shrink <- read.table(file = "At_deseq2_outfiles/day5_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

# Venn diagrams and gene sets.
comparison_types <- c("de", "up", "down")
l2fc_cutoffs <- 2 # c(0, 2)

for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    shoot_prefix = paste("sig_gene_sets/shoot_", compare, "/At_shoot_genes", sep = "")
    root_prefix = paste("sig_gene_sets/root_", compare, "/At_root_genes", sep = "")
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                                results2 = day3_results_SC_SI_shrink,
                              results3 = day5_results_SC_SI_shrink,
                              name1 = "Day1",
                              name2 = "Day3",
                              name3 = "Day5",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = shoot_prefix,
                              plot_outdir="plots/venn")
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_RC_RI_shrink,
                                results2 = day3_results_RC_RI_shrink,
                              results3 = day5_results_RC_RI_shrink,
                              name1 = "Day1",
                              name2 = "Day3",
                              name3 = "Day5",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = root_prefix,
                              plot_outdir="plots/venn")
  }
}


# Also do comparisons in between tissues.
for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    day1_prefix = paste("sig_gene_sets/day1_", compare, "/At_day1_genes", sep = "")
    day3_prefix = paste("sig_gene_sets/day3_", compare, "/At_day3_genes", sep = "")
    day5_prefix = paste("sig_gene_sets/day5_", compare, "/At_day5_genes", sep = "")
    
    deseq2_TwoWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                              results2 = day1_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day1_prefix,
                              plot_outdir="plots/venn")
    
    deseq2_TwoWayVenn_and_set(results1 = day3_results_SC_SI_shrink,
                              results2 = day3_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day3_prefix,
                              plot_outdir="plots/venn")
    
    deseq2_TwoWayVenn_and_set(results1 = day5_results_SC_SI_shrink,
                              results2 = day5_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day5_prefix,
                              plot_outdir="plots/venn")
  }
}

