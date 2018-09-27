setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/")

source("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")

# Load in DESeq2 output RDS objects:
day1_results_SC_SI_shrink <- readRDS(file = "At_deseq2/RDS_files/day1_results_SC_SI_shrink.rds")
day1_results_RC_RI_shrink <- readRDS(file = "At_deseq2/RDS_files/day1_results_RC_RI_shrink.rds")

day5_results_SC_SI_shrink <- readRDS(file = "At_deseq2/RDS_files/day5_results_SC_SI_shrink.rds")
day5_results_RC_RI_shrink <- readRDS(file = "At_deseq2/RDS_files/day5_results_RC_RI_shrink.rds")

dir.create("At_gene_sets", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_de", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_up", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_down", showWarnings = FALSE)
dir.create("At_gene_sets/root_de", showWarnings = FALSE)
dir.create("At_gene_sets/root_up", showWarnings = FALSE)
dir.create("At_gene_sets/root_down", showWarnings = FALSE)

dir.create("plots/venn", showWarnings = FALSE)

# Venn diagrams and gene sets.
comparison_types <- c("de", "up", "down")
l2fc_cutoffs <- c(0, 2)

for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    shoot_prefix = paste("At_gene_sets/shoot_", compare, "/At_shoot_genes", sep = "")
    root_prefix = paste("At_gene_sets/root_", compare, "/At_root_genes", sep = "")
    
    deseq2_TwoWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                              results2 = day5_results_SC_SI_shrink,
                                name1 = "Day1",
                                name2 = "Day5",
                                compare_type = compare,
                                padj_cut = 0.1,
                                l2fc_cut = l2fc,
                                prefix = shoot_prefix)
    
    deseq2_TwoWayVenn_and_set(results1 = day1_results_RC_RI_shrink,
                                results2 = day5_results_RC_RI_shrink,
                                name1 = "Day1",
                                name2 = "Day5",
                                compare_type = compare,
                                padj_cut = 0.1,
                                l2fc_cut = l2fc,
                                prefix = root_prefix)
  }
}


# Also do comparisons in between 
for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    day1_prefix = paste("At_gene_sets/day1_", compare, "/At_day1_genes", sep = "")
    day5_prefix = paste("At_gene_sets/day5_", compare, "/At_day5_genes", sep = "")
    
    deseq2_TwoWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                              results2 = day1_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day1_prefix)
    
    deseq2_TwoWayVenn_and_set(results1 = day5_results_SC_SI_shrink,
                              results2 = day5_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day5_prefix)
  }
}

