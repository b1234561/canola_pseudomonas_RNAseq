# Figure of (a) summary heatmap and (b) main venn diagrams

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("scripts/canola_pseudomonas_R_code.R")


library(cowplot)
library(ggplotify)
library(pheatmap)
library(gridExtra)

# Read in log2fold ratios.
ratios_AT <- read.table("At_deseq2_outfiles/At_homolog_deseq2_log2fold.txt",
                        header=TRUE, 
                        sep="\t",
                        stringsAsFactors = FALSE,
                        quote="",
                        comment.char = "")

rownames(ratios_AT) <- ratios_AT$At_gene

# Get info columns and remove from original df.
ratios_AT_info <- ratios_AT[, c(1, 2)]
ratios_AT <- ratios_AT[, -c(1, 2)]

# Make all names unique.
ratios_AT_info$unique_name <- make.unique(ratios_AT_info$NAME, sep = ".")

# Read in significant genes in at least 1 tissue on at least 1 day.
AT_sig_genes <- read.table("sig_gene_sets/root_shoot_de.txt", header=FALSE, stringsAsFactors = FALSE)$V1

# Subset to the rows idnetified above and set any absolute values greater than 2 to have a max abs value of 2 (to make it easier to visualize).
ratios_AT_set <- ratios_AT[AT_sig_genes, ]
ratios_AT_set[ratios_AT_set > 2] <- 2
ratios_AT_set[ratios_AT_set < -2] <- -2

# Re-order so that root samples are first.
ratios_AT_set <- ratios_AT_set[, c("r1", "r3", "r5", "s1", "s3", "s5")]

# Give columns clearer names
colnames(ratios_AT_set) <- c("Day 1 Root", "Day 3 Root", "Day 5 Root", "Day 1 Shoot", "Day 3 Shoot", "Day 5 Shoot")

panelA <- pheatmap(t(ratios_AT_set),
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   clustering_method = "average",
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   show_colnames = FALSE,
                   treeheight_col = 0,
                   gaps_row=c(3, 3, 3))


# Panels B and C - Venn diagrams for each tissue of overall DE genes (lfc > 2) by day

# Panel B - shoot Venn diagram
shoot1 <- read.table("sig_gene_sets/day1_de/At_day1_genes_Shoot_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
shoot3 <- read.table("sig_gene_sets/day3_de/At_day3_genes_Shoot_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
shoot5 <- read.table("sig_gene_sets/day5_de/At_day5_genes_Shoot_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

shoot_n13 = shoot1[which(shoot1 %in% shoot3)]
shoot_n35 = shoot3[which(shoot3 %in% shoot5)]
shoot_n15 = shoot1[which(shoot1 %in% shoot5)]
shoot_n135 = shoot_n13[which(shoot_n13 %in% shoot5)]

grid.newpage()
shoot_venn <- draw.triple.venn( area1=length(shoot1),
                                area2=length(shoot3),
                                area3=length(shoot5),
                                n12=length(shoot_n13),
                                n23=length(shoot_n35),
                                n13=length(shoot_n15),
                                n123=length(shoot_n135),
                                category=c("Day 1", "Day 3", "Day 5"),
                                fill=c("#1f78b4", "#e31a1c", "#33a02c"),
                                cat.fontfamily = rep("sans", 3),
                                fontfamily = rep("sans", 7),
                                cat.cex = rep(0.85, 3),
                                cex = rep(0.85, 7))

shoot_venn <- grid.arrange(gTree(children=shoot_venn), top="Shoot")

# Panel C - root Venn diagram
root1 <- read.table("sig_gene_sets/day1_de/At_day1_genes_Root_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
root3 <- read.table("sig_gene_sets/day3_de/At_day3_genes_Root_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
root5 <- read.table("sig_gene_sets/day5_de/At_day5_genes_Root_de_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

root_n13 = root1[which(root1 %in% root3)]
root_n35 = root3[which(root3 %in% root5)]
root_n15 = root1[which(root1 %in% root5)]
root_n135 = root_n13[which(root_n13 %in% root5)]

grid.newpage()
root_venn <- draw.triple.venn(area1=length(root1),
                              area2=length(root3),
                              area3=length(root5),
                              n12=length(root_n13),
                              n23=length(root_n35),
                              n13=length(root_n15),
                              n123=length(root_n135),
                              category=c("Day 1", "Day 3", "Day 5"),
                              fill=c("#1f78b4", "#e31a1c", "#33a02c"),
                              cat.fontfamily = rep("sans", 3),
                              fontfamily = rep("sans", 7),
                              cat.cex = rep(0.85, 3),
                              cex = rep(0.85, 7))

root_venn <- grid.arrange(gTree(children=root_venn), top="Root")

venn_grid <- plot_grid(root_venn, shoot_venn, labels=c('B', 'C'))

pdf(file = "plots/main/Figure2.pdf", width=7.3, height=7.3, onefile=FALSE)

plot_grid(plot.new(), as.grob(panelA), venn_grid, nrow=3, ncol=1, labels=c('A'), rel_heights = c(0.07, 1, 1))

dev.off()
