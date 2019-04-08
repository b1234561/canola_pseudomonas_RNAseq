### Commands to compare which genes were identified as DE in root / shoot in PA14 infected samples compared to similar studies in Arabidopsis.

### First compare our PA14 shoot infection samples with the leaf DC3000 paper (Vogel 2016).
### For this analysis comparisons will be made between our day 1 shoot with their day 2 leaf sample, and our day 5 shoot with their day 7 leaf sample.

### Next compare our day 1 PA14 root infection samples with the root WCS417 paper at 6 hours (Stringlis 2018).

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/")

# Function to take in a list of character vectors and return
# a distance matrix of the Jaccard distance between all pairwise comparisons.
char_vec_jaccard_dist <- function(list_in) {
  
  out_dist <- data.frame(matrix(NA, nrow=length(names(list_in)), ncol=length(names(list_in))))
  colnames(out_dist) <- names(list_in)
  rownames(out_dist) <- names(list_in)
  
  # Loop through all elements and perform all pairwise comparisons.
  for(e1 in names(list_in)) {
    for(e2 in names(list_in)) {
      gene_set1 <- list_in[[e1]]
      gene_set2 <- list_in[[e2]]
      out_dist[e1, e2] <- 1 - length(which(gene_set1 %in% gene_set2))/length(unique(c(gene_set1, gene_set2)))
    }
  }
  
  return(as.dist(out_dist))
  
}

shoot_day1_up <- read.table("At_gene_sets/shoot_up/At_shoot_genes_d1_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE, header=FALSE)$V1
shoot_day5_up <- read.table("At_gene_sets/shoot_up/At_shoot_genes_d5_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE, header=FALSE)$V1

shoot_day1_down <- read.table("At_gene_sets/shoot_down/At_shoot_genes_d1_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE, header=FALSE)$V1
shoot_day5_down <- read.table("At_gene_sets/shoot_down/At_shoot_genes_d5_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE, header=FALSE)$V1

root_day1_up <- read.table("At_gene_sets/root_up/At_root_genes_d1_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE, header=FALSE)$V1
root_day5_up <- read.table("At_gene_sets/root_up/At_root_genes_d5_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE, header=FALSE)$V1

root_day1_down <- read.table("At_gene_sets/root_down/At_root_genes_d1_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE, header=FALSE)$V1
root_day5_down <- read.table("At_gene_sets/root_down/At_root_genes_d5_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE, header=FALSE)$V1

# Note that significant cut-offs used below are based on what was reported in the Vogel paper.
vogel_de <- read.table("other_pseudomonas_RNAseq/Vogel_2016/T1_T2_DC3000_DE.txt", header=T, sep="\t", stringsAsFactors = FALSE)
rownames(vogel_de) <- vogel_de$AGI
vogel_de_day2 <- vogel_de[which(abs(vogel_de$T1_Ax_Pst.vs.Ax_CTL_log2FC) > 1 & vogel_de$T1_Ax_Pst.vs.Ax_CTL_FDR < 0.05), ]
vogel_de_day7 <- vogel_de[which(abs(vogel_de$T2_Ax_Pst.vs.Ax_CTL_log2FC) > 1 & vogel_de$T2_Ax_Pst.vs.Ax_CTL_FDR < 0.05), ]

vogel_down_day2 <- vogel_de_day2[which(vogel_de_day2$T1_Ax_Pst.vs.Ax_CTL_log2FC < -1), ]
vogel_up_day2 <- vogel_de_day2[which(vogel_de_day2$T1_Ax_Pst.vs.Ax_CTL_log2FC > 1), ]

vogel_down_day7 <- vogel_de_day7[which(vogel_de_day7$T2_Ax_Pst.vs.Ax_CTL_log2FC < -1), ]
vogel_up_day7 <- vogel_de_day7[which(vogel_de_day7$T2_Ax_Pst.vs.Ax_CTL_log2FC > 1), ]


stringlis_de <- read.table("other_pseudomonas_RNAseq/Stringlis_2018/WCS317_6h_DE.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(stringlis_de) <- stringlis_de$AGI_code
stringlis_up <- stringlis_de[which(stringlis_de$log2.fold.change > 1), ]
stringlis_down <- stringlis_de[which(stringlis_de$log2.fold.change < 1), ]

# Calculate pairwise jaccard index between all groupings.
up_down_sets <- list(Canola_shoot_d1_up=shoot_day1_up,
                Canola_shoot_d5_up=shoot_day5_up,
                Canola_shoot_d1_down=shoot_day1_down,
                Canola_shoot_d5_down=shoot_day5_down,
                
                At_leaf_d2_up=vogel_up_day2$AGI,
                At_leaf_d7_up=vogel_up_day7$AGI,
                At_leaf_d2_down=vogel_down_day2$AGI,
                At_leaf_d7_down=vogel_down_day7$AGI,
                
                Canola_root_d1_up=root_day1_up,
                Canola_root_d5_up=root_day5_up,
                Canola_root_d1_down=root_day1_down,
                Canola_root_d5_down=root_day5_down,
                
                At_root_6h_up=stringlis_up$AGI_code,
                At_root_6h_down=stringlis_down$AGI_code)
                
                           
up_down_sets_jaccard <- char_vec_jaccard_dist(up_down_sets)

up_down_sets_jaccard_hclust <- hclust(up_down_sets_jaccard)
plot(up_down_sets_jaccard_hclust)


de_sets <- list(Canola_shoot_d1_de=c(shoot_day1_up, shoot_day1_down),
                     Canola_shoot_d5_de=c(shoot_day5_up, shoot_day5_down),
                     At_leaf_d2_de=vogel_de_day2$AGI,
                     At_leaf_d7_de=vogel_de_day7$AGI,
                     Canola_root_d1_de=c(root_day1_up, root_day1_down),
                     Canola_root_d5_de=c(root_day5_up, root_day5_down),
                     At_root_6h_de=stringlis_de$AGI_code)


de_sets_jaccard <- char_vec_jaccard_dist(de_sets)

de_sets_jaccard_hclust <- hclust(de_sets_jaccard)

new_order <- de_sets_jaccard_hclust$labels[de_sets_jaccard_hclust$order]

de_sets_jaccard_matrix <- as.matrix(de_sets_jaccard)

library(gplots)
gplots::heatmap.2(de_sets_jaccard_matrix,
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = as.dendrogram(de_sets_jaccard_hclust),
                  Colv = as.dendrogram(de_sets_jaccard_hclust),
                  trace="none",
                  margins =c(15, 15),
                  denscol = "grey",
                  density.info = "density")



# Plot scatterplots of log2fold enrichment of overlapping genes.
early_leaf_overlap <- de_sets[["Canola_shoot_d1_de"]][which(de_sets[["Canola_shoot_d1_de"]]  %in% de_sets[["At_leaf_d2_de"]])]
late_leaf_overlap <- de_sets[["Canola_shoot_d5_de"]][which(de_sets[["Canola_shoot_d5_de"]]  %in% de_sets[["At_leaf_d7_de"]])]

early_root_overlap <- de_sets[["Canola_root_d1_de"]][which(de_sets[["Canola_root_d1_de"]]  %in% de_sets[["At_root_6h_de"]])]

deseq2_out <- read.table("tables/At_homolog_deseq2_log2fold.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, quote="", comment.char="")

par(mfrow=c(1, 3))
cor.test(deseq2_out[early_leaf_overlap, "s1"],
         vogel_de_day2[early_leaf_overlap, "T1_Ax_Pst.vs.Ax_CTL_log2FC"], method="spearman")
plot(deseq2_out[early_leaf_overlap, "s1"],
     vogel_de_day2[early_leaf_overlap, "T1_Ax_Pst.vs.Ax_CTL_log2FC"],
     main="Spearman R=0.282 P=0.055",
     xlab="Canola Shoot Day 1",
     ylab="Arabidopsis Leaf Day 2",
     pch=16)
abline(v=0, h=0, lwd=2, lty=2)


cor.test(deseq2_out[late_leaf_overlap, "s5"],
         vogel_de_day7[late_leaf_overlap, "T2_Ax_Pst.vs.Ax_CTL_log2FC"], method="spearman")
plot(deseq2_out[late_leaf_overlap, "s5"],
     vogel_de_day7[late_leaf_overlap, "T2_Ax_Pst.vs.Ax_CTL_log2FC"],
     xlab="Canola Shoot Day 5",
     ylab="Arabidopsis Leaf Day 7",
     main="Spearman R=0.514 P<2.2e-16",
     pch=16)
abline(v=0, h=0, lwd=2, lty=2)


cor.test(deseq2_out[early_root_overlap, "r1"],
         stringlis_de[early_root_overlap, "log2.fold.change"], method="spearman")
plot(deseq2_out[early_root_overlap, "r1"],
     stringlis_de[early_root_overlap, "log2.fold.change"],
     xlab="Canola Root Day 1",
     ylab="Arabidopsis Root 6h",
     main="Spearman R=0.692 P<2.2e-16",
     pch=16)
abline(v=0, h=0, lwd=2, lty=2)
