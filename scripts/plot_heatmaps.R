library(ggplot2)
library(reshape2)
library(scales)

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/tables/")

# Define function to cluster a set of genes.
cluster_genes <- function(df_in, genes_in) {
  genes_in_subset <- genes_in[which((genes_in %in% rownames(df_in)))]
  df_in_dist <- dist(df_in[genes_in_subset,])
  return(hclust(df_in_dist, method="complete"))
}

# Read in log2fold ratios.
ratios_AT <- read.table("At_homolog_deseq2_log2fold.txt",
                           header=T, 
                           sep="\t",
                           stringsAsFactors = FALSE,
                           quote="",
                           comment.char = "")

# Filter out rows that don't have at least an absolute ratio of 1 in at least 2 samples.
ratio_AT_high_abs_rows <- which(rowSums(abs(ratios_AT[, c(3,4,5,6)]) > 1) >= 2)
ratios_AT_filt <- ratios_AT[ratio_AT_high_abs_rows,]

# Remove extra columns and rows.
ratios_AT_filt <- ratios_AT_filt[, -2]
rownames(ratios_AT_filt) <- ratios_AT_filt$At_gene
ratios_AT_filt_hclust <- cluster_genes(ratios_AT_filt[,-1], rownames(ratios_AT_filt))

ratios_AT_filt_melt <- melt(ratios_AT_filt, variable.name = "sample")
colnames(ratios_AT_filt_melt) <- c("At_gene", "sample", "log2fold")

# Also melt unfiltered table for below:
ratios_AT_melt <- melt(ratios_AT[,-2], variable.name = "sample")
colnames(ratios_AT_melt) <- c("At_gene", "sample", "log2fold")
rownames(ratios_AT) <- ratios_AT$At_gene

ratios_AT_filt_melt_set <- ratios_AT_filt_melt
ratios_AT_filt_melt_set$log2fold[which(ratios_AT_filt_melt_set$log2fold > 2)] <- 2
ratios_AT_filt_melt_set$log2fold[which(ratios_AT_filt_melt_set$log2fold <  -2)] <- -2

ratios_AT_melt_set <- ratios_AT_melt
ratios_AT_melt_set$log2fold[which(ratios_AT_melt_set$log2fold > 2)] <- 2
ratios_AT_melt_set$log2fold[which(ratios_AT_melt_set$log2fold <  -2)] <- -2

ratios_AT_filt_melt_set$At_gene <- factor(ratios_AT_filt_melt_set$At_gene, levels = ratios_AT_filt_hclust$labels[ratios_AT_filt_hclust$order])

ggplot(ratios_AT_filt_melt_set, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme(axis.text.y=element_blank())

##########################################################
### Now make heatmaps for smaller subsets of interest. ###
##########################################################

# Read in small gene sets of interest.
death_genes <- read.table("../At_gene_sets/functional_groups_of_interest/Death_AT_genes.txt", header=F, stringsAsFactors = FALSE)
growth_genes <- read.table("../At_gene_sets/functional_groups_of_interest/Growth_AT_genes.txt", header=F, stringsAsFactors = FALSE)
immune_genes <- read.table("../At_gene_sets/functional_groups_of_interest/immune_AT_genes.txt", header=F, stringsAsFactors = FALSE)
primary_metabolism_genes <- read.table("../At_gene_sets/functional_groups_of_interest/primary_metabolism_AT_genes.txt", header=F, stringsAsFactors = FALSE)
secondary_metabolism_genes <- read.table("../At_gene_sets/functional_groups_of_interest/secondary_metabolism_AT_genes.txt", header=F, stringsAsFactors = FALSE)
UPR_genes <- read.table("../At_gene_sets/functional_groups_of_interest/UPR_AT_genes.txt", header=F, stringsAsFactors = FALSE)

death_genes$V1 <- toupper(death_genes$V1)
growth_genes$V1 <- toupper(growth_genes$V1)
immune_genes$V1 <- toupper(immune_genes$V1)
primary_metabolism_genes$V1 <- toupper(primary_metabolism_genes$V1)
secondary_metabolism_genes$V1 <- toupper(secondary_metabolism_genes$V1)
UPR_genes$V1 <- toupper(UPR_genes$V1)

# plot death genes.
log_AT_death_htclust <- cluster_genes(ratios_AT[,-c(1,2)], death_genes$V1)
log_AT_death_order <- log_AT_death_htclust$labels[log_AT_death_htclust$order]
log_AT_melt_set_death_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_death_order),]
log_AT_melt_set_death_subset$At_gene <- factor(log_AT_melt_set_death_subset$At_gene, levels = log_AT_death_order)
ggplot(log_AT_melt_set_death_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana death genes")

# plot growth genes.
log_AT_growth_htclust <- cluster_genes(ratios_AT[,-c(1,2)], growth_genes$V1)
log_AT_growth_order <- log_AT_growth_htclust$labels[log_AT_growth_htclust$order]
log_AT_melt_set_growth_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_growth_order),]
log_AT_melt_set_growth_subset$At_gene <- factor(log_AT_melt_set_growth_subset$At_gene, levels = log_AT_growth_order)
ggplot(log_AT_melt_set_growth_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana growth genes")


# plot immune genes.
log_AT_immune_htclust <- cluster_genes(ratios_AT[,-c(1,2)], immune_genes$V1)
log_AT_immune_order <- log_AT_immune_htclust$labels[log_AT_immune_htclust$order]
log_AT_melt_set_immune_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_immune_order),]
log_AT_melt_set_immune_subset$At_gene <- factor(log_AT_melt_set_immune_subset$At_gene, levels = log_AT_immune_order)
ggplot(log_AT_melt_set_immune_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana immune genes")

# plot primary_metabolism genes.
log_AT_primary_metabolism_htclust <- cluster_genes(ratios_AT[,-c(1,2)], primary_metabolism_genes$V1)
log_AT_primary_metabolism_order <- log_AT_primary_metabolism_htclust$labels[log_AT_primary_metabolism_htclust$order]
log_AT_melt_set_primary_metabolism_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_primary_metabolism_order),]
log_AT_melt_set_primary_metabolism_subset$At_gene <- factor(log_AT_melt_set_primary_metabolism_subset$At_gene, levels = log_AT_primary_metabolism_order)
ggplot(log_AT_melt_set_primary_metabolism_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana primary metabolism genes")

# plot secondary_metabolism genes.
log_AT_secondary_metabolism_htclust <- cluster_genes(ratios_AT[,-c(1,2)], secondary_metabolism_genes$V1)
log_AT_secondary_metabolism_order <- log_AT_secondary_metabolism_htclust$labels[log_AT_secondary_metabolism_htclust$order]
log_AT_melt_set_secondary_metabolism_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_secondary_metabolism_order),]
log_AT_melt_set_secondary_metabolism_subset$At_gene <- factor(log_AT_melt_set_secondary_metabolism_subset$At_gene, levels = log_AT_secondary_metabolism_order)
ggplot(log_AT_melt_set_secondary_metabolism_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana secondary metabolism genes")

# plot UPR genes.
log_AT_UPR_htclust <- cluster_genes(ratios_AT[,-c(1,2)], UPR_genes$V1)
log_AT_UPR_order <- log_AT_UPR_htclust$labels[log_AT_UPR_htclust$order]
log_AT_melt_set_UPR_subset <- ratios_AT_melt_set[which(ratios_AT_melt_set$At_gene %in% log_AT_UPR_order),]
log_AT_melt_set_UPR_subset$At_gene <- factor(log_AT_melt_set_UPR_subset$At_gene, levels = log_AT_UPR_order)
ggplot(log_AT_melt_set_UPR_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana UPR genes")
