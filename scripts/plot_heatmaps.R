library(ggplot2)
library(reshape2)
library(scales)

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/tables/")

# Define function to cluster a set of genes.
cluster_genes <- function(df_in, genes_in) {
  genes_in_subset <- genes_in[which((genes_in %in% rownames(df_in)))]
  df_in_dist <- dist(df_in[genes_in_subset,])
  return(hclust(df_in_dist))
}

# Plot subset of all genes in heatmap.
clustered_AT <- read.table("At_homolog_deseq2_log2fold_filt_clustered.txt",
                           header=T, 
                           sep="\t",
                           stringsAsFactors = FALSE,
                           quote="",
                           comment.char = "")

# Remove extra columns and rows.
clustered_AT <- clustered_AT[-1, -c(2,3)]
clustered_AT_to_recluster <- clustered_AT[,-1]
rownames(clustered_AT_to_recluster) <- clustered_AT$At_gene
clustered_AT_to_recluster_hclust <- cluster_genes(clustered_AT_to_recluster, rownames(clustered_AT_to_recluster))

clustered_AT_melt <- melt(clustered_AT, variable.name = "sample")
colnames(clustered_AT_melt) <- c("At_gene", "sample", "log2fold")

clustered_AT_melt_set <- clustered_AT_melt
clustered_AT_melt_set$log2fold[which(clustered_AT_melt_set$log2fold > 2)] <- 2
clustered_AT_melt_set$log2fold[which(clustered_AT_melt_set$log2fold <  -2)] <- -2

clustered_AT_melt_set$At_gene <- factor(clustered_AT_melt_set$At_gene, levels = clustered_AT_to_recluster_hclust$labels[clustered_AT_to_recluster_hclust$order])

ggplot(clustered_AT_melt_set, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme(axis.text.y=element_blank())

##########################################################
### Now make heatmaps for smaller subsets of interest. ###
##########################################################

# Read in full table of log2fold changes.
log_AT <- read.table("At_homolog_deseq2_log2fold.txt",
                    header=T, 
                    sep="\t",
                    stringsAsFactors = FALSE,
                    quote="",
                    comment.char = "")

rownames(log_AT) <- log_AT$At_gene
log_AT <- log_AT[,-2]
log_AT_to_recluster <- log_AT[,-1]

log_AT_melt <- melt(log_AT, variable.name = "sample")
colnames(log_AT_melt) <- c("At_gene", "sample", "log2fold")

log_AT_melt_set <- log_AT_melt
log_AT_melt_set$log2fold[which(log_AT_melt_set$log2fold > 2)] <- 2
log_AT_melt_set$log2fold[which(log_AT_melt_set$log2fold <  -2)] <- -2

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
log_AT_death_htclust <- cluster_genes(log_AT_to_recluster, death_genes$V1)
log_AT_death_order <- log_AT_death_htclust$labels[log_AT_death_htclust$order]
log_AT_melt_set_death_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_death_order),]
log_AT_melt_set_death_subset$At_gene <- factor(log_AT_melt_set_death_subset$At_gene, levels = log_AT_death_order)
ggplot(log_AT_melt_set_death_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana death genes")

# plot growth genes.
log_AT_growth_htclust <- cluster_genes(log_AT_to_recluster, growth_genes$V1)
log_AT_growth_order <- log_AT_growth_htclust$labels[log_AT_growth_htclust$order]
log_AT_melt_set_growth_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_growth_order),]
log_AT_melt_set_growth_subset$At_gene <- factor(log_AT_melt_set_growth_subset$At_gene, levels = log_AT_growth_order)
ggplot(log_AT_melt_set_growth_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme_minimal() +
  ggtitle("A. thaliana growth genes")

# plot immune genes.
log_AT_immune_htclust <- cluster_genes(log_AT_to_recluster, immune_genes$V1)
log_AT_immune_order <- log_AT_immune_htclust$labels[log_AT_immune_htclust$order]
log_AT_melt_set_immune_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_immune_order),]
log_AT_melt_set_immune_subset$At_gene <- factor(log_AT_melt_set_immune_subset$At_gene, levels = log_AT_immune_order)
ggplot(log_AT_melt_set_immune_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme_minimal() +
  ggtitle("A. thaliana immune genes")

# plot primary_metabolism genes.
log_AT_primary_metabolism_htclust <- cluster_genes(log_AT_to_recluster, primary_metabolism_genes$V1)
log_AT_primary_metabolism_order <- log_AT_primary_metabolism_htclust$labels[log_AT_primary_metabolism_htclust$order]
log_AT_melt_set_primary_metabolism_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_primary_metabolism_order),]
log_AT_melt_set_primary_metabolism_subset$At_gene <- factor(log_AT_melt_set_primary_metabolism_subset$At_gene, levels = log_AT_primary_metabolism_order)
ggplot(log_AT_melt_set_primary_metabolism_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme_minimal() +
  ggtitle("A. thaliana primary metabolism genes")

# plot secondary_metabolism genes.
log_AT_secondary_metabolism_htclust <- cluster_genes(log_AT_to_recluster, secondary_metabolism_genes$V1)
log_AT_secondary_metabolism_order <- log_AT_secondary_metabolism_htclust$labels[log_AT_secondary_metabolism_htclust$order]
log_AT_melt_set_secondary_metabolism_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_secondary_metabolism_order),]
log_AT_melt_set_secondary_metabolism_subset$At_gene <- factor(log_AT_melt_set_secondary_metabolism_subset$At_gene, levels = log_AT_secondary_metabolism_order)
ggplot(log_AT_melt_set_secondary_metabolism_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme_minimal() +
  ggtitle("A. thaliana secondary metabolism genes")

# plot UPR genes.
log_AT_UPR_htclust <- cluster_genes(log_AT_to_recluster, UPR_genes$V1)
log_AT_UPR_order <- log_AT_UPR_htclust$labels[log_AT_UPR_htclust$order]
log_AT_melt_set_UPR_subset <- log_AT_melt_set[which(log_AT_melt_set$At_gene %in% log_AT_UPR_order),]
log_AT_melt_set_UPR_subset$At_gene <- factor(log_AT_melt_set_UPR_subset$At_gene, levels = log_AT_UPR_order)
ggplot(log_AT_melt_set_UPR_subset, aes(y=At_gene, x=sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme_minimal() +
  ggtitle("A. thaliana UPR genes")

