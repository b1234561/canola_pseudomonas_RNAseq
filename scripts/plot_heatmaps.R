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

rownames(ratios_AT) <- ratios_AT$At_gene

# Get info columns.
ratios_AT_info <- ratios_AT[, c(1, 2)]

# Make all names unique.
ratios_AT_info$unique_name <- make.unique(ratios_AT_info$NAME, sep = ".")

# Identify rows that have at least an absolute ratio of 1 in at least 2 samples.
ratio_AT_high_abs_rows <- which(rowSums(abs(ratios_AT[, c(3,4,5,6)]) > 1) >= 2)

# Remove extra column.
ratios_AT <- ratios_AT[, -2]

# Also set any absolute values greater than 2 to have a max abs value of 2 (to make it easier to visualize).
ratios_AT_set <- ratios_AT[, -1]
ratios_AT_set[ratios_AT_set > 2] <- 2
ratios_AT_set[ratios_AT_set < -2] <- -2
ratios_AT_set$At_gene <- ratios_AT$At_gene

# Create heatmap based on filtered set of all genes.
ratios_AT_filt_set <- ratios_AT_set[ratio_AT_high_abs_rows,]

ratios_AT_filt_set_hclust <- cluster_genes(ratios_AT_filt_set[,-which(colnames(ratios_AT_filt_set) == "At_gene")],
                                           rownames(ratios_AT_filt_set))

ratios_AT_filt_set_melt <- melt(ratios_AT_filt_set, variable.name = "sample")
colnames(ratios_AT_filt_set_melt) <- c("At_gene", "Sample", "log2fold")
ratios_AT_filt_set_melt$Gene <- ratios_AT_info[ratios_AT_filt_set_melt$At_gene, "unique_name"]
ratios_AT_filt_set_melt$Gene <- factor(ratios_AT_filt_set_melt$Gene, levels = ratios_AT_info[ratios_AT_filt_set_hclust$labels[ratios_AT_filt_set_hclust$order], "unique_name"])

ggplot(ratios_AT_filt_set_melt, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  theme(axis.text.y=element_blank())

##########################################################
### Now make heatmaps for smaller subsets of interest. ###
##########################################################

# First melt set of ALL genes (unfiltered).
ratios_AT_set_melt <- melt(ratios_AT_set, variable.name = "sample")
colnames(ratios_AT_set_melt) <- c("At_gene", "Sample", "log2fold")
ratios_AT_set_melt$Gene <- ratios_AT_info[ratios_AT_set_melt$At_gene, "unique_name"]


# Read in small gene sets of interest.
death_genes <- read.table("../At_gene_sets/functional_groups_of_interest/Death_AT_genes.txt", header=F, stringsAsFactors = FALSE)
growth_genes <- read.table("../At_gene_sets/functional_groups_of_interest/Growth_AT_genes.txt", header=F, stringsAsFactors = FALSE)
immune_genes <- read.table("../At_gene_sets/functional_groups_of_interest/immune_AT_genes.txt", header=F, stringsAsFactors = FALSE)
primary_metabolism_genes <- read.table("../At_gene_sets/functional_groups_of_interest/primary_metabolism_AT_genes.txt", header=F, stringsAsFactors = FALSE)
secondary_metabolism_genes <- read.table("../At_gene_sets/functional_groups_of_interest/secondary_metabolism_AT_genes.txt", header=F, stringsAsFactors = FALSE)
UPR_genes <- read.table("../At_gene_sets/functional_groups_of_interest/UPR_AT_genes.txt", header=F, stringsAsFactors = FALSE)
ABA_genes <- read.table("../At_gene_sets/functional_groups_of_interest/ABA.txt", header=F, stringsAsFactors = FALSE)
ET_genes <- read.table("../At_gene_sets/functional_groups_of_interest/ET.txt", header=F, stringsAsFactors = FALSE)
JA_genes <- read.table("../At_gene_sets/functional_groups_of_interest/JA.txt", header=F, stringsAsFactors = FALSE)
SA_genes <- read.table("../At_gene_sets/functional_groups_of_interest/SA.txt", header=F, stringsAsFactors = FALSE)

# Convert to uppercase in case any genes are lowercase:
death_genes$V1 <- toupper(death_genes$V1)
growth_genes$V1 <- toupper(growth_genes$V1)
immune_genes$V1 <- toupper(immune_genes$V1)
primary_metabolism_genes$V1 <- toupper(primary_metabolism_genes$V1)
secondary_metabolism_genes$V1 <- toupper(secondary_metabolism_genes$V1)
UPR_genes$V1 <- toupper(UPR_genes$V1)
ABA_genes$V1 <- toupper(ABA_genes$V1)
ET_genes$V1 <- toupper(ET_genes$V1)
JA_genes$V1 <- toupper(JA_genes$V1)
SA_genes$V1 <- toupper(SA_genes$V1)

# plot death genes.
log_AT_death_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], death_genes$V1)
log_AT_death_order <- log_AT_death_htclust$labels[log_AT_death_htclust$order]
log_AT_death_order_names <- ratios_AT_info[log_AT_death_order, "unique_name"]
log_AT_melt_set_death_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% death_genes$V1),]
log_AT_melt_set_death_subset$Gene <- factor(log_AT_melt_set_death_subset$Gene, levels = log_AT_death_order_names)
log_AT_melt_set_death_subset$At_gene <- factor(log_AT_melt_set_death_subset$At_gene, levels = log_AT_death_order)

ggplot(log_AT_melt_set_death_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana death genes")

ggplot(log_AT_melt_set_death_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana death genes") + ylab("Gene")

# plot growth genes.
log_AT_growth_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], growth_genes$V1)
log_AT_growth_order <- log_AT_growth_htclust$labels[log_AT_growth_htclust$order]
log_AT_growth_order_names <- ratios_AT_info[log_AT_growth_order, "unique_name"]
log_AT_melt_set_growth_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% growth_genes$V1),]
log_AT_melt_set_growth_subset$Gene <- factor(log_AT_melt_set_growth_subset$Gene, levels = log_AT_growth_order_names)
log_AT_melt_set_growth_subset$At_gene <- factor(log_AT_melt_set_growth_subset$At_gene, levels = log_AT_growth_order)

ggplot(log_AT_melt_set_growth_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana growth genes")

ggplot(log_AT_melt_set_growth_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana growth genes") + ylab("Gene")


# plot immune genes.
log_AT_immune_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], immune_genes$V1)
log_AT_immune_order <- log_AT_immune_htclust$labels[log_AT_immune_htclust$order]
log_AT_immune_order_names <- ratios_AT_info[log_AT_immune_order, "unique_name"]
log_AT_melt_set_immune_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% immune_genes$V1),]
log_AT_melt_set_immune_subset$Gene <- factor(log_AT_melt_set_immune_subset$Gene, levels = log_AT_immune_order_names)
log_AT_melt_set_immune_subset$At_gene <- factor(log_AT_melt_set_immune_subset$At_gene, levels = log_AT_immune_order)

ggplot(log_AT_melt_set_immune_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana immune genes")

ggplot(log_AT_melt_set_immune_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana immune genes") + ylab("Gene")

# plot primary_metabolism genes.
log_AT_primary_metabolism_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], primary_metabolism_genes$V1)
log_AT_primary_metabolism_order <- log_AT_primary_metabolism_htclust$labels[log_AT_primary_metabolism_htclust$order]
log_AT_primary_metabolism_order_names <- ratios_AT_info[log_AT_primary_metabolism_order, "unique_name"]
log_AT_melt_set_primary_metabolism_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% primary_metabolism_genes$V1),]
log_AT_melt_set_primary_metabolism_subset$Gene <- factor(log_AT_melt_set_primary_metabolism_subset$Gene, levels = log_AT_primary_metabolism_order_names)
log_AT_melt_set_primary_metabolism_subset$At_gene <- factor(log_AT_melt_set_primary_metabolism_subset$At_gene, levels = log_AT_primary_metabolism_order)

ggplot(log_AT_melt_set_primary_metabolism_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana primary metabolism genes")

ggplot(log_AT_melt_set_primary_metabolism_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana primary metabolism genes") + ylab("Gene")


# plot secondary_metabolism genes.
log_AT_secondary_metabolism_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], secondary_metabolism_genes$V1)
log_AT_secondary_metabolism_order <- log_AT_secondary_metabolism_htclust$labels[log_AT_secondary_metabolism_htclust$order]
log_AT_secondary_metabolism_order_names <- ratios_AT_info[log_AT_secondary_metabolism_order, "unique_name"]
log_AT_melt_set_secondary_metabolism_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% secondary_metabolism_genes$V1),]
log_AT_melt_set_secondary_metabolism_subset$Gene <- factor(log_AT_melt_set_secondary_metabolism_subset$Gene, levels = log_AT_secondary_metabolism_order_names)
log_AT_melt_set_secondary_metabolism_subset$At_gene <- factor(log_AT_melt_set_secondary_metabolism_subset$At_gene, levels = log_AT_secondary_metabolism_order)

ggplot(log_AT_melt_set_secondary_metabolism_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana secondary metabolism genes")

ggplot(log_AT_melt_set_secondary_metabolism_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana secondary metabolism genes") + ylab("Gene")


# plot UPR genes.
log_AT_UPR_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], UPR_genes$V1)
log_AT_UPR_order <- log_AT_UPR_htclust$labels[log_AT_UPR_htclust$order]
log_AT_UPR_order_names <- ratios_AT_info[log_AT_UPR_order, "unique_name"]
log_AT_melt_set_UPR_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% UPR_genes$V1),]
log_AT_melt_set_UPR_subset$Gene <- factor(log_AT_melt_set_UPR_subset$Gene, levels = log_AT_UPR_order_names)
log_AT_melt_set_UPR_subset$At_gene <- factor(log_AT_melt_set_UPR_subset$At_gene, levels = log_AT_UPR_order)

ggplot(log_AT_melt_set_UPR_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana UPR genes")

ggplot(log_AT_melt_set_UPR_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana UPR genes") + ylab("Gene")

# plot ABA genes
log_AT_ABA_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], ABA_genes$V1)
log_AT_ABA_order <- log_AT_ABA_htclust$labels[log_AT_ABA_htclust$order]
log_AT_ABA_order_names <- ratios_AT_info[log_AT_ABA_order, "unique_name"]
log_AT_melt_set_ABA_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% ABA_genes$V1),]
log_AT_melt_set_ABA_subset$Gene <- factor(log_AT_melt_set_ABA_subset$Gene, levels = log_AT_ABA_order_names)
log_AT_melt_set_ABA_subset$At_gene <- factor(log_AT_melt_set_ABA_subset$At_gene, levels = log_AT_ABA_order)

ggplot(log_AT_melt_set_ABA_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana ABA genes")

ggplot(log_AT_melt_set_ABA_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana ABA genes") + ylab("Gene")


# plot ET genes
log_AT_ET_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], ET_genes$V1)
log_AT_ET_order <- log_AT_ET_htclust$labels[log_AT_ET_htclust$order]
log_AT_ET_order_names <- ratios_AT_info[log_AT_ET_order, "unique_name"]
log_AT_melt_set_ET_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% ET_genes$V1),]
log_AT_melt_set_ET_subset$Gene <- factor(log_AT_melt_set_ET_subset$Gene, levels = log_AT_ET_order_names)
log_AT_melt_set_ET_subset$At_gene <- factor(log_AT_melt_set_ET_subset$At_gene, levels = log_AT_ET_order)

ggplot(log_AT_melt_set_ET_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana ET genes")

ggplot(log_AT_melt_set_ET_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana ET genes") + ylab("Gene")


# plot JA genes
log_AT_JA_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], JA_genes$V1)
log_AT_JA_order <- log_AT_JA_htclust$labels[log_AT_JA_htclust$order]
log_AT_JA_order_names <- ratios_AT_info[log_AT_JA_order, "unique_name"]
log_AT_melt_set_JA_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% JA_genes$V1),]
log_AT_melt_set_JA_subset$Gene <- factor(log_AT_melt_set_JA_subset$Gene, levels = log_AT_JA_order_names)
log_AT_melt_set_JA_subset$At_gene <- factor(log_AT_melt_set_JA_subset$At_gene, levels = log_AT_JA_order)

ggplot(log_AT_melt_set_JA_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana JA genes")

ggplot(log_AT_melt_set_JA_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana JA genes") + ylab("Gene")

# plot SA genes
log_AT_SA_htclust <- cluster_genes(ratios_AT_set[,-which(colnames(ratios_AT_set) == "At_gene")], SA_genes$V1)
log_AT_SA_order <- log_AT_SA_htclust$labels[log_AT_SA_htclust$order]
log_AT_SA_order_names <- ratios_AT_info[log_AT_SA_order, "unique_name"]
log_AT_melt_set_SA_subset <- ratios_AT_set_melt[which(ratios_AT_set_melt$At_gene %in% SA_genes$V1),]
log_AT_melt_set_SA_subset$Gene <- factor(log_AT_melt_set_SA_subset$Gene, levels = log_AT_SA_order_names)
log_AT_melt_set_SA_subset$At_gene <- factor(log_AT_melt_set_SA_subset$At_gene, levels = log_AT_SA_order)

ggplot(log_AT_melt_set_SA_subset, aes(y=Gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana SA genes")

ggplot(log_AT_melt_set_SA_subset, aes(y=At_gene, x=Sample)) +
  geom_tile(aes(fill = log2fold)) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  ggtitle("A. thaliana SA genes") + ylab("Gene")
