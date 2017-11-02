library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas/")

source("/home/gavin/projects/pseudomonas/canola_pseudomonas/scripts/canola_pseudomonas_R_code.R")

# Read in At summed HTseq files and combine samples and filenames into single DF.
At_files <- grep("HTseq_.+_At_sum.txt", list.files("At_HTseq_files/"), value=TRUE)

samples <- gsub("HTseq_", "", At_files)
samples <- gsub("_At_sum.txt", "", samples)

sample_table <- data.frame(sample=samples, file=At_files, condition=NA, stringsAsFactors=FALSE)

# Read in RNAseq mapping file as well and add condition to DF.
# My edit of this file was to convert "ul" (w special char) -> "microlitre" in header.
rnaseq_map <- read.table("Harvard_output/canola_rnaseq_metadata_edit.txt",
                         sep="\t",
                         header=T,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote="")

rownames(rnaseq_map) <- rnaseq_map$HS.Seq.name

Bnapus_map <- read.table("tables/Bnapus_merged_func_annot.txt",
                         sep="\t",
                         header=T,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote="")

Bnapus_map_AT <- Bnapus_map[-which(is.na(Bnapus_map$Athaliana_gene)),]

At_to_descrip <- Bnapus_map_AT[,c("Athaliana_gene", "Athaliana_description")]
rownames(At_to_descrip) <- make.names(At_to_descrip$Athaliana_gene, unique=TRUE)

# Get corresponding condition for each sample from rnaseq_map.
sample_table$condition <- rnaseq_map[sample_table$sample, "Sample._Name"]

# Remove quotations from conditions, convert commands to hyphens,
# and remove replicate number from after C and I.
sample_table$condition <- gsub("\"", "", sample_table$condition)
sample_table$condition <- gsub(", ", "-", sample_table$condition)
sample_table$condition <- gsub("-C\\d-", "-C-", sample_table$condition)
sample_table$condition <- gsub("-I\\d-", "-I-", sample_table$condition)

# Convert all columns to be factors.
sample_table$sample <- as.factor(sample_table$sample)
sample_table$condition <- as.factor(sample_table$condition)
sample_table$file <- as.factor(sample_table$file)

# Build DESeq2 dataset from AT-mapped HTSeqCount formatted files.
sample_table_datset <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, 
                                         directory="At_HTseq_files", 
                                         design = ~ condition)

# Run DESeq2 on all samples.
sample_deseq2 <- DESeq(sample_table_datset)

# Run day 1 comparisons:
sample_day1_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S-I-D1", "S-C-D1"))
sample_day1_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R-I-D1", "R-C-D1"))

# Run day 3 comparisons:
sample_day3_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S-I-D3", "S-C-D3"))
sample_day3_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R-I-D3", "R-C-D3"))

# Run day 5 comparisons:
sample_day5_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S-I-D5", "S-C-D5"))
sample_day5_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R-I-D5", "R-C-D5"))

day1_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
                                    contrast=c("condition", "S-I-D1", "S-C-D1"), 
                                    res=sample_day1_results_SC_SI)
day1_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "R-I-D1", "R-C-D1"), 
                                       res=sample_day1_results_RC_RI)

day3_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "S-I-D3", "S-C-D3"), 
                                       res=sample_day3_results_SC_SI)
day3_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "R-I-D3", "R-C-D3"), 
                                       res=sample_day3_results_RC_RI)

day5_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "S-I-D5", "S-C-D5"), 
                                       res=sample_day5_results_SC_SI)
day5_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "R-I-D5", "R-C-D5"), 
                                       res=sample_day5_results_RC_RI)

# Create output directory.
dir.create("At_deseq2", showWarnings = FALSE)

write.table(day1_results_SC_SI_shrink, file="At_deseq2/day1_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day1_results_RC_RI_shrink, file="At_deseq2/day1_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day3_results_SC_SI_shrink, file="At_deseq2/day3_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day3_results_RC_RI_shrink, file="At_deseq2/day3_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day5_results_SC_SI_shrink, file="At_deseq2/day5_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day5_results_RC_RI_shrink, file="At_deseq2/day5_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#Combine all deseq2 output tables together.
combined_deseq2_out <- cbind(day1_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             day3_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             day5_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             day1_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                             day3_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                             day5_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE])

colnames(combined_deseq2_out) <- c("s1", "s3", "s5", "r1", "r3", "r5")

# Remove rows that are all NAs
combined_deseq2_out_nonNA <- combined_deseq2_out[-which(rowSums(is.na(combined_deseq2_out)) == 6), ]

combined_deseq2_out_nonNA$At_gene <- rownames(combined_deseq2_out_nonNA)
combined_deseq2_out_nonNA$NAME <- At_to_descrip[rownames(combined_deseq2_out_nonNA), "Athaliana_description"]

combined_deseq2_out_nonNA <- combined_deseq2_out_nonNA[,c("At_gene", "NAME", "s1", "s3", "s5", "r1", "r3", "r5")]

write.table(combined_deseq2_out_nonNA, file="tables/At_homolog_deseq2_log2fold.txt", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Run variance stabilizing transformation blindly to compare samples for QA.
vst_blind <- vst(sample_deseq2, blind=TRUE)

# Plot sample distances.
sampleDists <- dist(t(assay(vst_blind)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst_blind$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

qual_col <- c( "#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99" , "#e31a1c" ,
               "#fdbf6f" , "#ff7f00" , "#cab2d6" , "#6a3d9a" , "#ffff99" , "#b15928" )

# Plot PCA for same samples.
pcaData <- plotPCA(vst_blind, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=qual_col)

# Exploratory plots: (example)
plotMA(day5_results_SC_SI_shrink, ylim=c(-2,2))

dir.create("At_gene_sets", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_de", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_up", showWarnings = FALSE)
dir.create("At_gene_sets/shoot_down", showWarnings = FALSE)
dir.create("At_gene_sets/root_de", showWarnings = FALSE)
dir.create("At_gene_sets/root_up", showWarnings = FALSE)
dir.create("At_gene_sets/root_down", showWarnings = FALSE)

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
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                                results2 = day3_results_SC_SI_shrink,
                                results3 = day5_results_SC_SI_shrink,
                                name1 = "d1",
                                name2 = "d3",
                                name3 = "d5",
                                compare_type = compare,
                                padj_cut = 0.1,
                                l2fc_cut = l2fc,
                                prefix = shoot_prefix)
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_RC_RI_shrink,
                                results2 = day3_results_RC_RI_shrink,
                                results3 = day5_results_RC_RI_shrink,
                                name1 = "d1",
                                name2 = "d3",
                                name3 = "d5",
                                compare_type = compare,
                                padj_cut = 0.1,
                                l2fc_cut = l2fc,
                                prefix = root_prefix)
  }
}


### Re-run DESeq2 to find core genes with shared immune response in both shoots and roots.
sample_table$file <- as.factor(sample_table$file)

sample_table_core <- sample_table

sample_table_core$site <- NA
sample_table_core[grep("S-", sample_table_core$condition), "site"] <- "shoot"
sample_table_core[grep("R-", sample_table_core$condition), "site"] <- "root"
sample_table_core$site <- as.factor(sample_table_core$site)

sample_table_core$day <- NA
sample_table_core[grep("D1", sample_table_core$condition), "day"] <- "day1"
sample_table_core[grep("D3", sample_table_core$condition), "day"] <- "day3"
sample_table_core[grep("D5", sample_table_core$condition), "day"] <- "day5"
sample_table_core$day <- as.factor(sample_table_core$day)

sample_table_core$condition <- as.character(sample_table_core$condition)
sample_table_core$condition <- gsub("^R-", "", sample_table_core$condition)
sample_table_core$condition <- gsub("^S-", "", sample_table_core$condition)
sample_table_core$condition <- gsub("-D\\d$", "", sample_table_core$condition)
sample_table_core$condition <- as.factor(sample_table_core$condition)

sample_table_core_dataset <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table_core, 
                                                  directory="At_HTseq_files", 
                                                  design = ~ day + site + condition)

core_deseq2 <- DESeq(sample_table_core_dataset)

sampleresults_IvsC <- results(core_deseq2, contrast=c("condition", "I", "C"))
sampleresults_IvsC_shrink <- lfcShrink(core_deseq2, 
                                    contrast=c("condition", "I", "C"), 
                                    res=sampleresults_IvsC)

# Write out genes differentially expressed in ALL infected vs control samples.
core_de <- rownames(sampleresults_IvsC_shrink)[which(sampleresults_IvsC_shrink$padj < 0.1 & abs(sampleresults_IvsC_shrink$log2FoldChange) > 2)]
core_up <- rownames(sampleresults_IvsC_shrink)[which(sampleresults_IvsC_shrink$padj < 0.1 & sampleresults_IvsC_shrink$log2FoldChange > 2)]
core_down <- rownames(sampleresults_IvsC_shrink)[which(sampleresults_IvsC_shrink$padj < 0.1 & sampleresults_IvsC_shrink$log2FoldChange < -2)]

write_out_vec(core_de, c("At_gene_sets/At_core_IvsC_de.txt"))
write_out_vec(core_up, c("At_gene_sets/At_core_IvsC_up.txt"))
write_out_vec(core_down, c("At_gene_sets/At_core_IvsC_down.txt"))
