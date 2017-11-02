library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("VennDiagram")
library("colorRamps")

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas/")

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
            sep="\t", row.names=TRUE, col.names=TRUE)
write.table(day1_results_RC_RI_shrink, file="At_deseq2/day1_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=TRUE)

write.table(day3_results_SC_SI_shrink, file="At_deseq2/day3_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=TRUE)
write.table(day3_results_RC_RI_shrink, file="At_deseq2/day3_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=TRUE)

write.table(day5_results_SC_SI_shrink, file="At_deseq2/day5_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=TRUE)
write.table(day5_results_RC_RI_shrink, file="At_deseq2/day5_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=TRUE)


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

# Define function that will draw three way Venn diagram.
ThreeWayVenn <- function(set1, set2, set3, name1, name2, name3) {
      
  venn.plot <- draw.triple.venn( 
                area1=length(set1),
                area2=length(set2),
                area3=length(set3),
                n12=length(which(set1 %in% set2)),
                n23=length(which(set2 %in% set3)),
                n13=length(which(set1 %in% set3)),
                n123=length(which(set1[which(set1 %in% set2)] %in% set3)),
                category=c(name1, name2, name3),
                fill=c("#1f78b4", "#33a02c", "#e31a1c"))

  grid.draw(venn.plot);
  grid.newpage();
}

# Plot Venn diagrams of overlapping At genes DE (and then broken down by up and down-regulated).
de_day1_sc_si <- rownames(day1_results_SC_SI_shrink)[which(day1_results_SC_SI_shrink$padj < 0.1)]
de_day3_sc_si <- rownames(day3_results_SC_SI_shrink)[which(day3_results_SC_SI_shrink$padj < 0.1)]
de_day5_sc_si <- rownames(day5_results_SC_SI_shrink)[which(day5_results_SC_SI_shrink$padj < 0.1)]

up_day1_sc_si <- rownames(day1_results_SC_SI_shrink)[which(day1_results_SC_SI_shrink$padj < 0.1 & day1_results_SC_SI_shrink$log2FoldChange > 0)]
up_day3_sc_si <- rownames(day3_results_SC_SI_shrink)[which(day3_results_SC_SI_shrink$padj < 0.1 & day3_results_SC_SI_shrink$log2FoldChange > 0)]
up_day5_sc_si <- rownames(day5_results_SC_SI_shrink)[which(day5_results_SC_SI_shrink$padj < 0.1 & day5_results_SC_SI_shrink$log2FoldChange > 0)]

down_day1_sc_si <- rownames(day1_results_SC_SI_shrink)[which(day1_results_SC_SI_shrink$padj < 0.1 & day1_results_SC_SI_shrink$log2FoldChange < 0)]
down_day3_sc_si <- rownames(day3_results_SC_SI_shrink)[which(day3_results_SC_SI_shrink$padj < 0.1 & day3_results_SC_SI_shrink$log2FoldChange < 0)]
down_day5_sc_si <- rownames(day5_results_SC_SI_shrink)[which(day5_results_SC_SI_shrink$padj < 0.1 & day5_results_SC_SI_shrink$log2FoldChange < 0)]

de_day1_rc_ri <- rownames(day1_results_RC_RI_shrink)[which(day1_results_RC_RI_shrink$padj < 0.1)]
de_day3_rc_ri <- rownames(day3_results_RC_RI_shrink)[which(day3_results_RC_RI_shrink$padj < 0.1)]
de_day5_rc_ri <- rownames(day5_results_RC_RI_shrink)[which(day5_results_RC_RI_shrink$padj < 0.1)]

up_day1_rc_ri <- rownames(day1_results_RC_RI_shrink)[which(day1_results_RC_RI_shrink$padj < 0.1 & day1_results_RC_RI_shrink$log2FoldChange > 0)]
up_day3_rc_ri <- rownames(day3_results_RC_RI_shrink)[which(day3_results_RC_RI_shrink$padj < 0.1 & day3_results_RC_RI_shrink$log2FoldChange > 0)]
up_day5_rc_ri <- rownames(day5_results_RC_RI_shrink)[which(day5_results_RC_RI_shrink$padj < 0.1 & day5_results_RC_RI_shrink$log2FoldChange > 0)]

down_day1_rc_ri <- rownames(day1_results_RC_RI_shrink)[which(day1_results_RC_RI_shrink$padj < 0.1 & day1_results_RC_RI_shrink$log2FoldChange < 0)]
down_day3_rc_ri <- rownames(day3_results_RC_RI_shrink)[which(day3_results_RC_RI_shrink$padj < 0.1 & day3_results_RC_RI_shrink$log2FoldChange < 0)]
down_day5_rc_ri <- rownames(day5_results_RC_RI_shrink)[which(day5_results_RC_RI_shrink$padj < 0.1 & day5_results_RC_RI_shrink$log2FoldChange < 0)]

grid.newpage();
ThreeWayVenn(de_day1_sc_si, de_day3_sc_si, de_day5_sc_si, "Shoot Day1 DE", "Shoot Day3 DE", "Shoot Day5 DE")
ThreeWayVenn(up_day1_sc_si, up_day3_sc_si, up_day5_sc_si, "Shoot Day1 up", "Shoot Day3 up", "Shoot Day5 up")
ThreeWayVenn(down_day1_sc_si, down_day3_sc_si, down_day5_sc_si, "Shoot Day1 down", "Shoot Day3 down", "Shoot Day5 down")

ThreeWayVenn(de_day1_rc_ri, de_day3_rc_ri, de_day5_rc_ri, "Root Day1 DE", "Root Day3 DE", "Root Day5 DE")
ThreeWayVenn(up_day1_rc_ri, up_day3_rc_ri, up_day5_rc_ri, "Root Day1 up", "Root Day3 up", "Root Day5 up")
ThreeWayVenn(down_day1_rc_ri, down_day3_rc_ri, down_day5_rc_ri, "Root Day1 down", "Root Day3 down", "Root Day5 down")