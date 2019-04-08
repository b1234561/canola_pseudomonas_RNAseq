library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/")

source("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")

# Read in At summed HTseq files and combine samples and filenames into single DF.
At_files <- grep(".txt$", list.files("mmquant_out_AT/"), value=TRUE)

samples <- gsub("_mmquant_out_AT.txt$", "", At_files)

sample_table <- data.frame(sample=samples, file=At_files, condition=NA, stringsAsFactors=FALSE)

# Read in RNAseq mapping file as well and add condition to DF.
# My edit of this file was to convert "ul" (w special char) -> "microlitre" in header.
rnaseq_map <- read.table("canola_rnaseq_metadata_edit.txt",
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
sample_table$condition <- gsub(", ", "_", sample_table$condition)
sample_table$condition <- gsub("_C\\d_", "_C_", sample_table$condition)
sample_table$condition <- gsub("_I\\d_", "_I_", sample_table$condition)

# Remove D3 samples since we believe these have technical issues.
sample_table <- sample_table[-grep("_D3", sample_table$condition),]

# Convert all columns to be factors.
sample_table$sample <- as.factor(sample_table$sample)
sample_table$condition <- as.factor(sample_table$condition)
sample_table$file <- as.factor(sample_table$file)

# Build DESeq2 dataset from AT-mapped mmquant formatted files.
sample_table_datset <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, 
                                         directory="mmquant_out_AT", 
                                         design = ~ condition)

# Run DESeq2 on all samples.
sample_deseq2 <- DESeq(sample_table_datset)

# Write normalized table with ALL samples for certain plots:
sample_deseq2_counts <- counts(sample_deseq2, normalized=TRUE)
sample_deseq2_counts <- data.frame(sweep(sample_deseq2_counts, 2, colSums(sample_deseq2_counts), FUN="/")) * 100
sample_table$unique <- make.unique(as.character(sample_table$condition))
rownames(sample_table) <- sample_table$sample
colnames(sample_deseq2_counts) <- sample_table[colnames(sample_deseq2_counts), "unique"]
sample_deseq2_counts$At_gene <- rownames(sample_deseq2_counts)
sample_deseq2_counts <- sample_deseq2_counts[, sort(colnames(sample_deseq2_counts))]

write.table(x = sample_deseq2_counts, file="tables/At_homolog_deseq2_percents.txt", col.names=TRUE, row.names=FALSE,
            quote=FALSE, sep="\t")
            

# Run day 1 comparisons:
sample_day1_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S_I_D1", "S_C_D1"))
sample_day1_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R_I_D1", "R_C_D1"))

# Run day 3 comparisons:
#sample_day3_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S_I_D3", "S_C_D3"))
#sample_day3_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R_I_D3", "R_C_D3"))

# Run day 5 comparisons:
sample_day5_results_SC_SI <- results(sample_deseq2, contrast=c("condition", "S_I_D5", "S_C_D5"))
sample_day5_results_RC_RI <- results(sample_deseq2, contrast=c("condition", "R_I_D5", "R_C_D5"))

day1_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
                                    contrast=c("condition", "S_I_D1", "S_C_D1"), 
                                    res=sample_day1_results_SC_SI)
day1_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "R_I_D1", "R_C_D1"), 
                                       res=sample_day1_results_RC_RI)

#day3_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
#                                       contrast=c("condition", "S_I_D3", "S_C_D3"), 
#                                       res=sample_day3_results_SC_SI)
#day3_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
#                                       contrast=c("condition", "R_I_D3", "R_C_D3"), 
#                                       res=sample_day3_results_RC_RI)

day5_results_SC_SI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "S_I_D5", "S_C_D5"), 
                                       res=sample_day5_results_SC_SI)
day5_results_RC_RI_shrink <- lfcShrink(sample_deseq2, 
                                       contrast=c("condition", "R_I_D5", "R_C_D5"), 
                                       res=sample_day5_results_RC_RI)

# Create output directory.
dir.create("At_deseq2", showWarnings = FALSE)

write.table(day1_results_SC_SI_shrink, file="At_deseq2/day1_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day1_results_RC_RI_shrink, file="At_deseq2/day1_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# write.table(day3_results_SC_SI_shrink, file="At_deseq2/day3_AT_deseq2_SC-SI.txt", 
#             sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# write.table(day3_results_RC_RI_shrink, file="At_deseq2/day3_AT_deseq2_RC-RI.txt", 
#             sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day5_results_SC_SI_shrink, file="At_deseq2/day5_AT_deseq2_SC-SI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day5_results_RC_RI_shrink, file="At_deseq2/day5_AT_deseq2_RC-RI.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### Also, write out shrink files, sample_deseq2, and sample_table as RDS files.
dir.create("At_deseq2/RDS_files", showWarnings = FALSE)
saveRDS(object = day1_results_SC_SI_shrink, file = "At_deseq2/RDS_files/day1_results_SC_SI_shrink.rds")
saveRDS(object = day1_results_RC_RI_shrink, file = "At_deseq2/RDS_files/day1_results_RC_RI_shrink.rds")

#saveRDS(object = day3_results_SC_SI_shrink, file = "At_deseq2/RDS_files/day3_results_SC_SI_shrink.rds")
#saveRDS(object = day3_results_RC_RI_shrink, file = "At_deseq2/RDS_files/day3_results_RC_RI_shrink.rds")

saveRDS(object = day5_results_SC_SI_shrink, file = "At_deseq2/RDS_files/day5_results_SC_SI_shrink.rds")
saveRDS(object = day5_results_RC_RI_shrink, file = "At_deseq2/RDS_files/day5_results_RC_RI_shrink.rds")

saveRDS(object = sample_table, file = "At_deseq2/RDS_files/sample_table.rds")
saveRDS(object = sample_deseq2, file = "At_deseq2/RDS_files/sample_deseq2.rds")


#Combine all deseq2 output tables together.
combined_deseq2_out <- cbind(day1_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             #day3_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             day5_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                             day1_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                             #day3_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                             day5_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE])

colnames(combined_deseq2_out) <- c("s1", "s5", "r1","r5")

# Remove rows that are all NAs
combined_deseq2_out_nonNA <- combined_deseq2_out[-which(rowSums(is.na(combined_deseq2_out)) == 4), ]

combined_deseq2_out_nonNA$At_gene <- rownames(combined_deseq2_out_nonNA)
combined_deseq2_out_nonNA$NAME <- At_to_descrip[rownames(combined_deseq2_out_nonNA), "Athaliana_description"]

combined_deseq2_out_nonNA <- combined_deseq2_out_nonNA[,c("At_gene", "NAME", "s1", "s5", "r1", "r5")]

write.table(combined_deseq2_out_nonNA, file="tables/At_homolog_deseq2_log2fold.txt", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
