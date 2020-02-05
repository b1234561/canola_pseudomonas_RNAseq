### Commands to run DESeq2 on RNA-seq data mapped to Arabidopsis homologs.

rm(list=ls(all.names=TRUE))

library("DESeq2")

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")

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

# Get sequencing lane number for each sample
lane_map <- read.table("mapfile_lane.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE)
lane_map$lane_char <- NA
lane_map[which(lane_map$sequencing_lane == 4), "lane_char"] <- "four"
lane_map[which(lane_map$sequencing_lane == 5), "lane_char"] <- "five"

Bnapus_map <- read.table("tables/Bnapus_merged_func_annot.txt",
                         sep="\t",
                         header=T,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote="")

Bnapus_map_AT <- Bnapus_map[-which(is.na(Bnapus_map$Athaliana_gene)), ]

At_to_descrip <- Bnapus_map_AT[, c("Athaliana_gene", "Athaliana_description")]
rownames(At_to_descrip) <- make.names(At_to_descrip$Athaliana_gene, unique=TRUE)

# Get corresponding condition for each sample from rnaseq_map.
sample_table$condition <- rnaseq_map[sample_table$sample, "Sample._Name"]

# Remove quotations from conditions, convert commands to hyphens,
# and remove replicate number from after C and I.
sample_table$condition <- gsub("\"", "", sample_table$condition)
sample_table$condition <- gsub(", ", "_", sample_table$condition)
sample_table$condition <- gsub("_C\\d_", "_C_", sample_table$condition)
sample_table$condition <- gsub("_I\\d_", "_I_", sample_table$condition)

# Convert all columns to be factors.
sample_table$sample <- as.factor(sample_table$sample)
sample_table$condition <- as.factor(sample_table$condition)
sample_table$file <- as.factor(sample_table$file)
sample_table$lane_char <- as.factor(lane_map[sample_table$sample, "lane_char"])


# Build DESeq2 dataset from AT-mapped mmquant formatted files and then run DESeq2.
sample_table_dataset <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, 
                                                             directory="mmquant_out_AT/", 
                                                             design = ~ lane_char + condition)
sample_dataset <- DESeq(sample_table_dataset, parallel=TRUE)



# Run per day comparisons.
sample_day1_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D1", "S_C_D1"))
sample_day1_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D1", "R_C_D1"))

sample_day3_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D3", "S_C_D3"))
sample_day3_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D3", "R_C_D3"))

sample_day5_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D5", "S_C_D5"))
sample_day5_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D5", "R_C_D5"))



# Shrink log-fold change values.
day1_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D1", "S_C_D1"), 
                                               res=sample_day1_results_SC_SI)

day1_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D1", "R_C_D1"), 
                                               res=sample_day1_results_RC_RI)

day3_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D3", "S_C_D3"), 
                                               res=sample_day3_results_SC_SI)

day3_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D3", "R_C_D3"), 
                                               res=sample_day3_results_RC_RI)

day5_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D5", "S_C_D5"), 
                                               res=sample_day5_results_SC_SI)

day5_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D5", "R_C_D5"), 
                                               res=sample_day5_results_RC_RI)



### Write out intermediate files after shrinking l2fc values.
write.table(day1_results_SC_SI_shrink, file="At_deseq2_outfiles/day1_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day1_results_RC_RI_shrink, file="At_deseq2_outfiles/day1_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day3_results_SC_SI_shrink, file="At_deseq2_outfiles/day3_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day3_results_RC_RI_shrink, file="At_deseq2_outfiles/day3_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day5_results_SC_SI_shrink, file="At_deseq2_outfiles/day5_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day5_results_RC_RI_shrink, file="At_deseq2_outfiles/day5_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


sample_table_dataset <- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, 
                                                             directory="mmquant_out_AT/", 
                                                             design = ~ lane_char + condition)
sample_dataset <- DESeq(sample_table_dataset, parallel=TRUE)



# Also save main DESeq2 objects as RDS files.
saveRDS(object = sample_table_dataset, file = "At_deseq2_outfiles/sample_table.rds")
saveRDS(object = sample_dataset, file = "At_deseq2_outfiles/deseq2.rds")


# Create final dataframes with all samples.
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

combined_deseq2_out_nonNA <- combined_deseq2_out_nonNA[, c("At_gene", "NAME", "s1", "s3", "s5", "r1", "r3", "r5")]


write.table(combined_deseq2_out_nonNA, file="At_deseq2_outfiles/At_homolog_deseq2_log2fold.txt", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
