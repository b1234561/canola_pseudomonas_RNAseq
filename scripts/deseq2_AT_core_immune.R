### Re-run DESeq2 to find core genes with shared immune response in both shoots and roots.

library("DESeq2")

source("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")
setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/")

sample_table <- readRDS("At_deseq2/RDS_files/sample_table.rds")

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
                                                        directory="mmquant_out_AT", 
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