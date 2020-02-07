# Figure of concordance between qPCR and RNA-seq. Can make two scatterplots: one for root and one for shoot.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("scripts/canola_pseudomonas_R_code.R")

library(cowplot)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Read in qPCR results and RNA-seq (lfcshrink version)
# ABCG40 = AT1G15520
# BBE4 = AT1G26390
# CYP710A1 = AT2G34500
# RBCSF1 = AT5G38420 (which is called "2B" in A. thaliana rather than "F1", but was the top BLAST hit)

qPCR_At_homologs <- c("AT1G15520", "AT1G26390", "AT2G34500", "AT5G38420")

qPCR <- read.table("tables/qPCR_data.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

qPCR_day1 <- qPCR[which(qPCR$day == 1), ]
qPCR_day3 <- qPCR[which(qPCR$day == 3), ]
qPCR_day5 <- qPCR[which(qPCR$day == 5), ]

root_rna_day1 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(root_rna_day1) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(root_rna_day1) <- qPCR_At_homologs
root_rna_day1$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
root_rna_day1$Day <- "1"

root_day1_lfcshrink <- read.table("At_deseq2_outfiles/day1_results_RC_RI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
root_rna_day1$rnaseq_mean <- root_day1_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
root_rna_day1$rnaseq_se <- root_day1_lfcshrink[qPCR_At_homologs, "lfcSE"]

root_rna_day1$qpcr_mean <- qPCR_day1$root_qPCR_mean
root_rna_day1$qpcr_se <- qPCR_day1$root_qPCR_sd


root_rna_day3 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(root_rna_day3) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(root_rna_day3) <- qPCR_At_homologs
root_rna_day3$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
root_rna_day3$Day <- "3"

root_day3_lfcshrink <- read.table("At_deseq2_outfiles/day3_results_RC_RI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
root_rna_day3$rnaseq_mean <- root_day3_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
root_rna_day3$rnaseq_se <- root_day3_lfcshrink[qPCR_At_homologs, "lfcSE"]

root_rna_day3$qpcr_mean <- qPCR_day3$root_qPCR_mean
root_rna_day3$qpcr_se <- qPCR_day3$root_qPCR_sd


root_rna_day5 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(root_rna_day5) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(root_rna_day5) <- qPCR_At_homologs
root_rna_day5$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
root_rna_day5$Day <- "5"

root_day5_lfcshrink <- read.table("At_deseq2_outfiles/day5_results_RC_RI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
root_rna_day5$rnaseq_mean <- root_day5_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
root_rna_day5$rnaseq_se <- root_day5_lfcshrink[qPCR_At_homologs, "lfcSE"]

root_rna_day5$qpcr_mean <- qPCR_day5$root_qPCR_mean
root_rna_day5$qpcr_se <- qPCR_day5$root_qPCR_sd

root_rna <- rbind(root_rna_day1, root_rna_day3, root_rna_day5)
rownames(root_rna) <- NULL



shoot_rna_day1 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(shoot_rna_day1) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(shoot_rna_day1) <- qPCR_At_homologs
shoot_rna_day1$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
shoot_rna_day1$Day <- "1"

shoot_day1_lfcshrink <- read.table("At_deseq2_outfiles/day1_results_SC_SI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
shoot_rna_day1$rnaseq_mean <- shoot_day1_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
shoot_rna_day1$rnaseq_se <- shoot_day1_lfcshrink[qPCR_At_homologs, "lfcSE"]

shoot_rna_day1$qpcr_mean <- qPCR_day1$shoot_qPCR_mean
shoot_rna_day1$qpcr_se <- qPCR_day1$shoot_qPCR_sd


shoot_rna_day3 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(shoot_rna_day3) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(shoot_rna_day3) <- qPCR_At_homologs
shoot_rna_day3$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
shoot_rna_day3$Day <- "3"

shoot_day3_lfcshrink <- read.table("At_deseq2_outfiles/day3_results_SC_SI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
shoot_rna_day3$rnaseq_mean <- shoot_day3_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
shoot_rna_day3$rnaseq_se <- shoot_day3_lfcshrink[qPCR_At_homologs, "lfcSE"]

shoot_rna_day3$qpcr_mean <- qPCR_day3$shoot_qPCR_mean
shoot_rna_day3$qpcr_se <- qPCR_day3$shoot_qPCR_sd


shoot_rna_day5 <- data.frame(matrix(NA, nrow=4, ncol=6))
colnames(shoot_rna_day5) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
rownames(shoot_rna_day5) <- qPCR_At_homologs
shoot_rna_day5$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
shoot_rna_day5$Day <- "5"

shoot_day5_lfcshrink <- read.table("At_deseq2_outfiles/day5_results_SC_SI_shrink.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
shoot_rna_day5$rnaseq_mean <- shoot_day5_lfcshrink[qPCR_At_homologs, "log2FoldChange"]
shoot_rna_day5$rnaseq_se <- shoot_day5_lfcshrink[qPCR_At_homologs, "lfcSE"]

shoot_rna_day5$qpcr_mean <- qPCR_day5$shoot_qPCR_mean
shoot_rna_day5$qpcr_se <- qPCR_day5$shoot_qPCR_sd

shoot_rna <- rbind(shoot_rna_day1, shoot_rna_day3, shoot_rna_day5)
rownames(shoot_rna) <- NULL

root_scatterplot <- ggplot(root_rna, aes(x=rnaseq_mean, y=qpcr_mean, fill=Day, shape=Gene)) + 
                            geom_point(size=5, color="black") +
                            scale_fill_manual(values=c("white", "black", "grey")) +
                            scale_shape_manual(values=c(21, 22, 23, 24)) +
                            guides(fill = guide_legend(override.aes=list(shape=22))) +
                            geom_errorbar(aes(ymin=qpcr_mean - qpcr_se,
                                              ymax=qpcr_mean + qpcr_se),
                                          width=0.02,
                                          size=0.5) +
                            geom_errorbarh(aes(xmin=rnaseq_mean - rnaseq_se,
                                               xmax=rnaseq_mean + rnaseq_se),
                                           size=0.5) +
                            geom_vline(xintercept = 0, color="grey", linetype="dashed", size=0.5) +
                            geom_hline(yintercept = 0, color="grey", linetype="dashed", size=0.5) +
                            ggtitle("Root") +
                            ylab("RT-qPCR log-ratio") +
                            xlab("RNA-seq log-ratio") +
                            xlim(-5, 10) +
                            ylim(-5, 10) +
                            theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1),
                                  legend.background = element_blank(),
                                  legend.box.background = element_rect(colour = "black"))
  

shoot_scatterplot <- ggplot(shoot_rna, aes(x=rnaseq_mean, y=qpcr_mean, fill=Day, shape=Gene)) + 
  geom_point(size=5, color="black") +
  scale_fill_manual(values=c("white", "black", "grey")) +
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22))) +
  geom_errorbar(aes(ymin=qpcr_mean - qpcr_se,
                    ymax=qpcr_mean + qpcr_se),
                width=0.02,
                size=0.5) +
  geom_errorbarh(aes(xmin=rnaseq_mean - rnaseq_se,
                     xmax=rnaseq_mean + rnaseq_se),
                 size=0.5) +
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=0.5) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=0.5) +
  ggtitle("Shoot") +
  ylab("RT-qPCR log-ratio") +
  xlab("RNA-seq log-ratio") +
  xlim(-5, 10) +
  ylim(-5, 10) +
  theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))


pdf(file = "plots/main/Figure5_qPCR_vs_RNAseq.pdf", width=7.3, height=5, onefile=FALSE)
plot_grid(root_scatterplot, shoot_scatterplot,
          nrow=1,
          labels=c('A', 'B'))
dev.off()

