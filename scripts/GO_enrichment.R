### Commands to run GO enrichment on all gene sets.

rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(topGO)
library(Hmisc)

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

# First get set of background genes and mappings from At orthologs to GO terms.
gene_background <- rownames(read.table("At_deseq2_outfiles/all_samples/At_homolog_deseq2_log2fold_withall.txt",
                                       header=TRUE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE))

At_gene_to_GO_raw <- read.delim(file = "At_GO/ATH_to_GO.tsv", header=FALSE, sep="\t", stringsAsFactors = FALSE)

At_gene_to_GO <- list()

for(g in gene_background) {
  At_gene_to_GO[[g]] <- At_gene_to_GO_raw[which(At_gene_to_GO_raw$V1 == g), "V2"]
}

run_topGO_classic_enrich <- function(sig_genes, go_map, background, name, ontology_set="BP", min_sig=11, max_annotated=999, min_fdr=10^-3) {
  
  geneList <- factor(as.integer(background %in% sig_genes))
  names(geneList) <- background
  
  topGOdata <- new('topGOdata',
                   description=name,
                   ontology=ontology_set,
                   allGenes=geneList,
                   nodeSize = 10,
                   annot=annFUN.gene2GO,
                   gene2GO = go_map)
  
  classic_test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  
  resultFisher <- getSigGroups(topGOdata, classic_test)
  
  allRes <- GenTable(topGOdata, classic = resultFisher, ranksOf = "classic", topNodes=length(score(resultFisher)))
  allRes$fold <- allRes$Significant/allRes$Expected
  allRes$logfold <- log2(allRes$fold)
  allRes$abslogfold <- abs(allRes$logfold)
  allRes$p <- allRes$classic
  allRes$p[which(allRes$p == "< 1e-30")] <- "1e-30"
  allRes$p <- as.numeric(allRes$p)
  allRes$fdr <- p.adjust(allRes$p, "fdr")
  allRes$log10_fdr <- log10(allRes$fdr)

  # The above table would be useful to save for future reference and
  # the below commands will generate the top hits for visualization.
  allRes_abun <- allRes[which(allRes$Significant >= min_sig & allRes$Annotated <= max_annotated & allRes$fdr < min_fdr), ]
  allRes_abun <- allRes_abun[order(allRes_abun$abslogfold, decreasing = TRUE), ]
  
  return(list(topGOdata=topGOdata,
              resultFisher=resultFisher,
              all_results=allRes,
              filt_results=allRes_abun))
}


# root day 1 up
root_day1_up_sig <- read.table("sig_gene_sets/all_samples/day1_up/At_day1_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

root_day1_up_enrich <- run_topGO_classic_enrich(sig_genes = root_day1_up_sig,
                                                go_map = At_gene_to_GO,
                                                background = gene_background,
                                                name = "root_day1_up")

if(nrow(root_day1_up_enrich$filt_results) > 15) {
  root_day1_up_enrich$filt_results <- root_day1_up_enrich$filt_results[1:15, ]
}

root_day1_up_enrich$filt_results$capitalized <- capitalize(root_day1_up_enrich$filt_results$Term)

root_day1_up_enrich$filt_results$capitalized <- factor(root_day1_up_enrich$filt_results$capitalized,
                                                       levels=rev(root_day1_up_enrich$filt_results$capitalized))

root_day1_up_enrich_barplot <- ggplot(root_day1_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 1 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# root day 1 down
root_day1_down_sig <- read.table("sig_gene_sets/all_samples/day1_down/At_day1_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

root_day1_down_enrich <- run_topGO_classic_enrich(sig_genes = root_day1_down_sig,
                                                  go_map = At_gene_to_GO,
                                                  background = gene_background,
                                                  name = "root_day1_down")

if(nrow(root_day1_down_enrich$filt_results) > 15) {
  root_day1_down_enrich$filt_results <- root_day1_down_enrich$filt_results[1:15, ]
}


root_day1_down_enrich$filt_results$capitalized <- capitalize(root_day1_down_enrich$filt_results$Term)

root_day1_down_enrich$filt_results$capitalized <- factor(root_day1_down_enrich$filt_results$capitalized,
                                                         levels=rev(root_day1_down_enrich$filt_results$capitalized))

root_day1_down_enrich_barplot <- ggplot(root_day1_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 1 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# shoot day 1 up
shoot_day1_up_sig <- read.table("sig_gene_sets/all_samples/day1_up/At_day1_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

shoot_day1_up_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day1_up_sig,
                                                 go_map = At_gene_to_GO,
                                                 background = gene_background,
                                                 name = "shoot_day1_up")

if(nrow(shoot_day1_up_enrich$filt_results) > 15) {
  shoot_day1_up_enrich$filt_results <- shoot_day1_up_enrich$filt_results[1:15, ]
}


shoot_day1_up_enrich$filt_results$capitalized <- capitalize(shoot_day1_up_enrich$filt_results$Term)

shoot_day1_up_enrich$filt_results$capitalized <- factor(shoot_day1_up_enrich$filt_results$capitalized,
                                                        levels=rev(shoot_day1_up_enrich$filt_results$capitalized))

shoot_day1_up_enrich_barplot <- ggplot(shoot_day1_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 1 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# shoot day 1 down
shoot_day1_down_sig <- read.table("sig_gene_sets/all_samples/day1_down/At_day1_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

shoot_day1_down_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day1_down_sig,
                                                   go_map = At_gene_to_GO,
                                                   background = gene_background,
                                                   name = "shoot_day1_down")

if(nrow(shoot_day1_down_enrich$filt_results) > 15) {
  shoot_day1_down_enrich$filt_results <- shoot_day1_down_enrich$filt_results[1:15, ]
}


shoot_day1_down_enrich$filt_results$capitalized <- capitalize(shoot_day1_down_enrich$filt_results$Term)

shoot_day1_down_enrich$filt_results$capitalized <- factor(shoot_day1_down_enrich$filt_results$capitalized,
                                                          levels=rev(shoot_day1_down_enrich$filt_results$capitalized))

shoot_day1_down_enrich_barplot <- ggplot(shoot_day1_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 1 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))








# root day 3 up
root_day3_up_sig <- read.table("sig_gene_sets/all_samples/day3_up/At_day3_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

root_day3_up_enrich <- run_topGO_classic_enrich(sig_genes = root_day3_up_sig,
                                                go_map = At_gene_to_GO,
                                                background = gene_background,
                                                name = "root_day3_up")

if(nrow(root_day3_up_enrich$filt_results) > 15) {
  root_day3_up_enrich$filt_results <- root_day3_up_enrich$filt_results[1:15, ]
}

root_day3_up_enrich$filt_results$capitalized <- capitalize(root_day3_up_enrich$filt_results$Term)

root_day3_up_enrich$filt_results$capitalized <- factor(root_day3_up_enrich$filt_results$capitalized,
                                                       levels=rev(root_day3_up_enrich$filt_results$capitalized))

root_day3_up_enrich_barplot <- ggplot(root_day3_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 3 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# root day 3 down
root_day3_down_sig <- read.table("sig_gene_sets/all_samples/day3_down/At_day3_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

root_day3_down_enrich <- run_topGO_classic_enrich(sig_genes = root_day3_down_sig,
                                                  go_map = At_gene_to_GO,
                                                  background = gene_background,
                                                  name = "root_day3_down")

if(nrow(root_day3_down_enrich$filt_results) > 15) {
  root_day3_down_enrich$filt_results <- root_day3_down_enrich$filt_results[1:15, ]
}


root_day3_down_enrich$filt_results$capitalized <- capitalize(root_day3_down_enrich$filt_results$Term)

root_day3_down_enrich$filt_results$capitalized <- factor(root_day3_down_enrich$filt_results$capitalized,
                                                         levels=rev(root_day3_down_enrich$filt_results$capitalized))

root_day3_down_enrich_barplot <- ggplot(root_day3_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 3 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# shoot day 3 up
shoot_day3_up_sig <- read.table("sig_gene_sets/all_samples/day3_up/At_day3_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

shoot_day3_up_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day3_up_sig,
                                                 go_map = At_gene_to_GO,
                                                 background = gene_background,
                                                 name = "shoot_day3_up")

if(nrow(shoot_day3_up_enrich$filt_results) > 15) {
  shoot_day3_up_enrich$filt_results <- shoot_day3_up_enrich$filt_results[1:15, ]
}


shoot_day3_up_enrich$filt_results$capitalized <- capitalize(shoot_day3_up_enrich$filt_results$Term)

shoot_day3_up_enrich$filt_results$capitalized <- factor(shoot_day3_up_enrich$filt_results$capitalized,
                                                        levels=rev(shoot_day3_up_enrich$filt_results$capitalized))

shoot_day3_up_enrich_barplot <- ggplot(shoot_day3_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 3 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# shoot day 3 down
shoot_day3_down_sig <- read.table("sig_gene_sets/all_samples/day3_down/At_day3_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

shoot_day3_down_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day3_down_sig,
                                                   go_map = At_gene_to_GO,
                                                   background = gene_background,
                                                   name = "shoot_day3_down")

if(nrow(shoot_day3_down_enrich$filt_results) > 15) {
  shoot_day3_down_enrich$filt_results <- shoot_day3_down_enrich$filt_results[1:15, ]
}


shoot_day3_down_enrich$filt_results$capitalized <- capitalize(shoot_day3_down_enrich$filt_results$Term)

shoot_day3_down_enrich$filt_results$capitalized <- factor(shoot_day3_down_enrich$filt_results$capitalized,
                                                          levels=rev(shoot_day3_down_enrich$filt_results$capitalized))

shoot_day3_down_enrich_barplot <- ggplot(shoot_day3_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 3 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))









# root day 5 up
root_day5_up_sig <- read.table("sig_gene_sets/all_samples/day5_up/At_day5_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

root_day5_up_enrich <- run_topGO_classic_enrich(sig_genes = root_day5_up_sig,
                                                go_map = At_gene_to_GO,
                                                background = gene_background,
                                                name = "root_day5_up")

if(nrow(root_day5_up_enrich$filt_results) > 15) {
  root_day5_up_enrich$filt_results <- root_day5_up_enrich$filt_results[1:15, ]
}

root_day5_up_enrich$filt_results$capitalized <- capitalize(root_day5_up_enrich$filt_results$Term)

root_day5_up_enrich$filt_results$capitalized <- factor(root_day5_up_enrich$filt_results$capitalized,
                                                       levels=rev(root_day5_up_enrich$filt_results$capitalized))

root_day5_up_enrich_barplot <- ggplot(root_day5_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 5 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# root day 5 down
root_day5_down_sig <- read.table("sig_gene_sets/all_samples/day5_down/At_day5_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

root_day5_down_enrich <- run_topGO_classic_enrich(sig_genes = root_day5_down_sig,
                                                  go_map = At_gene_to_GO,
                                                  background = gene_background,
                                                  name = "root_day5_down")

if(nrow(root_day5_down_enrich$filt_results) > 15) {
  root_day5_down_enrich$filt_results <- root_day5_down_enrich$filt_results[1:15, ]
}


root_day5_down_enrich$filt_results$capitalized <- capitalize(root_day5_down_enrich$filt_results$Term)

root_day5_down_enrich$filt_results$capitalized <- factor(root_day5_down_enrich$filt_results$capitalized,
                                                         levels=rev(root_day5_down_enrich$filt_results$capitalized))

root_day5_down_enrich_barplot <- ggplot(root_day5_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Root Day 5 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# shoot day 5 up
shoot_day5_up_sig <- read.table("sig_gene_sets/all_samples/day5_up/At_day5_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

shoot_day5_up_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day5_up_sig,
                                                 go_map = At_gene_to_GO,
                                                 background = gene_background,
                                                 name = "shoot_day5_up")

if(nrow(shoot_day5_up_enrich$filt_results) > 15) {
  shoot_day5_up_enrich$filt_results <- shoot_day5_up_enrich$filt_results[1:15, ]
}


shoot_day5_up_enrich$filt_results$capitalized <- capitalize(shoot_day5_up_enrich$filt_results$Term)

shoot_day5_up_enrich$filt_results$capitalized <- factor(shoot_day5_up_enrich$filt_results$capitalized,
                                                        levels=rev(shoot_day5_up_enrich$filt_results$capitalized))

shoot_day5_up_enrich_barplot <- ggplot(shoot_day5_up_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 5 Up-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# shoot day 5 down
shoot_day5_down_sig <- read.table("sig_gene_sets/all_samples/day5_down/At_day5_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

shoot_day5_down_enrich <- run_topGO_classic_enrich(sig_genes = shoot_day5_down_sig,
                                                   go_map = At_gene_to_GO,
                                                   background = gene_background,
                                                   name = "shoot_day5_down")

if(nrow(shoot_day5_down_enrich$filt_results) > 15) {
  shoot_day5_down_enrich$filt_results <- shoot_day5_down_enrich$filt_results[1:15, ]
}


shoot_day5_down_enrich$filt_results$capitalized <- capitalize(shoot_day5_down_enrich$filt_results$Term)

shoot_day5_down_enrich$filt_results$capitalized <- factor(shoot_day5_down_enrich$filt_results$capitalized,
                                                          levels=rev(shoot_day5_down_enrich$filt_results$capitalized))

shoot_day5_down_enrich_barplot <- ggplot(shoot_day5_down_enrich$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dotted") +
  coord_flip() +
  xlab("GO annotation") +
  ylab("Fold enrichment") +
  labs(fill=expression('log'[10]*'(q)')) +
  scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, -3)) +
  ggtitle("Shoot Day 5 Down-regulated") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))





plot_grid(root_day1_up_enrich_barplot, root_day1_down_enrich_barplot, 
          root_day5_up_enrich_barplot, root_day5_down_enrich_barplot,
          shoot_day1_up_enrich_barplot, plot.new(),
          shoot_day5_up_enrich_barplot, shoot_day5_down_enrich_barplot,
          labels=c('A', 'B', 'C', 'D', 'E', '', 'F', 'G'),
          ncol = 2, nrow=4)
    
plot_grid(root_day3_up_enrich_barplot,
          shoot_day3_up_enrich_barplot,
          ncol = 2, nrow=1, labels=c('A', 'B'))
