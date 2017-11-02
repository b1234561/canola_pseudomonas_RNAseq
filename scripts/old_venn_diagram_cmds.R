### These were the original commands I was using to make the venn diagrams. 
### I wrote an updated function that also spits out the gene IDs overlapping in set. 
### Also it automatically figures out the sets based on specified cut-offs (see deseq2_ThreeWayVenn_and_set function)


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
