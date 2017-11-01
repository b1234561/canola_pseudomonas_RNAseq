# This script is used for converting the HTseq count output of reads mapped
# to canola to reads mapped to AT homologs.

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas/")

# First step is to read in all HTseq count files.
htseq_files <- list.files(pattern="_HTSeq.txt", 
                          path="Harvard_output/HT-result", 
                          full.names=TRUE,
                          recursive=TRUE)

# Read in mapping file of Bnapus genes to AT homologs (and other functional)
# categories.
Bnapus_map <- read.table("tables/Bnapus_merged_func_annot.txt",
                         sep="\t",
                         header=T,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote="")

# Limit to Bnapus genes which have an AT homolog (i.e. non-NA).
Bnapus_map_AT <- Bnapus_map[-which(is.na(Bnapus_map$Athaliana_gene)),]

# Set Bnapus gene ids to be rownames.
rownames(Bnapus_map_AT) <- Bnapus_map_AT$Bnapus_gene

# Create output directory.
dir.create("At_HTseq_files", showWarnings = FALSE)

# Loop over all HTseq files and sum up # reads mapped to genes with the same AT homolog.
for(htseq_file in htseq_files) {
  
  # Read in HTseq count file.
  htseq_in <- read.table(htseq_file,
                         header=FALSE,
                         stringsAsFactors = FALSE)
  
  colnames(htseq_in) <- c("Bnapus", "read_count")
  
  # Remove all gene ids not found in AT filtered mapfile.
  htseq_in_filt <- htseq_in[which(htseq_in$Bnapus %in% rownames(Bnapus_map_AT)),]
  
  # Add AT homolog as new column
  htseq_in_filt$At <- as.factor(Bnapus_map_AT[htseq_in_filt$Bnapus, "Athaliana_gene"])
  
  # Change Bnapus gene id to factor.
  htseq_in_filt$Bnapus <- as.factor(htseq_in_filt$Bnapus)
  
  # Sum read counts together based on matching At homolog.
  htseq_in_filt_sum <- aggregate(read_count ~ At, 
                                 data=htseq_in_filt, 
                                 FUN=sum)
  
  # Create output file that ends in "_At_sum.txt" rather than "_HTSeq.txt".
  outfile <- paste("At_HTseq_files",
                   gsub(pattern = "_HTSeq.txt$",
                        "_At_sum.txt", 
                        basename(htseq_file)),
                   sep="/")
  
  write.table(x = htseq_in_filt_sum, 
              file = outfile, 
              col.names = FALSE, 
              quote = FALSE, 
              row.names = FALSE)
  
}
