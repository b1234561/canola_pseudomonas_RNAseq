# Running BLASTn
The below commands were used to run blastn to annotate Brassica napus genes against Arabidopsis thaliana. Note that NCBI BLAST v2.6.0 was used for the below commands.

### The CDS sequences for both species were downloaded from ENSEMBL v36 (since these were the biggest files they are not part of this repository):
```
cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz
cds/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa.gz
```

### Both files were gunzipped:
```
gunzip cds/*gz
```

### Made BLAST database for A. thaliana (note that I specified the full path in this case, but normally this would just be installed globally):
```
/home/gavin/local/ncbi-blast-2.6.0+/bin/makeblastdb -in Arabidopsis_thaliana.TAIR10.cds.all.fa  -parse_seqids -dbtype nucl
```

### The command makeblastdb generate these database files (which were moved to At\_blast\_db):
```
Arabidopsis_thaliana.TAIR10.cds.all.fa.nhr
Arabidopsis_thaliana.TAIR10.cds.all.fa.nin
Arabidopsis_thaliana.TAIR10.cds.all.fa.nog
Arabidopsis_thaliana.TAIR10.cds.all.fa.nsd
Arabidopsis_thaliana.TAIR10.cds.all.fa.nsi
Arabidopsis_thaliana.TAIR10.cds.all.fa.nsq
```

### Ran BLASTn against this database using B. napus CDS file as query. Outputted in tablular format (outfmt 6), limited to hits with an E-value of 0.0001 or lower, ran over 20 threads.
```
/home/gavin/local/ncbi-blast-2.6.0+/bin/blastn -db Arabidopsis_thaliana.TAIR10.cds.all.fa -query Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -out blastn_out_napus_vs_thaliana_evalue0.0001.txt -evalue 0.0001 -outfmt "6" -num_threads 20
```

### Also tested out a higher E-value cut-off, which didn't make a difference:
```
/home/gavin/local/ncbi-blast-2.6.0+/bin/blastn -db Arabidopsis_thaliana.TAIR10.cds.all.fa -query Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -out blastn_out_napus_vs_thaliana_evalue0.1.txt -evalue 0.1 -outfmt "6" -num_threads 20
```

These blast output files were moved to _B.napus\_vs\_A.thaliana\_blast\_out_.

# Parsing functional information for each gene

### The headers of these CDS FASTA files contain information on each gene. This information was parsed into tables using the _parse\_ENSEMBL\_fasta\_header.py_.

```
python3 scripts/parse_ENSEMBL_fasta_header.py -f cds/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -o tables/Brassica_napus.AST_PRJEB5043_v1.cds.info.txt

python3 scripts/parse_ENSEMBL_fasta_header.py -f cds/Arabidopsis_thaliana.TAIR10.cds.all.fa -o tables/Arabidopsis_thaliana.TAIR10.cds.info.txt --one_transcript
```

Note that the --one\_transcript option above was used so that only one transcript was used as a representative for each gene.

# Parsing blast output in R

I then read this output blast file (for evalue=0.0001) and explored the data.
The commands run are in *scripts/explore_B.napus_vs_A.thaliana_blast_output.R*.
64996/101040 (64.33%) B. napus genes hit homologs in A. thaliana. The vast
majority of the time there was only one homolog hit for each B. napus gene (see
distribution of the number of hits in plots/num_unique_At_hits.pdf).

Two cleaned versions of the blast output were created as well. In both cases the
correct blast headers were added. In one case
(*tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean.txt*), only the top
transcript hit from each A. thaliana gene was shown. In the second case, only
the top gene hit was shown (*tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean_tophits.txt*).

### Difference in bitscore between top two hits
To give a measure of confidence for each top blsat hit the difference in bitscore
between the top hit and the next gene hit (when applicable) was outputted to
*tables/blastn_out_napus_vs_thaliana_top2_bitscore_diffs.txt*. The distribution
of these values is shown in *plots/top2_bitscore_diff_hist.pdf*. Clearly about 25%
of gene hits should be interpreted cautiously since they are of similar bitscore
to the next hit (<=100 bits).

### Combining annotations into single file
Next I combined all of the different annotations from the below files into one
table using the R commands in *scripts/merge_func_data_into_single_table.R*.

```
tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean_tophits.txt
tables/Brassica_napus.AST_PRJEB5043_v1.cds.info.txt
tables/Arabidopsis_thaliana.TAIR10.cds.info.txt
B.napus_prior_annotation/Brassica_napus_GO
B.napus_prior_annotation/Brassica_napus_IPR.withdescription
```

The merged annotation file is *tables/Bnapus_merged_func_annot.txt*.

Note that the files in *B.napus_prior_annotation* were downloaded from
http://www.genoscope.cns.fr/brassicanapus/ and were part of
[this paper](http://science.sciencemag.org/content/345/6199/950).

### Combining annotations with differentially expressed gene sets

Finally the raw log fold significant differentially expressed (DE) genes
provided by Jamie Cook (in *RNAseq_DE_genes/original_output/*) were combined
with these annotations in *scripts/merge_func_data_w_DE_genes.R*.

These combined files are in *RNAseq_DE_genes/annotated_output/*.

Note that the first column name was missing in
*RNAseq_DE_genes/original_output/RootDay3.txt* and had to be added by hand.

# Generating 16S Phylogeny for Pseudomonas Species

Downloaded unaligned 16S sequences from RDP 11.5 (current\_Bacteria\_unaligned.fa.gz) and training hidden markov models used by infernal to make alignment (in RDPinfernal1.1Traindata).

### Split sequences corresponding to Pseudomonas species into separate FASTAs

```
mkdir pseudomonas_16S_seq

python scripts/split_fasta_by_species.py -f RDP_release_11.5_align/current_Bacteria_unaligned.fa.gz -o pseudomonas_16S_seq/by_species -g Pseudomonas
```

Resulted in 15224 sequences from 206 Pseudomonas species.

### Filter 16S sequences 

Remove sequences contains Ns, less than 1500 nt, or bigger than 2500 nt.
```
cd pseudomonas_16S_seq

mkdir by_species_filt

parallel -j 60 'vsearch --fastx_filter {} --fastq_maxns 1 --fastq_minlen 1500 --fastq_maxlen 2500 -fastaout by_species_filt/{/.}_filt.fa' ::: by_species/*.fa
```

Retained only 1554 sequences.

### Cluster sequences into centroids of 99% per species

```
mkdir species_centroids

parallel -j 60 'vsearch --cluster_fast {} --id 0.99 --centroids species_centroids/{/.}_centroid.fa' ::: by_species_filt/*.fa
```

Retained 357 16S from 128 species.

### Make alignment with infernal

```
cat species_centroids/*fa >pseudomonas_filt_99id_centroids.fa

cmalign --cpu 70 --dnaout -g --noprob --outformat afa ../RDP_release_11.5_align/RDPinfernal1.1Traindata/bacteria_model.cm pseudomonas_filt_99id_centroids.fa > pseudomonas_filt_99id_centroids_RDP_aligned.afa  
```

Convert "." characters to "-".

```
sed 's/\./-/g' pseudomonas_filt_99id_centroids_RDP_aligned.afa > pseudomonas_filt_99id_centroids_RDP_aligned_no-period.afa
```

### Run RAxML to build phylogeny
```
raxmlHPC-PTHREADS-AVX -s pseudomonas_filt_99id_centroids_RDP_aligned_no-period.afa -n pseudomonas_filt_99id_centroids_RDP_aligned_raxml -m GTRGAMMA -p 51 -f d -T 70

mv RAxML_bestTree.pseudomonas_filt_99id_centroids_RDP_aligned_raxml RAxML_bestTree.pseudomonas_filt_99id_centroids_RDP_aligned_raxml.tre
```

I then downloaded locally and made plots in Dendroscope.
