# Running BLASTn
The below commands were used to run blastn to annotate Brassica napus genes against Arabidopsis thaliana. Note that NCBI BLAST v2.6.0 was used for the below commands.

### The CDS sequences for both species were downloaded from ENSEMBL v36 (since these were the biggest files they are not part of this repository):
```
cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz
cds/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa.gz
```

### Both files were gunzipped:
`
gunzip cds/*gz 
`

### Made BLAST database for A. thaliana (note that I specified the full path in this case, but normally this would just be installed globally):
`
/home/gavin/local/ncbi-blast-2.6.0+/bin/makeblastdb -in Arabidopsis_thaliana.TAIR10.cds.all.fa  -parse_seqids -dbtype nucl 
`

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
`
/home/gavin/local/ncbi-blast-2.6.0+/bin/blastn -db Arabidopsis_thaliana.TAIR10.cds.all.fa -query Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -out blastn_out_napus_vs_thaliana_evalue0.0001.txt -evalue 0.0001 -outfmt "6" -num_threads 20
`

### Also tested out a higher E-value cut-off, which didn't make a difference:
`
/home/gavin/local/ncbi-blast-2.6.0+/bin/blastn -db Arabidopsis_thaliana.TAIR10.cds.all.fa -query Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -out blastn_out_napus_vs_thaliana_evalue0.1.txt -evalue 0.1 -outfmt "6" -num_threads 20
`

These blast output files were moved to _B.napus\_vs\_A.thaliana\_blast\_out_.

# Parsing functional information for each gene

### The headers of these CDS FASTA files contain information on each gene. This information was parsed into tables using the _parse\_ENSEMBL\_fasta\_header.py_.

```
python3 scripts/parse_ENSEMBL_fasta_header.py -f cds/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -o tables/Brassica_napus.AST_PRJEB5043_v1.cds.info.txt

python3 scripts/parse_ENSEMBL_fasta_header.py -f cds/Arabidopsis_thaliana.TAIR10.cds.all.fa -o tables/Arabidopsis_thaliana.TAIR10.cds.info.txt --one_transcript
```

Note that the --one\_transcript option above was used so that only one transcript was used as a representative for each gene.

