# beech_micro_communities
![networks](https://github.com/LukDrey/beech_micro_communities/blob/main/combined_network.png)

This repository contains the data and pre-processing pipeline, as well as the R script for the analysis accompanying our paper: 

> Dreyling L, Schmitt I, Dal Grande F. 2022. Tree size drives diversity and community structure of microbial communities on the bark of beech (Fagus sylvatica). Frontiers in Forests and Global Change. doi: [10.3389/ffgc.2022.858382](https://www.frontiersin.org/articles/10.3389/ffgc.2022.858382/abstract).

## Contacts

**Lukas Dreyling**  
Doctoral Candidate  
[E-Mail](mailto:lukas.dreyling@senckenberg.de)  

**Imke Schmitt**  
Principal Investigator  
[E-Mail](mailto:imke.schmitt@senckenberg.de)  

**Francesco Dal Grande**  
Principal Investigator  
[E-Mail](mailto:francesco.dalgrande@unipd.it)  

## Contents

1. [Pre-Processing Pipeline](01_processing_pipeline.txt)
2. [Data](02_Data.zip)
3. [Analysis Script](03_beech_micro_communities.R)

Additionally you need to download the raw reads [here](XXXXXXXXXXXXX).  

If you lack the computing power to process the raw reads, the resulting ASV tables, FASTA files, taxonomy table and metadata are located [here](02_Data.zip).  

## Before starting

### You will need to have the following software installed.

#### Pre-Processing 
* fastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Cutadapt https://cutadapt.readthedocs.io/
* R https://www.r-project.org/
    - dada2 https://benjjneb.github.io/dada2/index.html
    - ShortRead https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html
    - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
* BLASTn https://www.ncbi.nlm.nih.gov/books/NBK279690/
* Seed2 http://www.biomed.cas.cz/mbu/lbwrf/seed/help.php

#### Analysis
* R https://www.r-project.org/
* Rstudio https://www.rstudio.com/
  - here https://here.r-lib.org/
  - tidyverse https://www.tidyverse.org/
  - decontam https://benjjneb.github.io/decontam/
  - phyloseq https://joey711.github.io/phyloseq/
  - LULU https://github.com/tobiasgf/lulu
  - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
  - microbiome https://github.com/microbiome/microbiome
