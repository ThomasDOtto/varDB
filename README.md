# varDB

This site holds the scripts and the data of a manuscript that will soon be submitted. The data are unpublished, but they can be used to explore biology. Please don't use the data or algorithm as part of a publication until are fully published. 

## Data

In the Data directory are the files of the assembled var genes and an overview table of the different samples. They are part of the pf3k project - https://www.malariagen.net/projects/pf3k. The full metadata can be found here: ftp://ngs.sanger.ac.uk/production/pf3k/release_5/pf3k_release_5_metadata_20170804.txt.gz. In this directory, all raw reads are available ftp://ngs.sanger.ac.uk/production/pf3k/release_5/BAM/.

We provide the var genes of 2400 assembled genomes (we tried all more, but several samples had not enough var genes to be part of the dataset). The file formats are nucleotide (nt) and amino acids (aa) sequences (length cut-off of 1kb). To this dataset, we also provide the domains and subdomains.

From this large dataset, we generated a "normalised" dataset with 60 samples from 12 countries. From this dataset, we further provide some files like the blast comparison and some Gephi files. 

## Assembly algorithms

The used assembly algorithms can be found here. We would like to state that they were implemented to work on an LSF architecture. 


