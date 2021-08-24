# Analysis scripts developed for genotyping STRs in linked-read data

This repository contains Python scripts developed to:
- extract in-repeat repeats (IRRs) using barcode from linked-read alignments
- estimate sizes of genomic intervals by calculating Jaccard index (JI) of barcode sharing

## Dependancies
- NumPy
- Pandas
- pybedtools
- pysam

## Usage
### IRR extraction
1. Identify barcodes of molecules that span target region and extract reads using barcodes

2. Identify IRRs from reads extracted in step 1

### Jaccard index size estimation
1. Profile random locations genome-wide to create database from comparisons in step 2

2. Match JI and barcode tallies of target regions against database to generate size estimates
