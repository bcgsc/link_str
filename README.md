# Analysis scripts developed for genotyping STRs in linked-read data

This repository contains Python scripts developed to:
- extract in-repeat repeats (IRRs) using barcode from linked-read alignments ([IRR extraction](irr/README.md))
- estimate sizes of genomic intervals by calculating Jaccard index (JI) of barcode sharing ([Distance estimate](jaccard_index/README.md))

## Dependancies
- NumPy
- Pandas
- pybedtools
- pysam

Author: [Readman Chiu](mailto:rchiu@bcgsc.ca)

:copyright: Canada's Michael Smith Genome Sciences Centre, BC Cancer
