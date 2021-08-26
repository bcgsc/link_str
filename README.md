# Analysis scripts developed for genotyping STRs in linked-read data

This repository contains Python scripts developed to:
- extract in-repeat repeats (IRRs) using barcode from linked-read alignments ([IRR extraction](irr/README.md))
- estimate sizes of genomic intervals by calculating Jaccard index (JI) of barcode sharing ([distance estimate](jaccard_index/README.md))

## Dependancies
- [NumPy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [pybedtools](https://daler.github.io/pybedtools/)
- [pysam](https://github.com/pysam-developers/pysam)
- [TRF](https://tandem.bu.edu/trf/trf.html) (for IRR extraction)
- [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (for IRR extraction)

Author: [Readman Chiu](mailto:rchiu@bcgsc.ca)

:copyright: Canada's Michael Smith Genome Sciences Centre, BC Cancer
