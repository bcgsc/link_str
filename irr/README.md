# Extract in-repeat reads (IRRs) using barcodes

There are 2 steps in the process:

1. Identify barcodes from the target region and extract reads given barcodes

 	Because this script needs to iterate through every read in the input FASTQ files it may take hours to finish. 

    ```
    python extract_irr.py <BAM> <coord> <out>
    ```
 
	Inputs:
	- `BAM` LongRanger processed BAM for "10x" or "tell_seq" reads, where barcodes were extracted from the `BX` tags. For "stlfr" BAMs, barcodes were extracted from the read names
	- `coord` "chrom start end". Make sure chromosome names agree with input BAM (i.e. with or without "chr")
	- `out` output file
    
	Some important arguments:
	- `--tech` "10x" (default) or "stlfr" or "tell_seq"
	- `--fqs` paired FASTQs, interleaved format accepted only for "10x"
	- `--w` window size on each side of the target region to examine alignments. Default: 500 kb
    
	Output TSV columns: 
 	- `barcode`
 	- `haplotype` 1 or 2 or "na"
 	- `read name`
 	- `read1 sequence` (excluding barcode sequence for 10x reads)
 	- `read2 sequence`

2. Screen extracted reads for IRRs

	In this step all the reads extracted from step 1 are examined if they are either **IRR pairs** (both mates are IRRs) or **IRR anchors** (1 mate is IRR and the other mate has non-repeat sequence that maps to the target region).
    
	We use [TRF](https://tandem.bu.edu/trf/trf.html) to checking sequences if reads have repeats and [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) to test if non-repeat sequences of **IRR anchor** candidates can be mapped to target region.
    
	```
 	python id_irr.py <extraction_output> <coord> <motif> <genome_fasta>
 	```
 	Inputs:
 	- `extraction output` output from `extract_irr.py`(step 1)
 	- `coord` "chrom start end". Make sure chromosome names agree with input BAM (i.e. with or without "chr")
 	- `motif` target repeat motif sequence
 	- `genome fasta` path of genome reference fasta - for extracting flanking sequences of target region to determine **IRR anchors**
 
 	Some important arguments:
 	- `--report` report of final tallies
 	- `trf_out` output TRF file (from previous run on same input), TRF will not be run if provided
 	- `--tech` "10x" (default) or "stlfr" or "tell_seq"
 	- `--debug` will not remove TRF output if this is used
 
 	Output report(`--report`):
 
 	Number of **IRR pairs** and **IRR anchors** will be reported for:
 	- total
 	- each barcode (line started with `bc`)
 	- each haplotype (line started with `hp`)
 
 	File format:
 	```
 	total: <# IRR pairs> <# IRR anchors>
 	bc <barcode>:<# IRR pairs> <# IRR anchors>
 	...
 	hp <haplotype>:<# IRR pairs> <# IRR anchors>
 	...
 	```

