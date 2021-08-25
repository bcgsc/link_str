# Extract in-repeat reads (IRRs) using barcodes

There are 2 steps in the process:

1. Identify barcodes from target region and extract reads given barcodes

 In this step barcodes of molecules that span the target region (and hence potentially capturing the IRRs) are identified and all reads with these barcodes are extracted.
 
 Because this script needs to iterate through every sequence in FASTQ file(s) it may take hours to finish. 

    ```
    python extract_irr.py <BAM> <coord> <out>
    ```
 
 Inputs:
 - `BAM` LongRanger processed BAM for "10x" or "tell_seq", where barcodes were extracted from `BX` tags. For "stlfr" BAMs, barcodes were extracted from read names
 - `coord` "chrom start end". Make sure chrom name is in the same format of the BAM (i.e. with or without "chr")
 - `out` output file
    
 Some important arguments:
 - `--tech` "10x" (default) or "stlfr" or "tell_seq"
 - `--fqs` paired FASTQs, interleaved format accepted for "10x"
 - `--w` window size on each side of the target region to examine alignments. Default: 500 kb
    
 Output columns (tab delimited): 
 - `barcode`
 - `haplotype` 1 or 2 or "na"
 - `read name`
 - `read1 sequence`
 - `read2 sequence`

2. Screen extracted reads for IRRs

 In this step all the reads extracted from the step 1 are examined if they are either **IRR pairs** (both mates are IRRs) or **IRR anchors** (1 mate is IRR and the other partially or not IRR that maps to the target region.
    
 We use [`TRF`](https://tandem.bu.edu/trf/trf.html) to checking sequences if reads have repeats and [`blastn`](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) to check if non-repeat sequences of **IRR anchor** candidates can be mapped to target region.
    
 ```
 python id_irr.py <extraction_output> <coord> <motif> <genome_fasta>
 ```
 Inputs:
 - `extraction output` output from `extract_irr.py`(step 1)
 - `coord` "chrom start end". Make sure chrom name is in the same format of the BAM (i.e. with or without "chr")
 - `motif` repeat motif sequence
 - `genome fasta` path of genome reference fasta - for extracting flanking sequences of target region to determine **IRR anchors**
 
 Some important arguments:
 - `--report` report of final tallies
 - `trf_out` `TRF` will not be re-run when output file (from previous run) provided
 - `--tech` "10x" (default) or "stlfr" or "tell_seq"
 - `--debug` run with `--debug` will not remove `TRF` output
 
 Output report:
 
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

