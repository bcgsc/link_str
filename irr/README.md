# Extract in-repeat reads (IRRs) using barcodes

There are 3 steps in the process:

1. Identify barcodes from the target region and extract reads given barcodes

 	Because this script needs to iterate through every read in the input FASTQs it may take hours to finish. 

    ```
    python extract_irr.py <BAM> <BED> <out> --fqs [fq1.gz fq2.gz]
    ```
	Inputs:
	- `BAM` LongRanger processed BAM for "10x" or "tell_seq" reads, where barcodes were extracted from the `BX` tags. For "stlfr" BAMs, barcodes were extracted from the read names
	- `BED` BED file of target loci. For each line:"chrom start end motif" (TAB-delimited). Make sure chromosome names agree with input BAM (i.e. with or without "chr" prefix)
	- `out` output file
    
	Some important arguments:
	- `--tech` "10x" (default) or "stlfr" or "tell_seq"
	- `--fqs` (required) paired gzipped FASTQs, interleaved format accepted only for "10x"
	- `--w` window size on each side of the target region to examine alignments. Default: 500 kb
    
	Output columns (space-delimited):
	- `coorde` chrom:start-end
 	- `barcode`
 	- `haplotype` 1 or 2 (extracted from the `HP` tags in the input BAM) or "na"
 	- `read name`
 	- `read1 sequence` (excluding barcode sequence for 10x reads)
 	- `read2 sequence`

2. Screen extracted reads for IRRs

	In this step all the reads extracted from step 1 are examined if they are either **IRR pairs** (both mates are IRRs) or **IRR anchors** (1 mate is IRR and the other mate has non-repeat sequence that maps to the target region).
    
	We use [TRF](https://tandem.bu.edu/trf/trf.html) to checking sequences if reads have repeats and [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) to test if non-repeat sequences of **IRR anchor** candidates can be mapped to target region.
    
	```
 	python id_irr.py <extraction_output> <BED> <genome_fasta>
 	```
 	Inputs:
 	- `extraction output` output from `extract_irr.py`(step 1)
 	- `BED` BED file of target loci. For each line:"chrom start end motif" (TAB-delimited). Make sure chromosome names agree with input BAM (i.e. with or without "chr" prefix)
 	- `genome fasta` path of genome reference fasta - for extracting flanking sequences of target region to determine **IRR anchors**
 
 	Some important arguments:
 	- `--report` report of final IRR tallies
 	- `trf_out` TRF output (from previous run on same input), TRF will not be run if provided
 	- `--tech` "10x" (default) or "stlfr" or "tell_seq"
 	- `--debug` will not remove TRF output if this is used
 
 	Output report(`--report`):
 
 	Number of **IRR pairs** and **IRR anchors** will be reported for:
 	- total
 	- each barcode (line started with `bc`)
 	- each haplotype (line started with `hp`)
 
 	File format:
 	```
 	<coord> total: <# IRR pairs> <# IRR anchors>
 	<coord> bc<barcode>:<# IRR pairs> <# IRR anchors>
 	...
 	<coord> hp <haplotype>:<# IRR pairs> <# IRR anchors>
 	...
 	```
	`coord`=chrom:start-end

3. Estimate size

	Interpreted from ExpansionHunter, this step estimate the repeat size from IRR counts, read length, and coverage

	```
	python get_irr_size.py <report from step 2> <BAM> <tech>
	```
	Inputs:
	- report from step 2 containing IRR count info
	- `BAM` same `BAM` used in step 1
	- `tech` "10x" (default) or "stlfr" or "tell_seq"

	additional parameters:
	- `--hp` haplotype to use (1 or 2)
	- `--het` use only half of coverage found because heterozygous allele expected

	Output columns(stdout):
	- `coord` chrom:start-end
	- `IRR_pairs` irr_pair counts
	- `IRR_anchors` irr_anchor counts
	- `IRRs` total irr counts
	- `coverage` local(+/- 1500bp) average coverage of locus using `get_depth.py`
	- `read length` read length (hard-coded: "10x"=128,151; "stlfr"=100; "tell_seq"=146) 
	- `size` size estimate
