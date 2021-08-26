# Distance estimate using Jaccard index of barcode sharing

This workflow is designed for estimating distances/sizes of genomic intervals by leveraging the inverse relationship between distance (_d_) and amount of barcode sharing - measured by the Jaccard index (_JI_).

_JI_ is calculated as the ratio of the intersection to the union of barcodes between the 2 flanking regions of a genomic interval.

We obtained a size estimate of a genomic (test) region in 2 steps:

1. Compile a database (model file) composed of _JI_ profiles of a lot of random locations with a range of sizes in the reference genome 

    A _JI_ profile for a given interval corresponds to the tuple: (_NA_, _NB_, _JAB_) where
    - _NA_ = # barcodes spanning the left flanking region
    - _NB_ = # barcodes spanning the right flanking region
    - _JAB_ = _JI_ of the given interval


    ```
    python model_span.py.py <BAM> <Dmin> <Dmax> <s> <m> <w> <out>
    ```
 
    Inputs:
    - `BAM` LongRanger processed BAM for "10x" or "tell_seq", where barcodes were extracted from `BX` tags. For "stlfr" BAMs, barcodes were extracted from read names
    - `Dmin` minimum _d_
    - `Dmax` maximum _d_
    - `s` step size (from `Dmin` to `Dmax`)
    - `m` number of locations per _d_
    - `w` size of flanking region (suggested `1000`)
    - `out` output file
    
    Some important arguments:
    - `--tech` "10x" (default) or "stlfr" or "tell_seq"
    - `--hap` "1" or "2", only use alignments of this haplotpye (specified by `HP` tag in BAM)
    - `--sw` window size on each side of the target region to identify barcodes and determine their spans Default: 200,000
    - `--x_only` pick positions from chrX only
    - `--auto_only` pick positions from authosomes only
    - `--nprocs` number of processes to use for parallelization
    - `--random_all` By default, `m` "seed" positions will be chosen and a number of intervals (based on `Dmin`, `Dmax`, and `s`) with the same starting coordinate will be created for each position. `random_all` will randomize all positions 
    - `--good_barcodes` a text file of barcodes (one-per-line) to only use
    - `--no_cov` don't determine or output coverage for each interval (**RECOMMENDED**)
    
    Output (model file) columns: (tab delimited): 
    - `JAB` - _JI_ at `pos`
    - `NA` # barcodes spanning the left flanking region
    - `NB` # barcodes spanning the right flanking region
    - `d` distance/size of interval
    - `pos` genomic coordinate of interval
    - `covA` coverage of left flanking region (won't be reported if `--no_cov`)
    - `covB` coverage of right flanking region (won't be reported if `--no_cov`)

2. Compute the _JI_ profile of the test region and extract the locations with the closet profile from the database

    This step will identify the intervals characterized in step 1 that have the closest (miniumum distance) _JI_ profiles to the test region. The sizes of these intervals will then be the basis of the final size estimate of the test region:
    - size estimate = median of sizes of closest intervals
    - range = IQR of sizes of closest intervals
    
    The "distance" in _JI_ profile between 2 intervals is computed as the sum of the absolute difference of each member of the tuple.
    
    ```
    python id_irr.py <BAM> <bed> <models> <out>
    ```
    Inputs:
    - `BAM` LongRanger processed BAM for "10x" or "tell_seq", where barcodes were extracted from `BX` tags. For "stlfr" BAMs, barcodes were extracted from read names
    - `bed` test region/intervals in BED format
    - `models` output from step 1 of the same `BAM` (can be more than one e.g. characterized with different interval sizes)
    - `out` output file
 
    Some important arguments:
    - `--hap` "1" or "2", only use alignments of this haplotpye (specified by `HP` tag in BAM)
    - `-E` maximum number intervals to keep from which to generate size estimates
    - `-F` maximum search distance used for finding "closest" intervals.
    - `-C` maximum difference between the coverage profiles (if coverage is taken into consideration). Coverage distance between 2 intervals is taken as the maximum of the difference between `covA` and `covB` of the 2 intervals
    - `--ignore_cov` don't compare coverage profiles (**RECOMMENDED**)
    - `--min_d` - don't consider intervals in database if their _d_ is less than `min_d`
 
    Output columns:
    - `chr` chromosome
    - `posA_ref` start coordinate given
    - `posB_ref` end coordinate given
    - `d_ref` size of interval
    - `gene` NA (deprecated)
    - `posA_exp` adjusted start coordinate (if `d_ref` is less than minimum size in given model)
    - `posB_exp` adjusted end coordinate (if `d_ref` is less than minimum size in given model)
    - `d_exp` size of adjusted interval
    - `covA` coverage of left flanking region
    - `covB` coverage of right flanking region
    - `median` median size of "closest" intervals
    - `Q1` Q1 of "closest" intervals
    - `Q3` Q3 of "closest" intervals
    - `IQR` IQR of "closest" intervals
    - `estimate` final size estimate - `median` with adjustment size added/substracted
    - `bound_lower` upper boundary of estimate - `Q1` with adjustment size added/substracted
    - `bound_upper` lower boundary of estimate - `Q3` with adjustment size added/substracted
    - `model` model file (step 1) used

