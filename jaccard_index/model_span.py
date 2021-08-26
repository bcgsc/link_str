import random
import pybedtools
import argparse
import pysam
from multiprocessing import Pool
import itertools
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
from operator import itemgetter
import configparser
import os

def get_flanks(start, end, flank_size): 
    left_flank = start - flank_size + 1, start
    right_flank = end, end + flank_size - 1
    return left_flank, right_flank
            
def worker_generate_data(args):
    bam_file, locs, w, sw, good_barcodes_file, hap, no_cov, is_random, tech = args
    bam = None
    if bam_file:
        bam = pysam.Samfile(bam_file, 'rb')

    good_barcodes = None
    if good_barcodes_file is not None:
        good_barcodes = get_good_barcodes(good_barcodes_file)

    data = generate_data(bam, locs, w, sw, good_barcodes, hap, no_cov, is_random, tech)
    
    return data
            
def multi_generate_data(bam_file, locations, w, sw, good_barcodes_file, hap, no_cov, is_random, tech, nprocs):
    batches = []
    size_per_batch = np.ceil(len(locations)/nprocs)

    locs = []
    slotted = 0
    for loc in locations:
        if len(locs) < size_per_batch:
            locs.append(loc)
        else:
            batches.append((bam_file, locs.copy(), w, sw, good_barcodes_file, hap, no_cov, is_random, tech))
            del locs[:]
    if locs:
        batches.append((bam_file, locs.copy(), w, sw, good_barcodes_file, hap, no_cov, is_random, tech))

    for i in range(len(batches)):
        b = batches[i]
            
    pool = Pool(processes=nprocs)
    try:
        batch_results = pool.map(worker_generate_data, batches)
    except KeyboardInterrupt:
        pool.terminate()
        return set()
    else:
        pool.close()
    pool.join()

    return list(itertools.chain(*batch_results))

def extract_barcodes(pos, barcode_bounds):
    barcodes = set()
    
    for bc in sorted(barcode_bounds, key=lambda bc: barcode_bounds[bc][0]):
        if barcode_bounds[bc][0] < pos[1] and barcode_bounds[bc][-1] > pos[0]:
            barcodes.add(bc)
        if barcode_bounds[bc][0] > pos[1]:
            break
            
    return barcodes    

def get_good_barcodes(infile):
    good_barcodes = set()
    with open(infile, 'r') as ff:
        for line in ff:
            good_barcodes.add(line.rstrip())

    return good_barcodes

def calc_jaccard(bc1, bc2):
    union = bc1 | bc2
    intersect = bc1 & bc2
    
    if len(union) == 0:
        return np.nan, np.nan, np.nan, np.nan
    else:
        j = float(len(intersect))/len(union)
        return j, len(bc1), len(bc2), len(intersect)
                
def get_distances(Dmax, Dmin, s):
    distances = []
    for d in range(Dmin, Dmax, s):
        distances.append(d)
        
    return distances

def prep_genome(genome, cfg):
    config = configparser.ConfigParser()
    config.read('/home/rchiu/work/git/tnr/genomes.cfg')

    if not genome in config.sections() or not 'chrom_sizes' in config[genome]:
        sys.exit('cannot find {} or "chrom_sizes" for {} in {}'.format(genome, genome, cfg))
    
    chrom_sizes = None
    try:
        chrom_sizes = pybedtools.BedTool(config[genome]['chrom_sizes'])
    except:
        sys.exit('cannot extract chromosome sizes from given input {}'.format(config[genome]['chrom_sizes']))

    gaps = None
    try: 
        gaps = pybedtools.BedTool(config[genome]['gaps'])
    except:
        sys.exit('cannot extract gap sizes from given input {}'.format(config[genome]['gaps']))

    segdups = None
    try:
        segdups = pybedtools.BedTool(config[genome]['segdups'])
    except:
        sys.exit('cannot extract segdups from given input {}'.format(config[genome]['segdups'])) 

    return chrom_sizes, gaps, segdups

def pick_locations(segments_list, Dmin, Dmax, s, m, w, use_seeds=True):
    locs = []
    if use_seeds:
        seeds = []
        while len(seeds) < m:
            still_needed = m - len(seeds)
            if still_needed >= len(segments_list):
                segments_pool = list(range(len(segments_list)))
            else:
                segments_pool = random.sample(list(range(len(segments_list))), still_needed)

            for i in segments_pool:
                chrom, start, end = segments_list[i]
                seed = random.randint(start, end)
                while (seed - w) < start or (seed + Dmax + w - 1) > end:
                    seed = random.randint(start, end)
                seeds.append((chrom, seed))
            
        for chrom, start in seeds:
            for size in range(Dmin, Dmax + 2, s):
                locs.append((chrom, start, start + size))
    
    else:
        n = (int((Dmax + 2 - Dmin)/s) + 1) * m
        sizes = list(range(Dmin, Dmax + 2, s)) * m
        random.shuffle(sizes)

        while len(locs) < n:
            still_needed = n - len(locs)
            if still_needed >= len(segments_list):
                segments_pool = list(range(len(segments_list)))
            else:
                segments_pool = random.sample(list(range(len(segments_list))), still_needed)

            for i in segments_pool:
                chrom, start, end = segments_list[i]
                size = sizes[len(locs)]
                loc = random.randint(start, end)
                while (loc - w) < start or (loc + size + w - 1) > end:
                    loc = random.randint(start, end)
                locs.append((chrom, loc, loc + size))

    locs_sorted = sorted(locs, key=itemgetter(0,1))
    #for loc in locs_sorted:
    #    print('ee', loc, loc[2] - loc[1])
    return locs_sorted

def select_locations(m, Dmin, Dmax, s, w, sw, genome, cfg, random_all, no_rmsk=False, no_segdup=False, x_only=False, auto_only=False):
    chrom_sizes, gaps, segdups = prep_genome(genome, cfg)
    
    n = (int((Dmax + 2 - Dmin)/s) + 1) * m
    sizes = list(range(Dmin, Dmax + 2, s)) * m
    random.shuffle(sizes)
    min_size = Dmax + 2 * w
    
    segments = chrom_sizes
    if gaps:
        segments = segments.subtract(gaps)
        
    if no_segdup:
        segments = segments.subtract(segdups)
        
    # filter by size
    segments = segments.filter(lambda s: len(s) >= min_size)
    
    chrom_sizes_dict = {}
    for c in chrom_sizes:
        chrom_sizes_dict[str(c[0])] = int(c[2])
        
    # convert bed object to list
    segments_list = []
    for segment in segments:
        if x_only and 'X' not in str(segment[0]).upper():
            continue
        if auto_only and not str(segment[0]).upper().replace('CHR', '').isdigit():
            continue

        if int(segment[1]) < sw or int(segment[2]) > chrom_sizes_dict[str(segment[0])] - sw:
            #print('skip', segment[0], segment[1], segment[2], chrom_sizes_dict[str(segment[0])])
            continue
        
        segments_list.append((str(segment[0]), int(segment[1]), int(segment[2])))

    locs_sorted = pick_locations(segments_list, Dmin, Dmax, s, m, w, not random_all)
    return locs_sorted

def extract_barcodes_simple(bam, chrom, pos, hap=None, good_barcodes=None):
    barcodes = set()
    for aln in bam.fetch(chrom, pos[0], pos[1]):
        if aln.is_unmapped:
            continue
        if not aln.has_tag('BX'):
            continue
        bc = aln.get_tag('BX').split('-')[0]

        if good_barcodes is not None and bc not in good_barcodes:
            continue

        if hap is not None:
            if aln.has_tag('HP'):
                hp = aln.get_tag('HP')
                if hp != hap:
                    continue
            else:
                continue
        barcodes.add(bc)

def extract_bc_aln(aln, tech='10x'):
    bc = None
    if tech.lower() == '10x':
        if aln.has_tag('BX'):
            bc = aln.get_tag('BX').split('-')[0]
    elif tech.lower() == 'stlfr':
        if '#' in aln.query_name:
            bc = aln.query_name.split('#')[1]

    return bc

def get_all_barcode_bounds(bam, chrom, pos, tech='10x', sw=200000, min_len=1000, hap=None, good_barcodes=None):
    spans = defaultdict(list)
    bounds = {}
    
    for aln in bam.fetch(chrom, pos[0]-sw, pos[1]+sw):
        if aln.is_unmapped:
            continue
        bc = extract_bc_aln(aln, tech=tech)
        if bc is None:
            continue
        '''
        if aln.is_unmapped:
            continue
        if not aln.has_tag('BX'):
            continue
        bc = aln.get_tag('BX').split('-')[0]
        '''
        if good_barcodes is not None and bc not in good_barcodes:
            continue
        
        if hap is not None:
            if aln.has_tag('HP'):
                hp = aln.get_tag('HP')
                if hp != hap:
                    continue
            else:
                continue
        
        spans[bc].extend([aln.reference_start, aln.reference_end])
        
    for bc in spans.keys():
        coords_sorted = sorted(spans[bc], key=int)

        if coords_sorted[-1] - coords_sorted[0] < min_len:
            continue
        
        bounds[bc] = coords_sorted[0], coords_sorted[-1]
        
    return bounds

def group_by_seed(locs):
    groups = defaultdict(list)
    for loc in locs:
        groups[(loc[0], loc[1])].append(loc)

    return groups

def generate_data(bam, locations, w, sw, good_barcodes, hap, no_cov, is_random, tech):
    barcode_bounds = {}
    cov1s = {}
    if not is_random:
        loc_groups = group_by_seed(locations)
        for chrom, start in loc_groups:
            end = max([loc[2] for loc in loc_groups[(chrom, start)]])
            flanks = get_flanks(start, end, w)
            barcode_bounds[(chrom, start)] = get_all_barcode_bounds(bam, chrom, [flanks[0][0], flanks[1][1]],
                                                                    sw=sw,
                                                                    hap=hap,
                                                                    good_barcodes=good_barcodes,
                                                                    tech=tech)
            if not no_cov:
                cov1s[(chrom, start)] = get_coverage(bam, chrom, flanks[0][0], flanks[0][1], hap)

    data = []
    for chrom, start, end in locations:
        bb = None
        if barcode_bounds:
            bb = barcode_bounds[(chrom, start)]
        cov1 = None
        if cov1s:
            cov1 = cov1s[(chrom, start)]
        cols = generate_data_per_location(bam, chrom, start, end, w, sw, good_barcodes, hap, no_cov, barcode_bounds=bb, cov1=cov1, tech=tech)
        if cols is not None:
            data.append(cols)
            
    return data

def get_coverage(bam, chrom, start, end, hap):
    cov = []
    for pc in bam.pileup(chrom, start, end):
        if hap is None:
            cov.append(pc.n)
        else:
            n = 0
            for pr in pc.pileups:
                if not pr.alignment.has_tag('HP') or pr.alignment.get_tag('HP') != hap:
                    continue
                n += 1
            cov.append(n)
    
    if cov:
        return np.average(cov)
    else:
        return None

def generate_data_per_location(bam, chrom, start, end, w, sw, good_barcodes, hap, no_cov, barcode_bounds=None, cov1=None, tech=None):
    data = []
    
    flanks = get_flanks(start, end, w)

    if barcode_bounds is None:
        barcode_bounds = get_all_barcode_bounds(bam, chrom, [flanks[0][0], flanks[1][1]],
                                                sw=sw, 
                                                hap=hap, 
                                                good_barcodes=good_barcodes,
                                                tech=tech)
    
    if not no_cov:
        if cov1 is None:
            cov1 = get_coverage(bam, chrom, flanks[0][0], flanks[0][1], hap)
        cov2 = get_coverage(bam, chrom, flanks[1][0], flanks[1][1], hap)

        if cov1 is None or cov2 is None:
            print('bad', chrom, start, end, cov1, cov2)
            return None

    bc1 = extract_barcodes(flanks[0], barcode_bounds)
    bc2 = extract_barcodes(flanks[1], barcode_bounds)
    #bc1 = extract_barcodes_simple(bam, chrom, flanks[0])
    #bc2 = extract_barcodes_simple(bam, chrom, flanks[1])

    j, bc_left, bc_right, intersect = calc_jaccard(bc1, bc2)
    cols = [j, bc_left, bc_right, intersect, end-start, '{}:{}-{}'.format(chrom, start, end)]
    if not no_cov:
        cols.extend(['{:.2f}'.format(cov1), '{:.2f}'.format(cov2)])

    return cols
            
def output(data, outfile, params, no_cov):
    cols = ['JAB', 'NA', 'NB', 'NAB', 'd', 'pos']
    if not no_cov:
        cols.extend(['covA', 'covB'])
    df = pd.DataFrame(data, columns=cols)
    if not no_cov:
        df.dropna(subset=['JAB', 'NA', 'NB', 'NAB', 'covA', 'covB'], inplace=True)
    else:
        df.dropna(subset=['JAB', 'NA', 'NB', 'NAB'], inplace=True)

    with open(outfile, 'w') as out:
        out.write('# %s\n' % params)
        df.sort_values(by=['d', 'pos']).to_csv(out, sep='\t', index=False)
            
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("Dmin", type=int, help="minimum distance. Default:0", default=0)
    parser.add_argument("Dmax", type=int, help="maximum distance")
    parser.add_argument("s", type=int, help="step size")
    parser.add_argument("m", type=int, help="number of locations per distance")
    parser.add_argument("w", type=int, help="flank window size")
    parser.add_argument("out", type=str, help="name of output file")
    parser.add_argument("--tech", type=str, help="10x(default), stlfr", default='10x')
    parser.add_argument("--genome", type=str, help="hg38 or hg19 or b37")
    parser.add_argument("--cfg", type=str, help="config file specifying paths of chr_sizes and gaps files")
    parser.add_argument("--sw", type=int, default=200000, help="barcode screening window size. Default=200000")
    parser.add_argument("--good_barcodes", type=str, help="list of good barcodes")
    parser.add_argument("--hap", type=int, help="haplotype 1 or 2")
    parser.add_argument("--no_segdup", action='store_true', help="do not pick from segdup regions")
    parser.add_argument("--nprocs", type=int, help="number of processes", default=1)
    parser.add_argument("--x_only", action='store_true', help="just chromosome X")
    parser.add_argument("--auto_only", action='store_true', help="just autosomes")
    parser.add_argument("--no_cov", action='store_true', help="no coverage")
    parser.add_argument("--random_all", action='store_true', help="randomize all locations")

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    if not args.cfg or not os.path.exists(args.cfg):
        sys.exit('Need to provide config file to extract chromosome sizes and gap coordinates')
   
    locations = select_locations(args.m, args.Dmin, args.Dmax, args.s, args.w, args.sw, args.genome, args.cfg, args.random_all,
                                 no_segdup=args.no_segdup, x_only=args.x_only, auto_only=args.auto_only)

    if args.nprocs <= 1:
        distances = get_distances(args.Dmax, args.Dmin, args.s)
        bam = pysam.AlignmentFile(args.bam, 'rb')
        good_barcodes = None
        if args.good_barcodes is not None:
            good_barcodes = get_good_barcodes(args.good_barcodes)
        data = generate_data(bam, locations, args.w, args.sw, good_barcodes, args.hap, args.no_cov, args.random_all, args.tech)
    else:
        data = multi_generate_data(args.bam, locations, args.w, args.sw, args.good_barcodes, args.hap, args.no_cov, args.random_all, args.tech, args.nprocs)
    
    params = 'm={} Dmin={} Dmax={} s={} w={} bam={} no_segdup={} sw={} tech={} good_barcodes={} hap={}'.format(args.m,
                                                                                                               args.Dmin,
                                                                                                               args.Dmax,
                                                                                                               args.s,
                                                                                                               args.w,
                                                                                                               args.bam,
                                                                                                               args.no_segdup,
                                                                                                               args.sw,
                                                                                                               args.tech,
                                                                                                               args.good_barcodes,
                                                                                                               args.hap)
    output(data, args.out, params, args.no_cov)

if __name__ == '__main__':
    main()
