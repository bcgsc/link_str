import pysam
import sys
from collections import defaultdict
import gzip
import argparse

def extract_reads_10x(fqs, barcode2coords, barcode2hps, out_file):
    collected = defaultdict(list)
    barcode_len = len(list(barcode2coords.keys())[0])
    read_name = None
    with open(out_file, 'w') as out:
        with gzip.open(fqs[0], 'rt') as fq1, gzip.open(fqs[1], 'rt') as fq2:
            for line1, line2 in zip(fq1, fq2):
                if line1[0] == '@':
                    read_name = line1.lstrip('@').rstrip().split(' ')[0].split('/')[0]
                elif line1[0] == '+':
                    read_name = None
                elif read_name is not None:
                    barcode = line1[:barcode_len]
                    seq1 = line1.rstrip()
                    seq2 = line2.rstrip()
                    if barcode in barcode2coords:
                        for coord in barcode2coords[barcode]:
                            out.write('{} {} {} {} {} {}\n'.format(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2))
                            #print(coord, barcode, hps[barcode], read_name, seq1, seq2)
                        collected[barcode].append((read_name, seq1, seq2))

def extract_reads_10x_interleaved(fq_gz, barcode2coords, barcode2hps, out_file, targets=None):
    collected = defaultdict(list)
    barcode_len = len(list(barcode2coords.keys())[0])
    read_name = None
    count = 1
    with open(out_file, 'w') as out:
        with gzip.open(fq_gz, 'rt') as fq:
            for line in fq:
                if count % 4 == 1:
                    read_name = line.lstrip('@').rstrip().split(' ')[0].split('/')[0]
                    mate = 1 if count==1 or mate == 2 else 2
                elif count % 4 == 2:
                    if mate == 1:
                        seq1 = line.rstrip()
                        barcode = line[:barcode_len]
                    else:
                        seq2 = line.rstrip()

                        #print('zz', read_name, barcode, seq1, seq2)
                        if barcode in barcode2coords:
                            if targets is not None and read_name not in targets:
                                continue
                            for coord in barcode2coords[barcode]:
                                out.write('{} {} {} {} {} {}\n'.format(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2))
                                #print(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2)
                            collected[barcode].append((read_name, seq1, seq2))
                
                count += 1
    
def extract_reads_tell_seq(fqs, barcode2coords, barcode2hps, out_file):
    # 3 fqs, last one barcode
    collected = defaultdict(list)
    read_name = None
    with open(out_file, 'w') as out:
        with gzip.open(fqs[0], 'rt') as fq1, gzip.open(fqs[1], 'rt') as fq2, gzip.open(fqs[2], 'rt') as fq3:
            for line1, line2, line3 in zip(fq1, fq2, fq3):
                if line1[0] == '@':
                    read_name = line1.lstrip('@').rstrip().split(' ')[1]
                elif line1[0] == '+':
                    read_name = None
                elif read_name is not None:
                    barcode = line3.rstrip()
                    seq1 = line1.rstrip()
                    seq2 = line2.rstrip()
                    #print('zz', read_name, barcode, seq1, seq2)
                    if barcode in barcode2coords:
                        for coord in barcode2coords[barcode]:
                            out.write('{} {} {} {} {} {}\n'.format(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2))
                            #print(coord, barcode, read_name, barcode2hps[coord][barcode], seq1, seq2)
                        collected[barcode].append((read_name, seq1, seq2)) 

def extract_reads_stlfr(fqs, barcode2coords, barcode2hps, out_file):
    collected = defaultdict(list)
    read_name = None
    with open(out_file, 'w') as out:
        with gzip.open(fqs[0], 'rt') as fq1, gzip.open(fqs[1], 'rt') as fq2:
            for line1, line2 in zip(fq1, fq2):
                if line1[0] == '@' and len(line1.rstrip().split()) > 1:
                    #print('aa', line1.rstrip(), len(line1.rstrip().split()))
                    read_name, barcode = line1.lstrip('@').rstrip().split()[0].split('/')[0].split('#')
                    #print('bb', read_name, barcode)
                elif line1[0] == '+':
                    read_name = None
                elif read_name is not None:
                    seq1 = line1.rstrip()
                    seq2 = line2.rstrip()
                    if barcode in barcode2coords:
                        for coord in barcode2coords[barcode]:
                            out.write('{} {} {} {} {} {}\n'.format(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2))
                            #print(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2)
                        collected[barcode].append((read_name, seq1, seq2))

def extract_reads_stlfr_unzipped(fqs, barcode2coords, barcode2hps, out_file):
    collected = defaultdict(list)
    read_name = None
    with open(out_file, 'w') as out:
        with open(fqs[0], 'r') as fq1, open(fqs[1], 'r') as fq2:
            for line1, line2 in zip(fq1, fq2):
                if line1[0] == '@' and len(line1.rstrip().split()) > 1:
                    #print('aa', line1.rstrip(), len(line1.rstrip().split()))
                    read_name, barcode = line1.lstrip('@').rstrip().split()[0].split('/')[0].split('#')
                    #print('bb', read_name, barcode)
                elif line1[0] == '+':
                    read_name = None
                elif read_name is not None:
                    seq1 = line1.rstrip()
                    seq2 = line2.rstrip()
                    if barcode in barcode2coords:
                        for coord in barcode2coords[barcode]:
                            out.write('{} {} {} {} {} {}\n'.format(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2))
                            #print(coord, barcode, barcode2hps[coord][barcode], read_name, seq1, seq2)
                        collected[barcode].append((read_name, seq1, seq2))

def extract_aln_bc(aln, tech):
    bc = None
    if tech == '10x':
        try:
            bc = aln.get_tag('BX').split('-')[0]
        except:
            bc = None

    elif tech == 'tell_seq':
        try:
            bc = aln.get_tag('BX').split('-')[0]
        except:
            bc = None

    elif tech == 'stlfr':
        try:
            bc = aln.query_name.split('#')[1]
            if bc == '0_0_0':
                bc = None
        except:
            bc = None

    return bc

def get_barcodes(bam, chrom, start, end, w, tech, min_len=0, min_frags=2, close=10000, targets=None):
    hps = {}
    barcodes = defaultdict(list)
    for aln in bam.fetch(chrom, start-w/2, end+w/2):
        if targets is not None and not aln.query_name in targets:
            continue

        barcode = extract_aln_bc(aln, tech=tech)
        if barcode is None:
            continue

        hp = 'na'
        if aln.has_tag('HP'):
            hp = aln.get_tag('HP')

        if aln.is_proper_pair and aln.query_alignment_length == aln.infer_read_length():
            barcodes[barcode].extend([aln.reference_start, aln.reference_end])
            hps[barcode] = hp

    kept = []
    for bc in barcodes:
        coords = sorted(barcodes[bc])
        if coords[-1] - coords[0] < min_len:
            continue
        if len(barcodes[bc])/2 >= min_frags:
            if (coords[0] <= start and coords[-1] >= end) or\
               (coords[-1] < start and coords[-1] > start - close) or\
               (coords[0] > end and coords[0] < end + close) or\
               (coords[0] >= start and coords[0] <= end) or\
               (coords[-1] >= start and coords[-1] <= end):
                #print('bc', start, end, bc, hps[bc], len(barcodes[bc]), barcodes[bc])
                kept.append(bc)

    return kept, hps

def output_barcodes(bc2_coords, out_file):
    with open(out_file, 'w') as out:
        for bc, coords in bc2_coords.items():
            for coord in coords:
                out.write('{} {}\n'.format(coord, bc))

def get_targets(infile):
    targets = set()
    with open(infile, 'r') as ff:
        for line in ff:
            targets.add(line.rstrip())

    return targets

def get_coords(bed_file):
    coords = []
    with open(bed_file, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            coords.append((cols[0], int(cols[1]), int(cols[2])))

    return coords

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("coords", type=str, help="coords bed file")
    parser.add_argument("out", type=str, help="output file")
    parser.add_argument("--fqs", type=str, nargs='+', help="fastqs")
    parser.add_argument("--out_bc", type=str, help="output of list of barcodes")
    parser.add_argument("--w", type=int, default=500000, help="window size on each side")
    parser.add_argument("--tech", type=str, default='10x', help="10x(default), tell_seq, or stlfr")
    parser.add_argument("--unzipped", action='store_true', help="fastq files unzipped (only for stlfr")
    parser.add_argument("--targets", type=str, help="target reads")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bam = pysam.Samfile(args.bam, 'rb')
    
    coords = get_coords(args.coords)

    targets = None
    if args.targets:
        targets = get_targets(args.targets)

    barcode2coords = defaultdict(list)
    barcode2hps = defaultdict(dict)
    for chrom, start, end in coords:
        coord = '{}:{}-{}'.format(chrom, start, end)
        barcodes, hps = get_barcodes(bam, chrom, start, end, args.w, args.tech, targets=targets)
        #print('qq', chrom, start, end, barcodes, hps)

        for bc in barcodes:
            barcode2coords[bc].append(coord)
            
        barcode2hps[coord] = hps

    if args.out_bc:
        output_barcodes(barcode2coords, args.out_bc)

    if args.fqs:
        if len(args.fqs) == 2 and args.tech == '10x':
            extract_reads_10x(args.fqs, barcode2coords, barcode2hps, args.out)
        elif len(args.fqs) == 3 and args.tech == 'tell_seq':
            extract_reads_tell_seq(args.fqs, barcode2coords, barcode2hps, args.out)
        elif len(args.fqs) == 2 and args.tech == 'stlfr':
            if args.unzipped:
                extract_reads_stlfr_unzipped(args.fqs, barcode2coords, barcode2hps, args.out)
            else:
                extract_reads_stlfr(args.fqs, barcode2coords, barcode2hps, args.out)
        elif len(args.fqs) == 1 and args.tech == '10x':
            extract_reads_10x_interleaved(args.fqs[0], barcode2coords, barcode2hps, args.out, targets=targets)

main()
