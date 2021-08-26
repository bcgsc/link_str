import sys
import subprocess
import os
from collections import defaultdict
import pysam
import argparse

def is_same_repeat(reps):
    if len(reps[0]) <= len(reps[1]):
        rep1, rep2 = reps[0], reps[1]
    else:
        rep2, rep1 = reps[0], reps[1]

    if rep1 == rep2:
        return True

    perms1 = []
    for i in range(len(rep1)):
        pat = rep1[i:] + rep1[:i]
        perms1.append(pat)

    for p1 in perms1:
        if p1 == rep2:
            return True

    return False

def type_trf_cols(cols):
    return list(map(int, cols[:3])) + [float(cols[3])] + list(map(int, cols[4:12])) + [float(cols[12])] + cols[13:]

def parse_trf(trf_output):
    results = defaultdict(list)
    with open(trf_output, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split()
            if not cols:
                continue
            if cols[0] == 'Sequence:':
                seq = cols[1]
            elif len(cols) == 15:
                results[seq].append(type_trf_cols(cols))

    return results

def screen_trf(all_results, repeat, seqs, edge=10, min_repeat_frac=0.8):
    repeat_rc = reverse_complement(repeat)

    irrs = {}
    bcs = {}
    for seq, results in all_results.items():
        for r in results:
            matched_repeat = None

            if is_same_repeat((r[13], repeat)):
                matched_repeat = repeat
            elif is_same_repeat((r[13], repeat_rc)):
                matched_repeat = repeat_rc

            if matched_repeat:
                rtype = None
                read, mate, bc, hap, rlen = seq.split('_')
                bcs[read] = bc
                if r[0] <= edge and r[1] >= int(rlen) - edge:
                    rtype = ('full', matched_repeat)
                elif r[1] - r[0] + 1 >= min_repeat_frac * int(rlen):
                    rtype = ('full', matched_repeat)
                elif r[0] <= edge:
                    rtype = ('part_start', matched_repeat)
                elif r[1] >= int(rlen) - edge:
                    rtype = ('part_end', matched_repeat)

                if rtype is not None:
                    if not read in irrs:
                        irrs[read] = {1:None, 2:None}

                    irrs[read][int(mate)] = (rtype, r)

    pairs = []
    anchors = {}
    for read in sorted(irrs.keys()):
        if irrs[read][1] and irrs[read][2]:
            if irrs[read][1][0][0] == 'full' and irrs[read][2][0][0] == 'full':
                if irrs[read][1][0][1] != irrs[read][2][0][1]:
                    #print('ww', read, bcs[read], 'irr_pair')
                    pairs.append(read)

            # one full and one partial - anchor
            elif irrs[read][1][0][0] == 'full' or irrs[read][2][0][0] == 'full':
                if irrs[read][1][0][0] == 'full':
                    m_full, m_part = 1,2
                else:
                    m_full, m_part = 2,1
                    
                # check if motifs are different
                if irrs[read][m_full][0][1] != irrs[read][m_part][0][1]:
                    anchors[read] = {}
                    anchors[read][m_full] = 'full'
                    if irrs[read][m_part][0][0] == 'part_start':
                        seq = seqs[(read, m_part)][int(irrs[read][m_part][1][1]):]
                    else:
                        seq = seqs[(read, m_part)][:int(irrs[read][m_part][1][0])-1]
                    # allow for sequencing errors
                    if set(seq).issubset(set(irrs[read][m_part][1][13])):
                        pairs.append(read)
                    else:
                        anchors[read][m_part] = seq
                        print('nn', read, irrs[read][m_part], anchors[read][m_part])

        elif irrs[read][1] and irrs[read][1][0][0] == 'full':
            seq = seqs[(read, 2)]
            anchors[read] = {1:'full', 2:seq}

        elif irrs[read][2] and irrs[read][2][0][0] == 'full':
            seq = seqs[(read, 1)]
            anchors[read] = {2:'full', 1:seq}

    return pairs, anchors

def make_trf_fasta(in_file, out_file, tech):
    # for checking anchor sequences

    if tech == '10x':
        trim_len = 23
    else:
        trim_len = 0
    seqs = {}
    barcodes = {}
    haps = {}
    with open(out_file, 'w') as out:
        with open(in_file, 'r') as ff:
            for line in ff:
                '''
                if line[:3] != 'ggg':
                    continue
                cols = line.split()[1:]
                '''
                cols = line.split()
                bc, hap, read = cols[:3]
                # these 2 lines for stlfr
                if '_' in read:
                    read = read.replace('_', ':')
                if '_' in bc:
                    bc = bc.replace('_', ':')

                barcodes[read] = bc
                haps[read] = hap
                seq1, seq2 = cols[-2:]
                if len(seq1) == len(seq2):
                    seq1 = seq1[trim_len:]

                h1 = '{}_1_{}_{}_{}'.format(read, bc, hap, len(seq1))
                h2 = '{}_2_{}_{}_{}'.format(read, bc, hap, len(seq2))

                out.write('>{}\n{}\n'.format(h1, seq1))
                out.write('>{}\n{}\n'.format(h2, seq2))
                seqs[(read, 1)] = seq1
                seqs[(read, 2)] = seq2

    return seqs, barcodes, haps

def make_blastn_query(anchors, out_file):
    with open(out_file, 'w') as out:
        for read in sorted(anchors.keys()):
            for mate in anchors[read]:
                if anchors[read][mate] != 'full':
                    out.write('>{}_{}\n{}\n'.format(read, len(anchors[read][mate]), anchors[read][mate]))

def make_blastn_subject(chrom, start, end, genome_fasta, out_file, size=500):
    fa = pysam.Fastafile(genome_fasta)
    with open(out_file, 'w') as out:
        out.write('>upper\n{}\n'.format(fa.fetch(chrom, int(start) - size, int(start))))
        out.write('>lower\n{}\n'.format(fa.fetch(chrom, int(end), int(end) + size)))

def parse_blastn(blastn_results):
    results = {}
    with open(blastn_results, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            if not cols[0] in results or\
               float(cols[-1]) > float(results[cols[0]][-1]):
                results[cols[0]] = cols[1:]
            
    return results

def run_blastn(query, subject, out_file):
    cmd = 'blastn -query {} -subject {} -outfmt 6 -word_size 10 > {}'.format(query, subject, out_file)
    FNULL = open(os.devnull, 'w')
    returncode = subprocess.call(cmd, shell=True)
    
    if os.path.exists(out_file):
        return parse_blastn(out_file)
    else:
        return None

def run_trf(fasta):
    trf_params = '2 5 5 80 10 10 500'
    cmd = 'trf {} {} -d -h'.format(fasta, trf_params)
    out_file = '{}.{}.dat'.format(fasta, trf_params.replace(' ', '.'))

    FNULL = open(os.devnull, 'w')
    returncode = subprocess.call(cmd, shell=True, stdout=FNULL, stderr=FNULL)

    if os.path.exists(out_file):
        return parse_trf(out_file), out_file
    else:
        return None, out_file

def screen_anchors(blastn_results, subject_size, min_mapped=0.8, edge=10):
    good = []
    for read_rlen, cols in blastn_results.items():
        read, rlen = read_rlen.split('_')
        rlen = int(rlen)
        sstart, send = int(cols[7]), int(cols[8])
        mapped = float(cols[2]) / rlen
        if mapped >= min_mapped:
            print('anchor_passed', read, mapped, cols[2], rlen, sstart, send)
            good.append(read)
        else:
            print('anchor_failed', read, mapped, cols[2], rlen, sstart, send)

    return good

def report(pairs, anchors, barcodes, haps, out_file=None):
    lines = []
    lines.append('total:{} {}'.format(len(pairs), len(anchors)))
    
    count_haps = {}
    for read in pairs:
        print('irr_pair', read, barcodes[read], haps[read])
        hp = haps[read]
        if not hp in count_haps:
            count_haps[hp] = {'pair':0, 'anchor':0}
        count_haps[hp]['pair'] += 1
    for read in anchors:
        hp = haps[read]
        if not hp in count_haps:
            count_haps[hp] = {'pair':0, 'anchor':0}
        count_haps[hp]['anchor'] += 1

    count_barcodes = {}
    for read in pairs:
        bc = barcodes[read]
        if not bc in count_barcodes:
            count_barcodes[bc] = {'pair':0, 'anchor':0}
        count_barcodes[bc]['pair'] += 1
    for read in anchors:
        bc = barcodes[read]
        if not bc in count_barcodes:
            count_barcodes[bc] = {'pair':0, 'anchor':0}
        count_barcodes[bc]['anchor'] += 1

    for bc in sorted(count_barcodes.keys()):
        lines.append('bc {}:{} {}'.format(bc, count_barcodes[bc]['pair'], count_barcodes[bc]['anchor']))
    for hp in sorted(count_haps.keys()):
        lines.append('hp {}:{} {}'.format(hp, count_haps[hp]['pair'], count_haps[hp]['anchor']))

    if out_file is None:
        for line in lines:
            print(line)
    else:
        with open(out_file, 'w') as out:
            for line in lines:
                out.write('{}\n'.format(line))

def reverse_complement(seq):
    """Reverse complements sequence string"""
    complement = str.maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("extract_output", type=str, help="extract_out")
    parser.add_argument("coord", type=str, nargs=3, help="STR locus coordinate: chromosome start end")
    parser.add_argument("motif", type=str, help="STR motif")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("--report", type=str, help="report file")
    parser.add_argument("--trf_out", type=str, help="trf output")
    parser.add_argument("--tech", type=str, default='10x', help="10x(default), tell_seq, or stlfr")
    parser.add_argument("--debug", action='store_true', help="debug mode i.e. keep trf output")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    chrom, start, end = args.coord[0], int(args.coord[1]), int(args.coord[2])
    flank_size = 500

    tmp_files = []
    out_fa = 'reads.fa'
    tmp_files.append(out_fa)
    seqs, barcodes, haps = make_trf_fasta(args.extract_output, out_fa, args.tech)

    if not args.trf_out:
        trf_results, trf_out = run_trf(out_fa)
        tmp_files.append(trf_out)
    else:
        trf_results = parse_trf(args.trf_out)
        tmp_files.append(args.trf_out)
    if trf_results:
        pairs, anchors = screen_trf(trf_results, args.motif.upper(), seqs)
    
        if anchors:
            query_fa = 'anchors.fa'
            subject_fa = 'flanks.fa'
            blastn_out = 'anchors.blastn'
            tmp_files.extend([query_fa, subject_fa, blastn_out])
            make_blastn_query(anchors, query_fa)
            make_blastn_subject(chrom, start, end, args.genome_fasta, subject_fa, size=flank_size)
            results = run_blastn(query_fa, subject_fa, blastn_out)

            anchors_passed = screen_anchors(results, flank_size)

            report(pairs, anchors_passed, barcodes, haps, out_file=args.report)

    # cleanup
    if not args.debug:
        for ff in tmp_files:
            if os.path.exists(ff):
                os.remove(ff)
main()   
