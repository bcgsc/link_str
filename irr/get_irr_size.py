import sys
from subprocess import check_output, Popen, PIPE
import re
import numpy as np
import argparse
import os

def get_irr_calls(infile, hp=None):
    calls = {}
    #print(hp)
    with open(infile, 'r') as ff:
        for line in ff:
            if hp is None and 'total:' in line:
                cols = line.rstrip().split()
                coord = cols[0]
                pairs = int(cols[-2])
                anchors = int(cols[-1])
                calls[coord] = {'pairs': pairs, 'anchors':anchors}
            elif hp is not None and 'hp:' in line:
                cols = line.rstrip().split(' ')
                coord = cols[0]
                if cols[1].split(':')[1] == str(hp):
                    pairs = int(cols[-2])
                    anchors = int(cols[-1])

                    if not coord in calls:
                        calls[coord] = {'pairs': pairs, 'anchors':anchors}
                    else:
                        calls[coord]['pairs'] += pairs
                        calls[coord]['anchors'] += anchors

                elif cols[1].split(':')[1] == 'na':
                    pairs = int(cols[-2])/2
                    anchors = int(cols[-1])/2

                    if not coord in calls:
                        calls[coord] = {'pairs': pairs, 'anchors':anchors}
                    else:
                        calls[coord]['pairs'] += pairs
                        calls[coord]['anchors'] += anchors

    for coord in calls.keys():
        calls[coord]['total'] = int(2*calls[coord]['pairs'] + calls[coord]['anchors'])

    return calls
    #return int(2*pairs + anchors)

def get_size(irrs, cov, rlen, repeat_len=None):
    sizes = []
    if type(rlen) is int:
        sizes.append(rlen * (1 + irrs/cov))
    else:
        for r in rlen:
            sizes.append(r * (1 + irrs/cov))

    if repeat_len is not None:
        sizes = [s/repeat_len for s in sizes]
 
    if len(sizes) == 1:
        return round(sizes[0])
    else:
        return '{},{}'.format(round(sizes[0]), round(sizes[1]))

def get_cov(bam, chrom, start, end, hp):
    cov = None
    scripts_dir = os.path.dirname(os. path. realpath(__file__))
    cmd = 'python {}/get_depth.py {} {} {} {}'.format(scripts_dir, bam, chrom, start, end)
    if hp is not None:
        cmd += ' --hp {}'.format(hp)
    with Popen(cmd, stdout=PIPE, stderr=None, shell=True) as process:
        output = process.communicate()[0].decode("utf-8")
        #print(output)
        for line in output.split('\n'):
            if 'depth' in line:
                cov = float(line.split(' ')[1])

    return cov

rlens = {'10x': (128,151),
         'stlfr': 100,
         'tell_seq': 146}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("irr_out", type=str, help="id_irr.py output")
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("tech", type=str, help="tech")
    parser.add_argument("--hp", type=str, help="hp")
    parser.add_argument("--het", action='store_true', help="het, reduce coverage by half")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    irr_calls = get_irr_calls(args.irr_out, hp=args.hp)
    for locus in sorted(irr_calls.keys()):
        #print(locus, irr_calls[locus])

        chrom, start, end = re.split('[:-]', locus)
        cov = get_cov(args.bam, chrom, start, end, args.hp)
        if args.het:
            cov = cov/2
        #print('cov', cov)

        size = get_size(irr_calls[locus]['total'], cov, rlens[args.tech])
        print(locus, irr_calls[locus]['pairs'], irr_calls[locus]['anchors'], irr_calls[locus]['total'], cov, rlens[args.tech], size)
main()

