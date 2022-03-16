import pysam
import numpy as np
import argparse

def get_depth(bam, chrom, start, end, target_hp=None, target_rg=None):
    depths = []
    for pc in bam.pileup(chrom, start, end):
        if target_hp is not None:
            depth = 0
            for pileupread in pc.pileups:
                if pileupread.alignment.has_tag('HP') and pileupread.alignment.get_tag('HP') == target_hp:
                   # print(pc.pos, pileupread.alignment.query_name, pileupread.alignment.get_tag('HP'), target_hp, pc.n)
                    if target_rg is not None:
                        if pileupread.alignment.has_tag('RG') and pileupread.alignment.get_tag('RG') == target_rg:
                            depth += 1
                    else:
                        depth += 1
        else:
            depth = pc.n

        depths.append(depth)

    print(depths)
    return np.median(depths)

def get_depth_coord(bam, coord, s=500, w=1000, hp=None, rg=None):
    depth_up = get_depth(bam, coord[0], coord[1]-s-w, coord[1]-s, target_hp=hp, target_rg=rg)
    depth_down = get_depth(bam, coord[0], coord[1]+s, coord[1]+s+w, target_hp=hp, target_rg=rg)

    print(depth_up, depth_down)
    if hp is not None and not depth_up and not depth_down:
        depth_up = get_depth(bam, coord[0], coord[1]-s-w, coord[1]-s, target_hp=None, target_rg=rg)
        depth_down = get_depth(bam, coord[0], coord[1]+s, coord[1]+s+w, target_hp=None, target_rg=rg)
        return np.mean((depth_up, depth_down))/2

    return np.mean((depth_up, depth_down))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("coord", type=str, nargs=3, help="chrom start end")
    parser.add_argument("--s", type=int, default=500, help="padding away from coord, default:500")
    parser.add_argument("--w", type=int, default=1000, help="span size to calculate coverage, default:1000")
    parser.add_argument("--hp", type=int, help="haplotype 1 or 2")
    parser.add_argument("--rg", type=str, help="read group")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam)
    coord = args.coord
    depth = get_depth_coord(bam, (coord[0], int(coord[1]), int(coord[2])), s=args.s, w=args.w, hp=args.hp, rg=args.rg)
    print('depth', depth)

main()
'''
bam = pysam.AlignmentFile('/projects/CIHR_LRPG/INDHU/linkedREADS/SIM/stLFR/hg38.ATXN10_4000.sort.mdup.rg.bam', 'rb')
coord = ('chr22', 45795355, 45795424)
get_depth_coord(bam, coord)
'''
