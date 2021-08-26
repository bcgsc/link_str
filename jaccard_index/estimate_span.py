import pandas as pd
import numpy as np
import argparse
from model_span import get_coverage, get_flanks, get_all_barcode_bounds, extract_barcodes, calc_jaccard, get_good_barcodes
import pysam
import sys
from glob import glob
import os

def get_coverages(bam, chrom, pos, flank_size, hap):
    cov_left = get_coverage(bam, chrom, pos[0] - flank_size, pos[0] - 1, hap)
    cov_right = get_coverage(bam, chrom, pos[1] + 1, pos[1] + flank_size, hap)

    if cov_left is not None and cov_right is not None:
        return round(cov_left, 2), round(cov_right, 2)
    elif cov_left is None and cov_right is None:
        return np.nan, np.nan
    elif cov_left is None:
        return np.nan, round(cov_right, 2)
    else:
        return round(cov_left, 2), np.nan

def extract_coverages(df, bam, flank_size, hap):
    df['covA'], df['covB'] = zip(*df.apply(lambda c: get_coverages(bam, c['chr'], (c['posA_exp'], c['posB_exp']), flank_size, hap), axis=1,))

def apply_min_d(min_d, d_ref, posA, posB):
    if d_ref < min_d:
        pad = np.ceil(float(min_d - d_ref) / 2)
        return posA - pad, posB + pad, posB + pad - (posA - pad) - 1
    else:
        return posA, posB, d_ref

def apply_min_ds(df, min_d):
    df['posA_exp'], df['posB_exp'], df['d_exp'] = zip(*df.apply(lambda c: apply_min_d(min_d,
                                                                                      c['d_ref'],
                                                                                      c['posA_ref'],
                                                                                      c['posB_ref'],
                                                                                      ),
                                                                axis=1
                                                   ))   

def parse_bed(bed_file):
    #cols = ['chr', 'posA_ref', 'posB_ref']
    df = pd.read_csv(bed_file, header=None, sep='\t')
    # 4th column (optional) is gene name, 5th column if provided should be motif, not used
    if df.shape[1] > 4:
        df = df.iloc[:,:4]
    cols = ['chr', 'posA_ref', 'posB_ref']
    if df.shape[1] == 3:
        df.columns = cols
        df['gene'] = np.nan
    elif df.shape[1] == 4:
        cols.append('gene')
        df.columns = cols
 
    df['chr'] = df['chr'].astype(str)
    df['d_ref'] = df['posB_ref'] - df['posA_ref'] - 1

    return df

def get_params(model_file):
    params = {}
    with open(model_file, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                params = dict([p.split('=') for p in line[1:].lstrip().rstrip().split()])
            break
    return params

def compute_jaccard(bam, pos1, pos2, w, sw, good_barcodes, hap, tech):
    flanks = get_flanks(pos1[1], pos2[1], w)
    chrom = str(pos1[0])
    
    barcode_bounds = get_all_barcode_bounds(bam, chrom, [flanks[0][0], flanks[1][1]],
                                            sw=sw, 
                                            hap=hap, 
                                            good_barcodes=good_barcodes,
                                            tech=tech)
    
    bc1 = extract_barcodes(flanks[0], barcode_bounds)
    bc2 = extract_barcodes(flanks[1], barcode_bounds)

    result = calc_jaccard(bc1, bc2)
    if not bc1 or not bc2:
        print('pp', pos1, pos2, bc1, bc2, result)
    return calc_jaccard(bc1, bc2)

def get_jaccard_indices(df, bam, w, sw, good_barcodes, hap, tech):
    df['JAB'], df['NA'], df['NB'], df['NAB'] = zip(*df.apply(lambda c: compute_jaccard(bam,
                                                                                       (c['chr'], c['posA_exp']),
                                                                                       (c['chr'], c['posB_exp']),
                                                                                       w,
                                                                                       sw,
                                                                                       good_barcodes,
                                                                                       hap,
                                                                                       tech,
                                                                                       ), 
                                                             axis=1,
                                                             ))
    

def get_estimates(df_coords, df_model, E, F, C, ignore_cov):
    df_coords['median'], df_coords['Q1'], df_coords['Q3'], df_coords['IQR'] = zip(*df_coords.apply(lambda c: get_estimate((c['JAB'], c['NA'], c['NB'], c['NAB'], c['covA'], c['covB']),
                                                                                                                          df_model,
                                                                                                                          E,
                                                                                                                          F,
                                                                                                                          C,
                                                                                                                          ignore_cov),
                                                                                                   axis=1,
                                                                                  ))
   
def find_closest(j_tuple, df_model, E, F, C, ignore_cov):
    df_closest = df_model.iloc[(df_model['JAB']-j_tuple[0]).abs().argsort()[:E]].copy()
    df_closest['f'] = df_closest.apply(lambda c: calculate_distance((c['NA'], c['NB'], c['NAB']),
                                                                    j_tuple[1:]),
                                       axis=1,
                                       )
    if 'covA' in df_model.columns and not ignore_cov:
        df_closest['c'] = df_closest.apply(lambda c: calculate_cov_diff(j_tuple[4:], (c['covA'], c['covB'])),
                                           axis=1,)
        df_filtered = df_closest[(df_closest['f'] < F) & (df_closest['c'] <= C)]
    else:
        df_filtered = df_closest[df_closest['f'] < F]

    print(j_tuple)
    print(df_filtered)
    #df_filtered.to_csv('tmp.tsv', sep='\t')

    return df_filtered
    
def get_estimate(j_tuple, df_model, E, F, C, ignore_cov):
    df = find_closest(j_tuple, df_model, E, F, C, ignore_cov)
    q1 =  df['d'].quantile(0.25)
    q3 = df['d'].quantile(0.75)
    return df['d'].median(), q1, q3, q3 - q1 

def calculate_distance(tuple1, tuple2):
    d1 = np.abs(tuple1[0] - tuple2[0]) + np.abs(tuple1[1] - tuple2[1]) + np.abs(tuple1[2] - tuple2[2])
    d2 = np.abs(tuple1[0] - tuple2[1]) + np.abs(tuple1[1] - tuple2[0]) + np.abs(tuple1[2] - tuple2[2])
    return min(d1, d2)

def calculate_cov_diff(covs1, covs2):
    c1 = sorted(covs1)
    c2 = sorted(covs2)
    diff1 = np.abs(c1[0] - c2[0])
    diff2 = np.abs(c1[1] - c2[1])

    return max(diff1, diff2)

def pick_best(dfs):
    df_concat = pd.concat(dfs)
    
    idx = df_concat.groupby(['chr', 'posA_exp', 'posB_exp'])['IQR'].transform(min) == df_concat['IQR']
     
    return df_concat[idx]

def calculate_estimate(df):
    df['d_exp'] = df['posB_exp'] - df['posA_exp'] - 1
    df['estimate'] = df['median'] - (df['d_exp'] - df['d_ref'])
    df['bound_lower'] = df['Q1'] - (df['d_exp'] - df['d_ref'])
    df['bound_upper'] = df['Q3'] - (df['d_exp'] - df['d_ref'])

def output(df, out_file):
    cols = ['chr', 'posA_ref', 'posB_ref', 'd_ref', 'gene',
            'posA_exp', 'posB_exp', 'd_exp', 'covA', 'covB',
            'median', 'Q1', 'Q3', 'IQR',
            'estimate', 'bound_lower', 'bound_upper', 'model']
    integer_cols = ['posA_exp', 'posB_exp', 'd_ref', 'd_exp', 'median', 'Q1', 'Q3', 'IQR', 'estimate', 'bound_lower', 'bound_upper']
    df[integer_cols] = df[integer_cols].fillna(-1)
    df[integer_cols] = df[integer_cols].astype('int')
    df[integer_cols] = df[integer_cols].astype('str')
    df[integer_cols] = df[integer_cols].replace('-1', np.nan)
    df[cols].to_csv(out_file, sep='\t', na_rep='NA', index=False)
        
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("bed", type=str, help="bed file")
    parser.add_argument("models", type=str, nargs='+', help="model file(s)")
    parser.add_argument("out", type=str, help="output file")
    parser.add_argument("--hap", type=int, help="haplotype 1 or 2")
    parser.add_argument("-E", type=int, default=200, help="size of nearest Jaccard index neighborhood. Default:200")
    parser.add_argument("-F", type=int, default=20, help="maximum search distance. Default:20")
    parser.add_argument("-C", type=float, default=1, help="maximum search distance. Default:1")
    parser.add_argument("--ignore_cov", action='store_true', help="ignore coverage matching")
    parser.add_argument("--min_d", type=int, help="minimum d")

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam, 'rb')
    
    # read coordinates
    df_bed = parse_bed(args.bed)

    models = []
    if len(args.models) == 1 and os.path.isdir(args.models[0]):
        models = sorted(glob('%s/*.mod' % args.models[0]))
    else:
        models = args.models

    df_results = []
    for model in models:
        print(model)
        # make clean df_bed
        df = df_bed.copy(deep=True)
        df['model'] = model
        
        # get parameters from model
        params = get_params(model)
        
        # read model
        df_model = pd.read_csv(model, skiprows=[0], sep='\t')

        # apply minimum distance based on model
        min_d = df_model['d'].min()
        if args.min_d and args.min_d > min_d:
            min_d = args.min_d
        apply_min_ds(df, min_d)

        extract_coverages(df, bam, int(params['w']), args.hap)

        good_barcodes = None
        if 'good_barcodes' in params and os.path.exists(params['good_barcodes']):
            good_barcodes = get_good_barcodes(params['good_barcodes'])

        # calcuate jaccard indices for each coordinate
        get_jaccard_indices(df, bam,
                            int(params['w']),
                            int(params['sw']),
                            good_barcodes,
                            args.hap,
                            params['tech'])

        df.dropna(subset=['JAB', 'NA', 'NB', 'NAB', 'covA', 'covB'], inplace=True)

        # collect estimates
        get_estimates(df, df_model, args.E, args.F, args.C, args.ignore_cov)

        df_results.append(df)

    # pick best estimate based on smallest IQR
    df_best = pick_best(df_results)
    
    # adjust estimate to real coordinate
    calculate_estimate(df_best)
    
    # include coordinates that don't have any results
    df_final = pd.merge(df_best, df_bed, how='right', on=['chr', 'posA_ref', 'posB_ref', 'd_ref', 'gene'])

    # output in tsv format
    output(df_final, args.out)

main()
