#!/usr/bin/env python3
import argparse
import sys
import os
import pandas as pd

if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Process mpileup VCF and known SNPs to output novel variant candidates.'
    )
    ap.add_argument('--v', required=True, help='VCF file (can be .vcf or .vcf.gz)')
    ap.add_argument('--snp', required=False, help='BED file containing known SNPs')
    ap.add_argument('--o', required=True, help='Output CSV file')
    ap.add_argument('--min_variant', type=float, default=0,
                    help='Minimum %variant allele frequency to keep (default=0, no filtering)')
    args = ap.parse_args()

    try:
        v = pd.read_csv(
            args.v,
            comment='#',
            sep='\t',
            header=None,
            names=[
                'chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'tags', 'taginfo'
            ],
            low_memory=False
        )
    except Exception as e:
        sys.exit(f"❌ Failed to read VCF file: {e}")

    try:
        snp = pd.read_csv(
            args.snp,
            sep='\t',
            header=None,
            names=['chr', 'pos-1', 'pos', 'strand'],
            low_memory=False
        ).drop(columns=['pos-1'])
    except Exception as e:
        sys.exit(f"❌ Failed to read SNP BED file: {e}")

    # Standardize data types
    v['chr'] = v['chr'].astype(str)
    snp['chr'] = snp['chr'].astype(str)
    v['pos'] = pd.to_numeric(v['pos'], errors='coerce').fillna(0).astype(int)
    snp['pos'] = pd.to_numeric(snp['pos'], errors='coerce').fillna(0).astype(int)

    # Extract DP4 from INFO field
    v['DP4'] = v['info'].str.extract(r'DP4=([\d,]+)')[0]
    v[['ref_fwd', 'ref_rev', 'alt_fwd', 'alt_rev']] = v['DP4'].str.split(',', expand=True)

    for col in ['ref_fwd', 'ref_rev', 'alt_fwd', 'alt_rev']:
        v[col] = pd.to_numeric(v[col], errors='coerce').fillna(0).astype(int)

    # Calculate coverage and variant reads
    v['coverage'] = v['ref_fwd'] + v['ref_rev'] + v['alt_fwd'] + v['alt_rev']
    v['alt_reads'] = v['alt_fwd'] + v['alt_rev']

    # Avoid divide-by-zero
    v = v[v['coverage'] > 0]

    # Compute %variant
    v['%variant'] = round((v['alt_reads'] / v['coverage']) * 100, 2)

    # Apply minimum %variant filter here
    v = v[v['%variant'] >= args.min_variant]

    # Clean dataframe
    v_cut = v[['chr', 'pos', 'ref', 'alt', 'coverage', 'alt_reads', '%variant']]

    # Remove known SNPs
    merged = v_cut.merge(snp, how='left', on=['chr', 'pos'])
    merged_clean = merged[~merged['strand'].isin(['+', '-'])].copy()
    merged_clean['#variant_reads'] = merged_clean['alt_reads']
    merged_clean = merged_clean.drop(columns=['strand', 'alt_reads'])

    # Save result
    try:
        merged_clean.to_csv(args.o, index=False)
        print(f"✅ Output written to {args.o}")
    except Exception as e:
        sys.exit(f"❌ Failed to write output: {e}")
