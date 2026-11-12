#!/usr/bin/env python3
"""500 kb TE hotspot enrichment for gene families, SD and recurrent CNV.

Inputs:
  --te       TSV: chr start end <one or more TE count columns>
  --genes    BED: chr start end gene_id (one row per gene locus)
  --families TSV: gene_id family (one gene can have multiple family rows)
  --sd       BED: chr start end [other columns]
  --cnv      BED: chr start end [other columns]; recurrent CNV REGION only

All files MUST use the same genome assembly and 0-based, half-open coordinates.
Dependencies: pandas, numpy, scipy, statsmodels. No bedtools required.
"""
import argparse
from collections import defaultdict
from math import ceil

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def bed3(path, extra=0):
    x = pd.read_csv(path, sep='\t', header=None, comment='#', dtype={0: str})
    if x.shape[1] < 3 + extra:
        raise ValueError(f'{path}: expected at least {3 + extra} BED columns')
    x = x.iloc[:, :3 + extra].copy()
    x.iloc[:, 1:3] = x.iloc[:, 1:3].apply(pd.to_numeric, errors='raise')
    x.columns = ['chr', 'start', 'end'] + (['gene_id'] if extra else [])
    x['start'] = x['start'].astype(np.int64)
    x['end'] = x['end'].astype(np.int64)
    return x[x['end'] > x['start']].reset_index(drop=True)


def gene_to_window(genes, windows):
    """Gene midpoint -> a single, non-overlapping genomic window."""
    out = []
    gene_groups = dict(tuple(genes.groupby('chr', sort=False)))
    for chrom, w in windows.groupby('chr', sort=False):
        g = gene_groups.get(chrom)
        if g is None or len(g) == 0:
            continue
        w = w.sort_values('start')
        starts = w['start'].to_numpy(dtype=np.int64)
        ends = w['end'].to_numpy(dtype=np.int64)
        ids = w['wid'].to_numpy(dtype=np.int64)
        mids = ((g['start'].to_numpy(dtype=np.int64) + g['end'].to_numpy(dtype=np.int64)) // 2)
        pos = np.searchsorted(starts, mids, side='right') - 1
        valid = pos >= 0
        valid[valid] &= mids[valid] < ends[pos[valid]]
        if valid.any():
            out.append(pd.DataFrame({'wid': ids[pos[valid]], 'gene_id': g['gene_id'].to_numpy()[valid]}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=['wid', 'gene_id'])


def interval_coverage(bed_path, windows):
    """Fraction of each window covered by the UNION of BED intervals."""
    ann = bed3(bed_path)
    out = np.zeros(len(windows), dtype=float)
    groups = dict(tuple(ann.groupby('chr', sort=False)))
    for chrom, w in windows.groupby('chr', sort=False):
        a = groups.get(chrom)
        if a is None or a.empty:
            continue
        # Merge intervals to avoid double counting overlapping annotations.
        intervals = a[['start', 'end']].sort_values(['start', 'end']).to_numpy(dtype=np.int64)
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1]:
                merged.append([int(s), int(e)])
            else:
                merged[-1][1] = max(merged[-1][1], int(e))
        w = w.sort_values('start')
        j = 0
        for row in w.itertuples(index=False):
            while j < len(merged) and merged[j][1] <= row.start:
                j += 1
            k, overlap = j, 0
            while k < len(merged) and merged[k][0] < row.end:
                overlap += max(0, min(row.end, merged[k][1]) - max(row.start, merged[k][0]))
                k += 1
            out[row.wid] = overlap / (row.end - row.start)
    return out


def hotspot_mask(counts, target_n):
    """Choose largest nonzero values, retaining ALL ties at cutoff."""
    counts = pd.to_numeric(counts, errors='raise').to_numpy(dtype=float)
    nonzero = counts[counts > 0]
    if nonzero.size == 0:
        return None, np.nan
    cutoff = np.sort(nonzero)[-min(target_n, nonzero.size)]
    return counts >= cutoff, float(cutoff)


def calc(hot, feature):
    hot = np.asarray(hot, dtype=bool)
    feature = np.asarray(feature, dtype=bool)
    a = int(np.count_nonzero(hot & feature))
    b = int(np.count_nonzero(hot & ~feature))
    c = int(np.count_nonzero(~hot & feature))
    d = int(np.count_nonzero(~hot & ~feature))
    odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    fe = (a / (a + b)) / ((a + c) / (a + b + c + d)) if (a + c) else np.nan
    # Finite plotting value, even when the raw OR is infinite.
    or_plot = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    return a, b, c, d, odds, fe, p, or_plot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--te', required=True)
    ap.add_argument('--genes', required=True)
    ap.add_argument('--families', required=True)
    ap.add_argument('--sd', required=True)
    ap.add_argument('--cnv', required=True)
    ap.add_argument('--top', type=float, default=0.05, help='Top fraction of TE bins (default 0.05)')
    ap.add_argument('--sd-cutoff', type=float, default=0.10, help='SD-covered fraction to call SD-rich')
    ap.add_argument('--cnv-cutoff', type=float, default=0.10, help='CNV-covered fraction to call CNV-associated')
    ap.add_argument('--min-bins', type=int, default=3, help='Minimum number of positive bins for a feature')
    ap.add_argument('--prefix', default='TE500kb')
    args = ap.parse_args()

    te = pd.read_csv(args.te, sep='\t', dtype={'chr': str})
    if not {'chr', 'start', 'end'}.issubset(te.columns):
        raise ValueError('TE file requires header columns: chr, start, end')
    te_cols = [c for c in te.columns if c not in ('chr', 'start', 'end')]
    if not te_cols:
        raise ValueError('TE file requires at least one TE-count column')
    te[['start', 'end']] = te[['start', 'end']].astype(np.int64)
    if (te.end <= te.start).any():
        raise ValueError('TE windows must have start < end')
    te = te.reset_index(drop=True)
    te['wid'] = np.arange(len(te))
    # Windows must be non-overlapping, otherwise gene midpoint assignment is ambiguous.
    for _, w in te.groupby('chr'):
        w = w.sort_values('start')
        if (w['start'].to_numpy()[1:] < w['end'].to_numpy()[:-1]).any():
            raise ValueError('Overlapping windows on the same chromosome')

    genes = bed3(args.genes, extra=1)
    genes['gene_id'] = genes['gene_id'].astype(str)
    family = pd.read_csv(args.families, sep='\t', dtype=str)
    if not {'gene_id', 'family'}.issubset(family.columns):
        raise ValueError('Family table requires columns: gene_id, family')
    family = family[['gene_id', 'family']].dropna().drop_duplicates()
    assigned = gene_to_window(genes, te)
    gf = assigned.merge(family, on='gene_id', how='inner').drop_duplicates(['wid', 'family', 'gene_id'])
    count_table = gf.groupby(['family', 'wid'], observed=True)['gene_id'].nunique().reset_index(name='count')
    features = {}
    for fam, sub in count_table.groupby('family', sort=False):
        cnt = np.zeros(len(te), dtype=np.int32)
        cnt[sub['wid'].to_numpy(dtype=int)] = sub['count'].to_numpy(dtype=int)
        features[('gene_family_present', fam)] = cnt >= 1
        features[('gene_family_multi_2plus', fam)] = cnt >= 10 #

    sd_frac = interval_coverage(args.sd, te)
    cnv_frac = interval_coverage(args.cnv, te)
    features[('genomic_region', 'SD_coverage')] = sd_frac >= args.sd_cutoff
    features[('genomic_region', 'recurrent_CNV_region')] = cnv_frac >= args.cnv_cutoff

    # Export per-window annotations, useful for debugging and plotting.
    ann = te[['chr', 'start', 'end', 'wid']].copy()
    ann['SD_fraction'] = sd_frac
    ann['recurrent_CNV_fraction'] = cnv_frac
    ann.to_csv(args.prefix + '_window_annotations.tsv', sep='\t', index=False)

    rows = []
    target_n = max(1, ceil(len(te) * args.top))
    for col in te_cols:
        hot, cutoff = hotspot_mask(te[col], target_n)
        if hot is None:
            continue
        for (ftype, fname), mark in features.items():
            if int(mark.sum()) < args.min_bins:
                continue
            a,b,c,d,odds,fe,p,or_plot = calc(hot, mark)
            rows.append(dict(TE_group=col, feature_type=ftype, feature=fname,
                             n_hot=int(hot.sum()), n_feature=int(mark.sum()),
                             TE_cutoff=cutoff, a=a, b=b, c=c, d=d,
                             OR=odds, FE=fe, pvalue=p, log2_OR_plot=np.log2(or_plot)))
    if not rows:
        raise ValueError('No valid tests; check inputs / --min-bins')
    res = pd.DataFrame(rows)
    res['FDR'] = multipletests(res.pvalue, method='fdr_bh')[1]
    res.sort_values(['FDR', 'pvalue', 'TE_group'], inplace=True)
    res.to_csv(args.prefix + '_enrichment.tsv', sep='\t', index=False)
    print(f'Windows: {len(te)}; TE groups: {len(te_cols)}; tests: {len(res)}')
    print(f'Gene-family genes mapped: {gf.gene_id.nunique()}; output: {args.prefix}_enrichment.tsv')
    print('2x2 table: a=TE-hot & feature+, b=TE-hot & feature-, c=not-hot & feature+, d=not-hot & feature-')


if __name__ == '__main__':
    main()
