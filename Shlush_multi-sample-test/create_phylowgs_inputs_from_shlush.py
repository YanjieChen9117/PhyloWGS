#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Create PhyloWGS-style ssm_data.txt and cnv_data.txt from preprocessed
Shlush mutation and copy-number tables.

This script is intentionally simple and geared toward testing. It makes
reasonable assumptions to map the Shlush schema to the PhyloWGS formats:

  - Mutation table columns (per row):
        contig, position, norm_read, tum_read, tum_vaf, gene_name, sample, variant
    are mapped as:
        * chromosome = contig with optional 'chr' prefix stripped
        * position   = position (1-based)
        * a          = norm_read (approx. reference reads in tumour)
        * d          = norm_read + tum_read (approx. total depth)
        * gene       = "<chromosome>_<position>" (e.g. "1_234432")
        * mu_r,mu_v  = 0.999, 0.5 (as in other PhyloWGS examples)

  - Copy-number table columns (per row):
        contig, start, stop, cn, sample, type
    are mapped as:
        * chromosome = contig with optional 'chr' prefix stripped
        * start, stop = start, stop
        * total_cn     = cn (integer)
        * allele-specific CN (major_cn, minor_cn) are *approximated* from total_cn:
              total_cn == 0 -> major_cn=0, minor_cn=0
              total_cn == 1 -> major_cn=1, minor_cn=0
              total_cn == 2 -> major_cn=1, minor_cn=1
              total_cn >= 3 -> major_cn=total_cn-1, minor_cn=1
        * cell_prev is set to 1.0 for all segments (clonal assumption).

  - PhyloWGS output formats:
        ssm_data.txt:
            id  gene    a   d   mu_r    mu_v
            s0  1_2344  52  57  0.999   0.5
        cnv_data.txt:
            cnv a   d   ssms    physical_cnvs
            c0  0   266773091   s2,0,1;...    chrom=1,start=...,end=...,major_cn=1,minor_cn=0,cell_prev=1

    Here we generate *single-sample* inputs, so a and d are scalars (not
    comma-separated vectors). The multi-sample extension (comma-separated
    per-sample values) can be added later, following parser/test/outputs/
    multisamp_cnvs/{ssm_data,cnv_data}.txt.
"""

from __future__ import print_function

import csv
import os
import math


READ_LENGTH = 100.0  # used only for a rough CNV depth estimate


def norm_chr(contig):
    """Normalize chromosome names to match PhyloWGS expectations."""
    if contig is None:
        return ''
    contig = str(contig)
    if contig.lower().startswith('chr'):
        contig = contig[3:]
    # Map common non-numeric chromosomes if present
    mapping = {'x': '23', 'y': '24', 'mt': '25', 'm': '25'}
    key = contig.lower()
    return mapping.get(key, contig)


def read_shlush_mutation_table(path):
    """
    Read a Shlush mutation_table CSV and return a list of records with:
        {
          'chrom': <string>,
          'pos': <int>,
          'a': <int>,   # reference reads (approx)
          'd': <int>,   # total depth (approx)
        }
    """
    muts = []
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        required = ['contig', 'position', 'norm_read', 'tum_read']
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError("Missing required mutation columns in %s: %r" %
                             (path, missing))
        for row in reader:
            try:
                chrom = norm_chr(row['contig'])
                pos = int(row['position'])
                norm_read = int(row['norm_read'])
                tum_read = int(row['tum_read'])
            except (ValueError, TypeError):
                # Skip malformed rows
                continue

            a = max(0, norm_read)
            d = max(1, norm_read + tum_read)

            muts.append({
                'chrom': chrom,
                'pos': pos,
                'a': a,
                'd': d,
            })
    return muts


def decompose_total_cn(total_cn):
    """
    Heuristic decomposition of total copy number into (major_cn, minor_cn).
    This is necessarily approximate because the Shlush table does not
    provide allele-specific CN. The rule is:

        0 -> (0, 0)  (homozygous deletion)
        1 -> (1, 0)  (LOH)
        2 -> (1, 1)  (diploid)
        >=3 -> (total_cn - 1, 1)  (simple gain)
    """
    if total_cn <= 0:
        return 0, 0
    if total_cn == 1:
        return 1, 0
    if total_cn == 2:
        return 1, 1
    return max(1, total_cn - 1), 1


def read_shlush_cnv_table(path):
    """
    Read a Shlush copy_number_table CSV and return a list of segments:
        {
          'chrom': <string>,
          'start': <int>,
          'end': <int>,
          'total_cn': <int>,
        }
    """
    segs = []
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        required = ['contig', 'start', 'stop', 'cn']
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError("Missing required CNV columns in %s: %r" %
                             (path, missing))
        for row in reader:
            try:
                chrom = norm_chr(row['contig'])
                start = int(row['start'])
                end = int(row['stop'])
                total_cn = int(row['cn'])
            except (ValueError, TypeError):
                continue
            if end < start:
                continue
            segs.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'total_cn': total_cn,
            })
    return segs


def estimate_mean_depth(muts):
    """Estimate mean depth d from SSMs; fall back to a modest default."""
    ds = [m['d'] for m in muts if m.get('d', 0) > 0]
    if not ds:
        return 100.0
    return float(sum(ds)) / float(len(ds))


def build_ssm_data(muts, sample_label, out_path):
    """
    Write a single-sample PhyloWGS ssm_data.txt-style table from mutations.
    """
    if not muts:
        with open(out_path, 'w') as fh:
            fh.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        return

    dirname = os.path.dirname(out_path)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(out_path, 'w') as fh:
        fh.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        for idx, m in enumerate(muts):
            ssm_id = "s%d" % idx
            gene = "%s_%d" % (m['chrom'], m['pos'])
            a = m['a']
            d = m['d']
            fh.write("%s\t%s\t%d\t%d\t0.999\t0.5\n" %
                     (ssm_id, gene, a, d))


def build_cnv_data(segs, muts, out_path):
    """
    Write a single-sample PhyloWGS cnv_data.txt-style table from CN segments,
    assigning SSMs by genomic overlap.
    """
    dirname = os.path.dirname(out_path)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(out_path, 'w') as fh:
        fh.write("cnv\ta\td\tssms\tphysical_cnvs\n")
        if not segs:
            return

        mean_dp = estimate_mean_depth(muts)

        # Pre-index SSMs by chromosome for quick lookup
        muts_by_chr = {}
        for idx, m in enumerate(muts):
            muts_by_chr.setdefault(m['chrom'], []).append((idx, m))

        for i, seg in enumerate(segs):
            chrom = seg['chrom']
            start = seg['start']
            end = seg['end']
            total_cn = seg['total_cn']
            major_cn, minor_cn = decompose_total_cn(total_cn)

            # Approximate total reads over the segment
            length_bp = max(1, end - start + 1)
            d_est = int(max(1, round((length_bp / READ_LENGTH) * mean_dp)))
            if major_cn + minor_cn > 0:
                a_est = int(round(d_est * float(minor_cn) /
                                  float(max(1, major_cn + minor_cn))))
            else:
                a_est = 0

            # Assign SSMs to this segment
            ssms_entries = []
            if chrom in muts_by_chr:
                for ssm_idx, m in muts_by_chr[chrom]:
                    if start <= m['pos'] <= end:
                        ssm_id = "s%d" % ssm_idx
                        ssms_entries.append("%s,%d,%d" %
                                            (ssm_id, minor_cn, major_cn))
            ssms_str = ";".join(ssms_entries)

            physical_cnvs = (
                "chrom=%s,start=%d,end=%d,major_cn=%d,minor_cn=%d,cell_prev=%s"
                % (chrom, start, end, major_cn, minor_cn, "1")
            )

            cnv_id = "c%d" % i
            fh.write("%s\t%d\t%d\t%s\t%s\n" %
                     (cnv_id, a_est, d_est, ssms_str, physical_cnvs))


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Create PhyloWGS ssm_data.txt and cnv_data.txt from "
                    "Shlush preprocessed mutation and copy-number tables.")
    parser.add_argument('--mutation-table', required=True,
                        help='Path to Shlush *_mutation_table.csv')
    parser.add_argument('--copy-number-table', required=True,
                        help='Path to Shlush *_copy_number_table.csv')
    parser.add_argument('--sample-id', required=True,
                        help='Sample identifier label (e.g. Patient-1_Dx)')
    parser.add_argument('--out-mutation-dir', required=True,
                        help='Output directory for PhyloWGS-style mutation_table')
    parser.add_argument('--out-cnv-dir', required=True,
                        help='Output directory for PhyloWGS-style copy_number_file')

    args = parser.parse_args()

    muts = read_shlush_mutation_table(args.mutation_table)
    segs = read_shlush_cnv_table(args.copy_number_table)

    out_ssm_path = os.path.join(
        args.out_mutation_dir, '%s_ssm_data.txt' % args.sample_id
    )
    out_cnv_path = os.path.join(
        args.out_cnv_dir, '%s_cnv_data.txt' % args.sample_id
    )

    build_ssm_data(muts, args.sample_id, out_ssm_path)
    build_cnv_data(segs, muts, out_cnv_path)

    print("Wrote SSM data to:", out_ssm_path)
    print("Wrote CNV data to:", out_cnv_path)


if __name__ == '__main__':
    main()

