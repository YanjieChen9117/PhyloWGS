#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Construct multi-sample PhyloWGS inputs for Patient-6 (Dx + Rx) and run
PhyloWGS (multievolve.py) on the combined data.

Inputs (already in single-sample PhyloWGS format, created previously):

  Shlush_multi-sample-test/mutation_table/Patient-6_Dx_ssm_data.txt
  Shlush_multi-sample-test/mutation_table/Patient-6_Rx_ssm_data.txt
  Shlush_multi-sample-test/copy_number_file/Patient-6_Dx_cnv_data.txt
  Shlush_multi-sample-test/copy_number_file/Patient-6_Rx_cnv_data.txt

Outputs (multi-sample format, 2 samples: Dx,Rx):

  Shlush_multi-sample-test/mutation_table/Patient-6_multi_ssm_data.txt
  Shlush_multi-sample-test/copy_number_file/Patient-6_multi_cnv_data.txt

The multi-sample formats follow the examples in:
  - parser/test/outputs/multisamp_cnvs/ssm_data.txt
  - parser/test/outputs/multisamp_cnvs/cnv_data.txt
"""

from __future__ import print_function

import csv
import os
import subprocess as sp


ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def read_single_ssm(path):
    """
    Read a single-sample ssm_data.txt into a mapping:
      gene -> (a, d)
    Also keep chromosome and position parsed from gene string "chr_pos".
    """
    by_gene = {}
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene = row['gene']
            try:
                chrom_str, pos_str = gene.split('_', 1)
                chrom = chrom_str
                pos = int(pos_str)
            except Exception:
                # Skip malformed entries
                continue
            try:
                a = int(row['a'])
                d = int(row['d'])
            except ValueError:
                continue
            by_gene[gene] = {
                'chrom': chrom,
                'pos': pos,
                'a': a,
                'd': d,
            }
    return by_gene


def build_multi_ssm(dx_path, rx_path, out_path):
    """
    Combine Dx and Rx single-sample SSMs into a two-sample ssm_data.txt.
    We take the union of loci (by `gene` column).
    For loci missing in one sample, we set a=0, d=1.
    """
    dx = read_single_ssm(dx_path)
    rx = read_single_ssm(rx_path)

    genes = sorted(set(dx.keys()) | set(rx.keys()))

    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(out_path, 'w') as fh:
        fh.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        for idx, gene in enumerate(genes):
            dx_rec = dx.get(gene)
            rx_rec = rx.get(gene)

            if dx_rec is None:
                a_dx, d_dx = 0, 1
            else:
                a_dx, d_dx = dx_rec['a'], dx_rec['d']

            if rx_rec is None:
                a_rx, d_rx = 0, 1
            else:
                a_rx, d_rx = rx_rec['a'], rx_rec['d']

            ssm_id = "s%d" % idx
            a_str = "%d,%d" % (a_dx, a_rx)
            d_str = "%d,%d" % (d_dx, d_rx)
            fh.write("%s\t%s\t%s\t%s\t0.999\t0.5\n" %
                     (ssm_id, gene, a_str, d_str))

    # Also return positional info for CNV assignment
    pos_info = {}
    for idx, gene in enumerate(genes):
        chrom_str, pos_str = gene.split('_', 1)
        pos_info[idx] = {
            'id': "s%d" % idx,
            'chrom': chrom_str,
            'pos': int(pos_str),
        }
    return pos_info


def parse_physical_cnvs(phys_str):
    """
    Parse physical_cnvs string of the form:
      chrom=1,start=1234,end=5678,major_cn=2,minor_cn=1,cell_prev=1
    or with multi-sample cell_prev:
      cell_prev=0.0|0.718
    Return dict with keys: chrom, start, end, major_cn, minor_cn, cell_prev.
    """
    parts = phys_str.strip().split(',')
    out = {}
    for p in parts:
        if '=' not in p:
            continue
        k, v = p.split('=', 1)
        out[k] = v
    # normalize types where needed later
    return out


def read_single_cnv(path):
    """
    Read a single-sample cnv_data.txt into:
      key -> { 'a': int, 'd': int, 'chrom': str, 'start': int,
               'end': int, 'major_cn': int, 'minor_cn': int,
               'cell_prev': str }
    where key = (chrom, start, end, major_cn, minor_cn).
    """
    segs = {}
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            phys = parse_physical_cnvs(row['physical_cnvs'])
            try:
                chrom = phys['chrom']
                start = int(phys['start'])
                end = int(phys['end'])
                major_cn = int(phys['major_cn'])
                minor_cn = int(phys['minor_cn'])
                a = int(row['a'])
                d = int(row['d'])
            except Exception:
                continue
            cell_prev = phys.get('cell_prev', '1')
            key = (chrom, start, end, major_cn, minor_cn)
            segs[key] = {
                'a': a,
                'd': d,
                'chrom': chrom,
                'start': start,
                'end': end,
                'major_cn': major_cn,
                'minor_cn': minor_cn,
                'cell_prev': cell_prev,
            }
    return segs


def build_multi_cnv(dx_path, rx_path, ssm_pos_info, out_path):
    """
    Combine Dx and Rx single-sample CNVs into a two-sample cnv_data.txt.

    - Union of CNV segments keyed by (chrom,start,end,major_cn,minor_cn).
    - For missing segments in a sample, set a=0, d=1, cell_prev=0.
    - Recompute `ssms` by overlapping SSM loci with segments.
    - `physical_cnvs` has cell_prev as "dx|rx".
    """
    dx = read_single_cnv(dx_path)
    rx = read_single_cnv(rx_path)

    keys = sorted(set(dx.keys()) | set(rx.keys()))

    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(out_path, 'w') as fh:
        fh.write("cnv\ta\td\tssms\tphysical_cnvs\n")

        for idx, key in enumerate(keys):
            chrom, start, end, major_cn, minor_cn = key
            dx_seg = dx.get(key)
            rx_seg = rx.get(key)

            if dx_seg is None:
                a_dx, d_dx, prev_dx = 0, 1, "0"
            else:
                a_dx = dx_seg['a']
                d_dx = dx_seg['d']
                prev_dx = dx_seg['cell_prev']

            if rx_seg is None:
                a_rx, d_rx, prev_rx = 0, 1, "0"
            else:
                a_rx = rx_seg['a']
                d_rx = rx_seg['d']
                prev_rx = rx_seg['cell_prev']

            a_str = "%d,%d" % (a_dx, a_rx)
            d_str = "%d,%d" % (d_dx, d_rx)

            # Assign SSMs to this segment by genomic overlap
            ssms_entries = []
            for s_idx, info in ssm_pos_info.items():
                if info['chrom'] != chrom:
                    continue
                if start <= info['pos'] <= end:
                    ssm_id = info['id']
                    ssms_entries.append("%s,%d,%d" % (ssm_id, minor_cn, major_cn))
            ssms_str = ";".join(ssms_entries)

            cell_prev = "%s|%s" % (prev_dx, prev_rx)
            phys = (
                "chrom=%s,start=%d,end=%d,major_cn=%d,minor_cn=%d,cell_prev=%s"
                % (chrom, start, end, major_cn, minor_cn, cell_prev)
            )

            cnv_id = "c%d" % idx
            fh.write("%s\t%s\t%s\t%s\t%s\n" %
                     (cnv_id, a_str, d_str, ssms_str, phys))


def run_multievolve(ssm_path, cnv_path, out_dir,
                    num_chains=4, burnin=500, mcmc=1000):
    """
    Run PhyloWGS multievolve.py on the provided multi-sample inputs.
    Also capture stdout/stderr into files named:
        phylowgs_stdout.txt
        phylowgs_stderr.txt
    inside the output directory, similar to PhyloWGS_runtime_analysis.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cmd = [
        'python2',
        os.path.join(ROOT, 'multievolve.py'),
        '--num-chains', str(num_chains),
        '--burnin-samples', str(burnin),
        '--mcmc-samples', str(mcmc),
        '--ssms', ssm_path,
        '--cnvs', cnv_path,
        '--output-dir', out_dir,
    ]
    print("Running:", " ".join(cmd))

    stdout_path = os.path.join(out_dir, 'phylowgs_stdout.txt')
    stderr_path = os.path.join(out_dir, 'phylowgs_stderr.txt')

    with open(stdout_path, 'w') as fh_out, open(stderr_path, 'w') as fh_err:
        proc = sp.Popen(cmd, cwd=ROOT, stdout=fh_out, stderr=fh_err)
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(
                "multievolve.py exited with code %d; "
                "see %s and %s for details" %
                (ret, stdout_path, stderr_path)
            )


def main():
    dx_ssm = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'mutation_table',
        'Patient-6_Dx_ssm_data.txt'
    )
    rx_ssm = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'mutation_table',
        'Patient-6_Rx_ssm_data.txt'
    )
    dx_cnv = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'copy_number_file',
        'Patient-6_Dx_cnv_data.txt'
    )
    rx_cnv = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'copy_number_file',
        'Patient-6_Rx_cnv_data.txt'
    )

    out_ssm = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'mutation_table',
        'Patient-6_multi_ssm_data.txt'
    )
    out_cnv = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'copy_number_file',
        'Patient-6_multi_cnv_data.txt'
    )

    print("Building multi-sample SSM for Patient-6 ...")
    ssm_pos_info = build_multi_ssm(dx_ssm, rx_ssm, out_ssm)

    print("Building multi-sample CNV for Patient-6 ...")
    build_multi_cnv(dx_cnv, rx_cnv, ssm_pos_info, out_cnv)

    chains_dir = os.path.join(
        ROOT, 'Shlush_multi-sample-test', 'Patient-6_multisample_chains'
    )
    print("Running PhyloWGS multievolve on Patient-6 multi-sample inputs ...")
    run_multievolve(out_ssm, out_cnv, chains_dir)

    print("Done.")
    print("SSM file:", out_ssm)
    print("CNV file:", out_cnv)
    print("Chains output dir:", chains_dir)


if __name__ == '__main__':
    main()

