#!/usr/bin/env python3
"""
Pipeline stage 1A: expand the experimental design into a permutation plan.

Input: the tested six-column community manifest plus target abundances, one
depth definition (coverage, total bases, or paired-read count), and seeds.
Output: one TSV row per abundance x depth x seed combination. Downstream PBS
array jobs use ``run_index`` as their 1-based row selector.

Create a flat RandomReadsMG permutation plan.

One row = one randomreads virtual-metagenome run.

This version keeps the execution columns used by run_rrmg_plan_array.sh:
  total_cov / total_bases / total_pairs
and adds design-check columns:
  community_size_bp, computed_total_bases, computed_total_pairs,
  computed_total_cov, target_genome_bp, target_expected_cov, etc.

Python 3.6 compatible.

``community_size_bp`` is a genome-equivalent size, not the sum of every FASTA.
Without an explicit reference/size it is the abundance-weighted mean genome
length. In ``read_fraction`` mode, target read pairs are proportional to the
requested target abundance.
"""

import argparse
import csv
import gzip
import itertools
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_list(s, cast=str):
    out = []
    for x in str(s).split(','):
        x = x.strip()
        if x:
            out.append(cast(x))
    return out


def parse_seeds(s):
    seeds = []
    for part in str(s).split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            a, b = part.split('-', 1)
            a_i, b_i = int(a), int(b)
            step = 1 if b_i >= a_i else -1
            seeds.extend(list(range(a_i, b_i + step, step)))
        else:
            seeds.append(int(part))
    return seeds


def safe_tag(s):
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return x.strip('_') or 'NA'


def norm_key(s):
    return re.sub(r"[\s_\-]+", "", str(s)).lower()


def open_text(path):
    path = Path(path)
    if str(path).endswith('.gz'):
        return gzip.open(str(path), 'rt', errors='ignore')
    return path.open('r', errors='ignore')


def fasta_len_bp(path, acgt_only=False):
    path = Path(path)
    if not path.exists() or not path.is_file():
        raise SystemExit("[ERROR] FASTA not found: {}".format(path))
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith('>'):
                continue
            s = line.strip().upper()
            if acgt_only:
                n += sum(1 for c in s if c in 'ACGT')
            else:
                n += len(s)
    if n <= 0:
        raise SystemExit("[ERROR] FASTA length is 0: {}".format(path))
    return n


def sniff_dict_reader(path):
    path = Path(path)
    with path.open(newline='') as fh:
        sample = fh.read(8192)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters='\t,;')
        except csv.Error:
            dialect = csv.excel_tab
        return list(csv.DictReader(fh, dialect=dialect))


def pick(row, names, default=''):
    m = {norm_key(k): k for k in row.keys()}
    for name in names:
        nk = norm_key(name)
        if nk in m:
            return row.get(m[nk], default)
    return default


def read_manifest(path):
    raw = sniff_dict_reader(path)
    rows = []
    for i, r in enumerate(raw, start=1):
        gid = pick(r, ['genome_ID', 'genome', 'id']) or 'Genome{}'.format(i)
        name = pick(r, ['name', 'Name', 'species', 'taxon']) or gid
        taxid = pick(r, ['taxid', 'NCBI ID', 'NCBI_ID', 'TaxID']) or 'NA'
        rel_s = pick(r, ['rel_abund', 'relative_abundance', 'abundance', 'rel']) or '0'
        fasta = pick(r, ['fasta', 'path_to_fasta_file', 'genome_path', 'path'])
        role = (pick(r, ['role']) or 'background').lower()
        try:
            rel = float(rel_s)
        except Exception:
            rel = 0.0
        rows.append({
            'genome_ID': gid,
            'name': name,
            'taxid': str(taxid),
            'rel_abund': max(0.0, rel),
            'fasta': fasta,
            'role': role,
        })
    if not rows:
        raise SystemExit('[ERROR] Empty manifest')
    return rows


def is_target_row(row, target_name='', target_taxid='', target_genome_id=''):
    if str(row.get('role', '')).lower() == 'target':
        return True
    if target_genome_id and str(row.get('genome_ID', '')) == str(target_genome_id):
        return True
    if target_taxid and str(target_taxid) not in {'NA', '0', ''} and str(row.get('taxid', '')) == str(target_taxid):
        return True
    if target_name and str(row.get('name', '')).strip().lower() == str(target_name).strip().lower():
        return True
    return False


def resolve_target_fasta(rows, target_name='', target_taxid='', target_genome_id='', target_fasta=''):
    if target_fasta:
        return target_fasta
    for r in rows:
        if is_target_row(r, target_name, target_taxid, target_genome_id):
            if r.get('fasta'):
                return r['fasta']
    return ''


def final_abundance_rows(rows, target_abund, target_name='', target_taxid='', target_genome_id='', target_fasta=''):
    """Mimic rrmg_virtual_metagenome_permutation.sh insertion/rescaling.

    Target abundance can be 0. In that case target is not included in the final
    community for simulation, but target coverage columns can still be computed
    from target_fasta if available.
    """
    target_abund = float(target_abund)
    if not (0 <= target_abund < 1):
        raise SystemExit('[ERROR] target abundance must be >=0 and <1, got {}'.format(target_abund))

    background = [dict(r) for r in rows if not is_target_row(r, target_name, target_taxid, target_genome_id)]
    bg_total = sum(r['rel_abund'] for r in background if r['rel_abund'] > 0)
    if bg_total <= 0 and target_abund < 1:
        raise SystemExit('[ERROR] background relative abundance total is 0')

    final = []
    for r in background:
        if r['rel_abund'] <= 0:
            continue
        rr = dict(r)
        rr['rel_abund'] = rr['rel_abund'] / bg_total * (1.0 - target_abund)
        rr['role'] = 'background'
        final.append(rr)

    if target_abund > 0:
        # Prefer the manifest target row if present; otherwise synthesize from args.
        target_candidates = [r for r in rows if is_target_row(r, target_name, target_taxid, target_genome_id)]
        t = dict(target_candidates[0]) if target_candidates else {}
        t['genome_ID'] = target_genome_id or t.get('genome_ID') or safe_tag(target_name or 'Target')
        t['name'] = target_name or t.get('name') or t['genome_ID']
        t['taxid'] = target_taxid or t.get('taxid') or 'NA'
        t['fasta'] = target_fasta or t.get('fasta') or ''
        t['role'] = 'target'
        t['rel_abund'] = target_abund
        if not t['fasta']:
            raise SystemExit('[ERROR] target FASTA is required for target_abundance > 0')
        final.append(t)

    # Correct numerical drift.
    total = sum(r['rel_abund'] for r in final)
    if total <= 0:
        raise SystemExit('[ERROR] final community abundance total is 0')
    for r in final:
        r['rel_abund'] = r['rel_abund'] / total
    return final


def weighted_mean_len(rows, acgt_only=False):
    total = 0.0
    for r in rows:
        if not r.get('fasta'):
            raise SystemExit('[ERROR] Missing FASTA path for {}'.format(r.get('genome_ID')))
        L = fasta_len_bp(r['fasta'], acgt_only=acgt_only)
        total += r['rel_abund'] * L
    if total <= 0:
        raise SystemExit('[ERROR] weighted community genome size is 0')
    return total


def format_float(x, digits=10):
    if x == '' or x is None:
        return ''
    return ('{0:.%df}' % digits).format(float(x))


def main():
    ap = argparse.ArgumentParser(description='Create rrmg_permutation_plan.tsv for PBS-array randomreads simulations with precomputed depth/coverage columns.')
    ap.add_argument('--manifest', required=True, help='Flat community manifest: genome_ID,name,taxid,rel_abund,fasta,role')
    ap.add_argument('--out', required=True, help='Output plan TSV')
    ap.add_argument('--sample-prefix', required=True, help='Prefix for run/sample names')
    ap.add_argument('--target-name', required=True, help='Target species name, e.g. Aspergillus fumigatus')
    ap.add_argument('--target-taxid', default='', help='Target TaxID')
    ap.add_argument('--target-genome-id', default='', help='Target genome_ID in manifest, optional')
    ap.add_argument('--target-fasta', default='', help='Target genome FASTA path, optional if target row is in manifest')
    ap.add_argument('--target-abundances', required=True, help='Comma-separated target relative abundances, e.g. 0,0.0001,0.001,0.01')
    ap.add_argument('--coverages', default='', help='Comma-separated total community coverages, e.g. 1,5,20')
    ap.add_argument('--total-bases', default='', help='Comma-separated total read bases instead of coverages')
    ap.add_argument('--total-pairs', default='', help='Comma-separated total PE read pairs instead of coverages/total-bases')
    ap.add_argument('--seeds', required=True, help='Comma-separated or ranged seeds, e.g. 1-100 or 1,2,3')
    ap.add_argument('--readlen', type=int, default=150, help='Single-end read length for PE reads; PE bases = pairs*readlen*2, default 150')
    ap.add_argument('--reference-fasta', '--community-reference-fasta', dest='reference_fasta', default='',
                    help='Optional reference/community FASTA to define community_size_bp for coverage calculation')
    ap.add_argument('--community-size-bp', default='', help='Optional numeric community genome-equivalent size in bp; ignored if --reference-fasta is set')
    ap.add_argument('--count-acgt-only', action='store_true', help='Count only A/C/G/T bases in FASTA length calculations; default counts all non-header characters including N')
    ap.add_argument('--abundance-mode', choices=['read_fraction', 'copy_number'], default='read_fraction',
                    help='Coverage calculation mode; should match rrmg_virtual_metagenome_permutation.sh, default read_fraction')
    ap.add_argument('--species-name-for-kraken', default='', help='Species name to extract from Kraken2 report; default target-name')
    args = ap.parse_args()

    depth_modes = sum(bool(x) for x in [args.coverages, args.total_bases, args.total_pairs])
    if depth_modes != 1:
        raise SystemExit('[ERROR] Provide exactly one of --coverages, --total-bases, or --total-pairs')
    if args.readlen <= 0:
        raise SystemExit('[ERROR] --readlen must be > 0')

    rows = read_manifest(args.manifest)
    abundances = parse_list(args.target_abundances, float)
    seeds = parse_seeds(args.seeds)

    if args.coverages:
        depths = [('total_cov', x) for x in parse_list(args.coverages, float)]
    elif args.total_bases:
        depths = [('total_bases', x) for x in parse_list(args.total_bases, float)]
    else:
        depths = [('total_pairs', x) for x in parse_list(args.total_pairs, float)]

    if not abundances or not seeds or not depths:
        raise SystemExit('[ERROR] Empty abundances/depths/seeds')

    fixed_reference_bp = None
    reference_fasta = args.reference_fasta
    reference_source = 'weighted_mean_from_manifest'
    if reference_fasta:
        fixed_reference_bp = fasta_len_bp(reference_fasta, acgt_only=args.count_acgt_only)
        reference_source = 'reference_fasta'
    elif args.community_size_bp:
        fixed_reference_bp = float(args.community_size_bp)
        if fixed_reference_bp <= 0:
            raise SystemExit('[ERROR] --community-size-bp must be > 0')
        reference_source = 'manual_community_size_bp'

    target_fasta_resolved = resolve_target_fasta(rows, args.target_name, args.target_taxid, args.target_genome_id, args.target_fasta)
    target_genome_bp = ''
    if target_fasta_resolved:
        target_genome_bp = fasta_len_bp(target_fasta_resolved, acgt_only=args.count_acgt_only)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        # execution columns used by run_rrmg_plan_array.sh
        'run_index', 'run_id', 'sample', 'manifest',
        'target_name', 'target_taxid', 'target_genome_id', 'target_fasta', 'target_abundance',
        'total_cov', 'total_bases', 'total_pairs', 'seed', 'species_name_for_kraken',
        # new design/check columns
        'readlen', 'abundance_mode', 'reference_fasta', 'reference_source', 'community_size_bp',
        'computed_total_bases', 'computed_total_pairs', 'computed_total_cov',
        'target_fasta_resolved', 'target_genome_bp', 'target_expected_bases',
        'target_expected_pairs', 'target_expected_cov',
    ]

    idx = 0
    with out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()

        for abund, (depth_type, depth_value), seed in itertools.product(abundances, depths, seeds):
            idx += 1

            final_rows = final_abundance_rows(
                rows, abund,
                target_name=args.target_name,
                target_taxid=args.target_taxid,
                target_genome_id=args.target_genome_id,
                target_fasta=target_fasta_resolved,
            )
            if fixed_reference_bp is not None:
                community_size_bp = float(fixed_reference_bp)
            else:
                community_size_bp = float(weighted_mean_len(final_rows, acgt_only=args.count_acgt_only))

            # Execution depth columns: keep exactly one populated, so old runner/simulator remain safe.
            total_cov = total_bases = total_pairs = ''
            depth_tag = ''
            if depth_type == 'total_cov':
                total_cov = depth_value
                computed_total_cov = float(depth_value)
                computed_total_bases = computed_total_cov * community_size_bp
                computed_total_pairs = computed_total_bases / (args.readlen * 2.0)
                depth_tag = 'cov{}X'.format(safe_tag(depth_value))
            elif depth_type == 'total_bases':
                total_bases = depth_value
                computed_total_bases = float(depth_value)
                computed_total_pairs = computed_total_bases / (args.readlen * 2.0)
                computed_total_cov = computed_total_bases / community_size_bp
                depth_tag = 'bases{}'.format(safe_tag(depth_value))
            else:
                total_pairs = depth_value
                computed_total_pairs = float(depth_value)
                computed_total_bases = computed_total_pairs * args.readlen * 2.0
                computed_total_cov = computed_total_bases / community_size_bp
                depth_tag = 'pairs{}'.format(safe_tag(depth_value))

            if args.abundance_mode == 'copy_number':
                target_expected_cov = computed_total_cov * float(abund)
                target_expected_bases = float(target_expected_cov) * float(target_genome_bp) if target_genome_bp != '' else ''
            else:
                target_expected_bases = computed_total_bases * float(abund)
                target_expected_cov = (target_expected_bases / float(target_genome_bp)) if target_genome_bp not in ('', 0) else ''
            target_expected_pairs = (target_expected_bases / (args.readlen * 2.0)) if target_expected_bases != '' else ''

            abund_tag = 'ta{}'.format(safe_tag(abund))
            run_id = '{}__{}__{}__Seed{}'.format(safe_tag(args.sample_prefix), abund_tag, depth_tag, seed)
            sample = '{}_{}_{}'.format(safe_tag(args.sample_prefix), abund_tag, depth_tag)

            w.writerow({
                'run_index': idx,
                'run_id': run_id,
                'sample': sample,
                'manifest': args.manifest,
                'target_name': args.target_name,
                'target_taxid': args.target_taxid,
                'target_genome_id': args.target_genome_id,
                'target_fasta': args.target_fasta,
                'target_abundance': format_float(abund, 10),
                'total_cov': format_float(total_cov, 10) if total_cov != '' else '',
                'total_bases': format_float(total_bases, 3) if total_bases != '' else '',
                'total_pairs': str(int(round(float(total_pairs)))) if total_pairs != '' else '',
                'seed': seed,
                'species_name_for_kraken': args.species_name_for_kraken or args.target_name,
                'readlen': args.readlen,
                'abundance_mode': args.abundance_mode,
                'reference_fasta': reference_fasta,
                'reference_source': reference_source,
                'community_size_bp': format_float(community_size_bp, 3),
                'computed_total_bases': format_float(computed_total_bases, 3),
                'computed_total_pairs': format_float(computed_total_pairs, 3),
                'computed_total_cov': format_float(computed_total_cov, 10),
                'target_fasta_resolved': target_fasta_resolved,
                'target_genome_bp': target_genome_bp,
                'target_expected_bases': format_float(target_expected_bases, 3) if target_expected_bases != '' else '',
                'target_expected_pairs': format_float(target_expected_pairs, 3) if target_expected_pairs != '' else '',
                'target_expected_cov': format_float(target_expected_cov, 10) if target_expected_cov != '' else '',
            })

    print('[OK] Wrote plan: {}'.format(out))
    print('[OK] Runs: {}'.format(idx))
    print('[OK] Submit PBS array as: #PBS -J 1-{}'.format(idx))
    if fixed_reference_bp is not None:
        print('[OK] community_size_bp used for plan: {:.3f} ({})'.format(float(fixed_reference_bp), reference_source))
    else:
        print('[OK] community_size_bp was computed per row as weighted mean genome length from manifest')
    if target_genome_bp != '':
        print('[OK] target_genome_bp: {} ({})'.format(target_genome_bp, target_fasta_resolved))
    else:
        print('[WARN] target genome length was not computed; target_expected_cov columns will be empty')


if __name__ == '__main__':
    main()
