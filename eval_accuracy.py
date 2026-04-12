#!/usr/bin/env python3
"""
eval_accuracy.py – Compare a GRASP GAF against a PanAligner baseline GAF.

Two modes:
  1. --mode agreement  (default)
     For each read present in BOTH files, check whether the primary-alignment
     path's first node is the same.  Reports agreement rate.

  2. --mode position
     Reads must have been produced by simulate_reads.py, so their name is
     <seqname>_<start>_<end>.  For each aligned read, check whether the
     node the aligner chose contains the true reference position.
     Requires --gfa to look up node coordinates (SN/SO tags).

Usage:
  # Agreement rate (no ground truth needed):
  python3 eval_accuracy.py --mode agreement \
      --baseline  baseline_panaligner.gaf   \
      --query     grasp_output.gaf

  # Position accuracy (simulated reads + GFA with SN/SO):
  python3 eval_accuracy.py --mode position  \
      --gfa  pangenome.gfa                  \
      --query grasp_output.gaf
"""

import argparse
import sys
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# GAF helpers
# ---------------------------------------------------------------------------

def parse_gaf(path: str):
    """
    Return a dict  read_name -> primary GAF record (dict with fields).
    Primary = record with mapq >= 0 that is NOT marked secondary (tp:A:S).
    If multiple primaries, keep the one with the highest mapq.
    """
    records = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue
            rec = {
                'name':     fields[0],
                'qlen':     int(fields[1]),
                'qs':       int(fields[2]),
                'qe':       int(fields[3]),
                'strand':   fields[4],
                'path':     fields[5],
                'plen':     int(fields[6]),
                'ps':       int(fields[7]),
                'pe':       int(fields[8]),
                'nmatch':   int(fields[9]),
                'alen':     int(fields[10]),
                'mapq':     int(fields[11]),
                'tags':     {t[:2]: t for t in fields[12:] if ':' in t},
            }
            records[rec['name']].append(rec)

    primary = {}
    for name, recs in records.items():
        # Exclude secondary alignments (tp:A:S tag)
        candidates = [r for r in recs if r['tags'].get('tp', 'tp:A:P') != 'tp:A:S']
        if not candidates:
            candidates = recs
        best = max(candidates, key=lambda r: r['mapq'])
        primary[name] = best
    return primary


def first_node(path_str: str) -> str:
    """Extract the first node/chromosome name from a GAF path string.

    Handles two path formats PanAligner emits:
      compact:     NC_001139.9               (full-graph alignment)
      non-compact: >NC_001139.9:574777-620473 (sub-graph alignment)
    In both cases we want to return 'NC_001139.9'.
    Also handles plain node IDs like >s1>s2<s3.
    """
    m = re.match(r'[><]?([^><:]+)', path_str)
    return m.group(1) if m else path_str


def path_nodes(path_str: str):
    """Return list of node names from a GAF path string."""
    return re.findall(r'[><]([^><]+)', path_str) or [path_str]


# ---------------------------------------------------------------------------
# GFA helpers (for position mode)
# ---------------------------------------------------------------------------

def load_node_coords(gfa_path: str):
    """
    Parse GFA and return  node_name -> (seqname, start, end)
    using SN (sequence name) and SO (sequence offset) tags, plus LN (length).
    Nodes lacking SN/SO are skipped.
    """
    coords = {}
    with open(gfa_path) as fh:
        for line in fh:
            if not line.startswith('S'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            node_name = parts[1]
            tags = {}
            for t in parts[3:]:
                if ':' in t:
                    k = t[:2]
                    tags[k] = t
            sn_tag = tags.get('SN')
            so_tag = tags.get('SO')
            ln_tag = tags.get('LN')
            if sn_tag and so_tag and ln_tag:
                sn = sn_tag.split(':', 2)[2]
                so = int(so_tag.split(':', 2)[2])
                ln = int(ln_tag.split(':', 2)[2])
                coords[node_name] = (sn, so, so + ln)
    return coords


# ---------------------------------------------------------------------------
# Mode: agreement
# ---------------------------------------------------------------------------

def mode_agreement(args):
    print(f'[eval] Loading baseline: {args.baseline}')
    baseline = parse_gaf(args.baseline)
    print(f'[eval] Loading query:    {args.query}')
    query    = parse_gaf(args.query)

    common = set(baseline) & set(query)
    if not common:
        print('[eval] ERROR: No reads in common between baseline and query GAF.')
        sys.exit(1)

    agree = 0
    for name in common:
        b_node = first_node(baseline[name]['path'])
        q_node = first_node(query[name]['path'])
        if b_node == q_node:
            agree += 1

    rate = agree / len(common) * 100
    print(f'[eval] Reads in common : {len(common)}')
    print(f'[eval] Agree (same 1st node): {agree}')
    print(f'[eval] Agreement rate  : {rate:.2f}%')
    print(f'AGREEMENT_RATE={rate:.4f}')

    # Also report mapping rates
    baseline_mapped = len(baseline)
    query_mapped    = len(query)
    print(f'[eval] Baseline mapped reads: {baseline_mapped}')
    print(f'[eval] Query mapped reads   : {query_mapped}')
    if baseline_mapped > 0:
        rel = query_mapped / baseline_mapped * 100
        print(f'[eval] Relative mapping rate: {rel:.2f}%')
        print(f'MAPPING_RATE_BASELINE={baseline_mapped}')
        print(f'MAPPING_RATE_QUERY={query_mapped}')
        print(f'RELATIVE_MAPPING_RATE={rel:.4f}')


# ---------------------------------------------------------------------------
# Mode: position  (simulated reads with seqname_start_end names)
# ---------------------------------------------------------------------------

def mode_position(args):
    if not args.gfa:
        sys.exit('[eval] --gfa is required for position mode')

    print(f'[eval] Loading GFA node coords: {args.gfa}')
    coords = load_node_coords(args.gfa)
    if not coords:
        print('[eval] WARNING: No SN/SO tags found. Cannot check positions.')
        sys.exit(1)

    print(f'[eval] Loading query: {args.query}')
    query = parse_gaf(args.query)

    total = 0
    correct = 0
    unmapped = 0

    name_re = re.compile(r'^(.+)_(\d+)_(\d+)$')

    for read_name, rec in query.items():
        m = name_re.match(read_name)
        if not m:
            continue  # not a simulated read
        true_seq   = m.group(1)
        true_start = int(m.group(2))
        true_end   = int(m.group(3))
        total += 1

        nodes = path_nodes(rec['path'])
        hit = False
        for node in nodes:
            if node not in coords:
                continue
            n_seq, n_start, n_end = coords[node]
            # Correct if the node's reference interval overlaps the read's origin
            if n_seq == true_seq and n_start < true_end and n_end > true_start:
                hit = True
                break
        if hit:
            correct += 1

    rate = correct / total * 100 if total > 0 else 0.0
    print(f'[eval] Simulated reads evaluated: {total}')
    print(f'[eval] Correctly mapped          : {correct}')
    print(f'[eval] Correct-mapping rate      : {rate:.2f}%')
    print(f'CORRECT_MAPPING_RATE={rate:.4f}')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode',     choices=['agreement', 'position'],
                    default='agreement')
    ap.add_argument('--baseline', help='Baseline GAF (PanAligner, for agreement mode)')
    ap.add_argument('--query',    required=True, help='GAF to evaluate (GRASP output)')
    ap.add_argument('--gfa',      help='GFA graph with SN/SO node tags (for position mode)')
    args = ap.parse_args()

    if args.mode == 'agreement':
        if not args.baseline:
            sys.exit('[eval] --baseline is required for agreement mode')
        mode_agreement(args)
    else:
        mode_position(args)


if __name__ == '__main__':
    main()
