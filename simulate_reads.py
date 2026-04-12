#!/usr/bin/env python3
"""
simulate_reads.py  –  Generate simulated long reads from a FASTA reference.

Read names encode the true position so eval_accuracy.py can later check
whether an alignment is correct:
  <seqname>_<start>_<end>

Usage:
  python3 simulate_reads.py \
      --fasta  reference.fna \
      --output reads.fasta   \
      --length 10000         \
      --step   5000          \
      --num    5000          \
      --seed   42
"""

import argparse
import random
import sys
from pathlib import Path


def read_fasta(path: str):
    """Return a dict {seqname: sequence_str}."""
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def add_errors(seq: str, error_rate: float, rng: random.Random) -> str:
    """Add substitution errors at *error_rate* fraction of bases."""
    if error_rate <= 0:
        return seq
    bases = list(seq)
    for i, b in enumerate(bases):
        if rng.random() < error_rate:
            bases[i] = rng.choice('ACGT')
    return ''.join(bases)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta',  required=True, help='Input reference FASTA')
    ap.add_argument('--output', required=True, help='Output simulated-reads FASTA')
    ap.add_argument('--length', type=int, default=10000, help='Read length (default 10000)')
    ap.add_argument('--step',   type=int, default=0,
                    help='Step between reads; 0=random sampling (default 0)')
    ap.add_argument('--num',    type=int, default=5000,
                    help='Max number of reads to generate (default 5000)')
    ap.add_argument('--error',  type=float, default=0.0,
                    help='Per-base substitution error rate 0-1 (default 0)')
    ap.add_argument('--seed',   type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    seqs = read_fasta(args.fasta)

    if not seqs:
        sys.exit(f'[ERROR] No sequences found in {args.fasta}')

    L = args.length
    reads_written = 0

    with open(args.output, 'w') as out:
        if args.step > 0:
            # Tiled sampling
            for seqname, seq in seqs.items():
                pos = 0
                while pos + L <= len(seq) and reads_written < args.num:
                    subseq = seq[pos:pos + L].upper()
                    if 'N' * (L // 2) not in subseq:   # skip runs of N
                        read = add_errors(subseq, args.error, rng)
                        out.write(f'>{seqname}_{pos}_{pos+L}\n{read}\n')
                        reads_written += 1
                    pos += args.step
                if reads_written >= args.num:
                    break
        else:
            # Random sampling
            items = [(n, s) for n, s in seqs.items() if len(s) >= L]
            if not items:
                sys.exit(f'[ERROR] All sequences are shorter than --length {L}')
            total_len = sum(len(s) for _, s in items)
            weights   = [len(s) / total_len for _, s in items]

            attempts = 0
            while reads_written < args.num and attempts < args.num * 10:
                attempts += 1
                seqname, seq = rng.choices(items, weights=weights, k=1)[0]
                if len(seq) < L:
                    continue
                pos = rng.randint(0, len(seq) - L)
                subseq = seq[pos:pos + L].upper()
                if subseq.count('N') > L * 0.1:   # skip >10% N
                    continue
                read = add_errors(subseq, args.error, rng)
                out.write(f'>{seqname}_{pos}_{pos+L}\n{read}\n')
                reads_written += 1

    print(f'[simulate_reads] Wrote {reads_written} reads to {args.output}')


if __name__ == '__main__':
    main()
