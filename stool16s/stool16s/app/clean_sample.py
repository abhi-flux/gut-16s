#!/usr/bin/env python3
"""
Clean a raw stool sequence file (FASTA or FASTQ) so classify.py can read it.

Usage:
  python clean_sample.py --input raw_sample.fastq --output stool_sample.fasta
"""
import argparse
import re
import sys
from collections import Counter

VALID = set("ACGT")


def read_records(path):
    """Yield (header, sequence, quality) from FASTA or FASTQ. Quality is None for FASTA."""
    with open(path) as f:
        first = f.readline()
        f.seek(0)
        if first.startswith(">"):
            header, seq_lines = None, []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_lines), None
                    header, seq_lines = line[1:], []
                else:
                    seq_lines.append(line)
            if header is not None:
                yield header, "".join(seq_lines), None
        elif first.startswith("@"):
            # Standard 4-line FASTQ records (sequence must be on one line)
            while True:
                h = f.readline().strip()
                if not h:
                    break
                seq = f.readline().strip()
                f.readline()  # the '+' line
                qual = f.readline().strip()
                yield h[1:], seq, qual
        else:
            sys.exit("Input is neither FASTA nor FASTQ. The first line must start with > or @.")


def clean_id(raw):
    parts = raw.split()
    name = parts[0] if parts else "seq"
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def mean_quality(qual):
    """Average Phred score (Sanger/Illumina 1.8+ encoding)."""
    return sum(ord(c) - 33 for c in qual) / len(qual)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="raw sample file (FASTA or FASTQ)")
    parser.add_argument("--output", default="stool_sample.fasta", help="cleaned FASTA file")
    parser.add_argument("--min-len", type=int, default=200, help="drop reads shorter than this")
    parser.add_argument("--max-len", type=int, default=2000, help="drop reads longer than this")
    parser.add_argument("--min-qual", type=float, default=20.0,
                        help="FASTQ only: drop reads with average quality below this")
    args = parser.parse_args()

    kept = []
    reasons = Counter()
    used_ids = Counter()

    for header, seq, qual in read_records(args.input):
        reasons["total_read"] += 1
        seq = seq.upper().replace("U", "T").replace(" ", "")

        if any(base not in VALID for base in seq):
            reasons["removed: non-ACGT bases (e.g. N)"] += 1
            continue
        if len(seq) < args.min_len:
            reasons["removed: too short"] += 1
            continue
        if len(seq) > args.max_len:
            reasons["removed: too long"] += 1
            continue
        if qual is not None and mean_quality(qual) < args.min_qual:
            reasons["removed: low quality"] += 1
            continue

        seq_id = clean_id(header)
        used_ids[seq_id] += 1
        if used_ids[seq_id] > 1:
            seq_id = f"{seq_id}_{used_ids[seq_id]}"  # make duplicate IDs unique

        kept.append((seq_id, seq))
        reasons["kept"] += 1

    with open(args.output, "w") as out:
        for seq_id, seq in kept:
            out.write(f">{seq_id}\n{seq}\n")

    print(f"Cleaned file written to {args.output}")
    print(f"Reads in input: {reasons.pop('total_read', 0)}")
    for reason, n in sorted(reasons.items()):
        print(f"  {reason}: {n}")


if __name__ == "__main__":
    main()