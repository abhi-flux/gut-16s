#!/usr/bin/env python3
"""
Classify 16S sequences in a stool sample as one of your target bacteria or Unknown.

Usage:
  python classify.py --targets targets.txt --sample stool_sample.fasta
  python classify.py --targets targets.txt --sample stool_sample.fasta --skip-fetch
"""
import argparse
import csv
import subprocess
import sys
import time
from collections import Counter

from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com"  # NCBI requires a real email
# Entrez.api_key = "YOUR_KEY"            # optional: makes downloads faster

GENUS_MIN_ID = 97.0      # identity needed for a genus-level label
SPECIES_MIN_ID = 98.7    # identity needed for a species-level label
COVERAGE_MIN = 80.0      # % of the sample sequence that must align
EVALUE_MAX = 1e-20
MAX_REFS_PER_GENUS = 15  # reference sequences to download for "spp." entries
REFS_PER_SPECIES = 3     # reference sequences to download per named species
MIN_LEN, MAX_LEN = 1200, 1700
REF_FASTA = "references.fasta"
DB_NAME = "refs_db"
UNKNOWN = "Unknown"

# Alternative names NCBI may use for the same organism
ALIASES = {
    "Eubacterium rectale": ["Eubacterium rectale", "Agathobacter rectalis"],
    "Clostridium difficile": ["Clostridium difficile", "Clostridioides difficile"],
    "Ruminococcus gnavus": ["Ruminococcus gnavus", "Mediterraneibacter gnavus"],
}


def read_targets(path):
    """Return {genus: {"species": set of full names, "any": bool}}."""
    targets = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            genus, _, rest = line.partition(":")
            genus = genus.strip()
            entry = {"species": set(), "any": False}
            for item in rest.split(","):
                item = item.strip()
                if not item:
                    continue
                if "spp" in item:
                    entry["any"] = True
                    continue
                parts = item.split()
                # Expand abbreviations like "F. prausnitzii" to "Faecalibacterium prausnitzii"
                if parts[0].endswith("."):
                    entry["species"].add(f"{genus} {' '.join(parts[1:])}")
                else:
                    entry["species"].add(item)
            if not entry["species"] and not entry["any"]:
                entry["any"] = True  # genus with no species listed = any species
            targets[genus] = entry
    return targets


def esearch(term, retmax):
    handle = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
    ids = Entrez.read(handle)["IdList"]
    handle.close()
    time.sleep(0.4)
    return ids


def efetch(ids):
    if not ids:
        return []
    handle = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="fasta", retmode="text")
    recs = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    time.sleep(0.4)
    return recs


def fetch_references(targets, out_fasta):
    """Download 16S sequences from NCBI for every target genus and species."""
    seen = set()
    count = 0
    with open(out_fasta, "w") as out:
        for genus, entry in targets.items():
            # Named species: search each species (and its alias names)
            for sp in sorted(entry["species"]):
                epithet = sp.split()[1]
                for name in ALIASES.get(sp, [sp]):
                    term = f'"{name}"[Organism] AND 16S ribosomal RNA[Title]'
                    for rec in efetch(esearch(term, REFS_PER_SPECIES)):
                        if not (MIN_LEN <= len(rec.seq) <= MAX_LEN):
                            continue
                        if "uncultured" in rec.description.lower():
                            continue
                        acc = rec.id.split(".")[0]
                        if acc in seen:
                            continue
                        seen.add(acc)
                        rec.id = f"{genus}__{epithet}__{acc}"
                        rec.description = ""
                        SeqIO.write(rec, out, "fasta")
                        count += 1
            # Whole-genus entries: download several species from this genus
            if entry["any"]:
                term = f'"{genus}"[Organism] AND 16S ribosomal RNA[Title]'
                for rec in efetch(esearch(term, MAX_REFS_PER_GENUS * 3)):
                    toks = rec.description.split()
                    if len(toks) < 3 or toks[1] != genus or toks[2].endswith("."):
                        continue
                    if not (MIN_LEN <= len(rec.seq) <= MAX_LEN):
                        continue
                    if "uncultured" in rec.description.lower():
                        continue
                    acc = rec.id.split(".")[0]
                    if acc in seen:
                        continue
                    seen.add(acc)
                    rec.id = f"{genus}__{toks[2]}__{acc}"
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    count += 1
            print(f"Fetched references for {genus}")
    print(f"Total reference sequences: {count}")


def build_db(ref_fasta):
    subprocess.run(
        ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", DB_NAME],
        check=True, capture_output=True,
    )


def run_blast(sample_fasta):
    """Return the best reference hit for each sample sequence."""
    cmd = [
        "blastn", "-query", sample_fasta, "-db", DB_NAME,
        "-outfmt", "6 qseqid sseqid pident length qlen bitscore",
        "-evalue", str(EVALUE_MAX), "-max_hsps", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    best = {}
    for line in result.stdout.splitlines():
        qid, sid, pident, length, qlen, bits = line.split("\t")
        genus, species, acc = sid.split("__")
        hit = {
            "genus": genus,
            "species": species,
            "identity": float(pident),
            "coverage": 100 * int(length) / int(qlen),
            "bitscore": float(bits),
        }
        if qid not in best or hit["bitscore"] > best[qid]["bitscore"]:
            best[qid] = hit
    return best


def classify(hit, targets):
    """Return (genus label, species call). Label is 'Unknown' if no rule matches."""
    if not hit or hit["coverage"] < COVERAGE_MIN or hit["genus"] not in targets:
        return UNKNOWN, ""
    entry = targets[hit["genus"]]
    full = f"{hit['genus']} {hit['species']}"
    if full in entry["species"] and hit["identity"] >= SPECIES_MIN_ID:
        return hit["genus"], full
    if entry["any"] and hit["identity"] >= GENUS_MIN_ID:
        return hit["genus"], ""
    return UNKNOWN, ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True, help="targets.txt file")
    parser.add_argument("--sample", required=True, help="stool sample in FASTA format")
    parser.add_argument("--out", default="results.csv")
    parser.add_argument("--skip-fetch", action="store_true", help="reuse existing references.fasta")
    args = parser.parse_args()

    targets = read_targets(args.targets)
    if not targets:
        sys.exit("No targets found in the targets file.")

    if not args.skip_fetch:
        fetch_references(targets, REF_FASTA)
    build_db(REF_FASTA)

    hits = run_blast(args.sample)
    sample_ids = [rec.id for rec in SeqIO.parse(args.sample, "fasta")]

    counts = Counter()
    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "label", "species_call", "best_match",
                         "identity_pct", "coverage_pct"])
        for qid in sample_ids:
            hit = hits.get(qid)
            label, species_call = classify(hit, targets)
            counts[label] += 1
            if hit:
                best = f"{hit['genus']} {hit['species']}"
                writer.writerow([qid, label, species_call, best,
                                 f"{hit['identity']:.2f}", f"{hit['coverage']:.1f}"])
            else:
                writer.writerow([qid, label, "", "", "", ""])

    print(f"\nResults saved to {args.out}")
    print(f"Total sequences: {len(sample_ids)}")
    for label, n in counts.most_common():
        print(f"  {label}: {n} ({100 * n / len(sample_ids):.1f}%)")


if __name__ == "__main__":
    main()