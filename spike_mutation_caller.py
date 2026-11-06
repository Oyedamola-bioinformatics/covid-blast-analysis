#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

SPIKE_START, SPIKE_END = 21563, 25384

DOMAINS = [
    ("SP", 1, 13), ("NTD", 14, 305), ("RBD", 319, 541),
    ("SD1/2", 542, 685), ("S2", 686, 1273),
]
RBM = (437, 508)
FURIN = (681, 685)


def domain_of(pos):
    if FURIN[0] <= pos <= FURIN[1]:
        return "furin_site"
    if RBM[0] <= pos <= RBM[1]:
        return "RBD/RBM"
    for name, s, e in DOMAINS:
        if s <= pos <= e:
            return name
    return "linker"


def nt_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -6, -2
    return a


def protein_aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5
    a.target_end_gap_score = 0
    a.query_end_gap_score = 0
    return a


def extract_spike(ref_cds, genome, aligner):
    best = None
    for strand, seq in (("+", genome), ("-", genome.reverse_complement())):
        aln = aligner.align(seq, ref_cds)[0]
        if best is None or aln.score > best[0].score:
            best = (aln, seq, strand)
    aln, seq, strand = best
    g_blocks, r_blocks = aln.aligned
    g_start, g_end = g_blocks[0][0], g_blocks[-1][1]
    r_start, r_end = r_blocks[0][0], r_blocks[-1][1]
    start = max(0, g_start - r_start)
    end = min(len(seq), g_end + (len(ref_cds) - r_end))
    return seq[start:end], strand


def translate(seq):
    seq = seq[: len(seq) - len(seq) % 3]
    return seq.translate()


def call_variants(ref_aa, qry_aa, aligner):
    aln = aligner.align(ref_aa, qry_aa)[0]
    r_row, q_row = str(aln[0]), str(aln[1])
    variants = []
    pos = 0
    i, n = 0, len(r_row)
    while i < n:
        r, q = r_row[i], q_row[i]
        if r != "-" and q != "-":
            pos += 1
            if r != q and q not in "X*" and r != "X":
                variants.append(dict(type="SNV", pos=pos, label=f"{r}{pos}{q}"))
            elif q == "*":
                variants.append(dict(type="STOP", pos=pos, label=f"{r}{pos}*"))
            i += 1
        elif q == "-":
            start = pos + 1
            j = i
            while j < n and q_row[j] == "-" and r_row[j] != "-":
                pos += 1
                j += 1
            label = f"del{start}" if pos == start else f"del{start}-{pos}"
            variants.append(dict(type="DEL", pos=start, label=label))
            i = j
        else:
            j = i
            while j < n and r_row[j] == "-":
                j += 1
            ins = q_row[i:j]
            variants.append(dict(type="INS", pos=pos, label=f"ins{pos}_{pos + 1}{ins}"))
            i = j
    return variants


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-r", "--reference", required=True)
    p.add_argument("-q", "--query", required=True)
    p.add_argument("-o", "--outdir", default="results")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_genome = next(SeqIO.parse(args.reference, "fasta")).seq.upper()
    ref_cds = ref_genome[SPIKE_START - 1 : SPIKE_END]
    ref_aa = str(translate(ref_cds)).rstrip("*")
    assert len(ref_aa) == 1273, f"Expected 1273 aa spike, got {len(ref_aa)}"

    nt_al, aa_al = nt_aligner(), protein_aligner()
    long_rows, summary_rows = [], []

    for rec in SeqIO.parse(args.query, "fasta"):
        genome = rec.seq.upper()
        spike_nt, strand = extract_spike(ref_cds, genome, nt_al)
        flags = []
        if strand == "-":
            flags.append("reverse_complemented")
        if len(spike_nt) % 3:
            flags.append("possible_frameshift")
        if spike_nt.count("N") / max(len(spike_nt), 1) > 0.01:
            flags.append("high_N_content")

        qry_aa = str(translate(spike_nt)).rstrip("*")
        variants = call_variants(ref_aa, qry_aa, aa_al)

        for v in variants:
            long_rows.append(dict(sample=rec.id, mutation=v["label"], type=v["type"],
                                  position=v["pos"], domain=domain_of(v["pos"])))
        summary_rows.append(dict(sample=rec.id, n_mutations=len(variants),
                                 mutations=";".join(v["label"] for v in variants),
                                 spike_len_aa=len(qry_aa), flags=",".join(flags) or "none"))
        print(f"{rec.id}: {len(variants)} spike changes [{','.join(flags) or 'ok'}]")

    for name, rows in (("spike_mutations_long.tsv", long_rows),
                       ("spike_mutations_summary.tsv", summary_rows)):
        if rows:
            with open(outdir / name, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
    print(f"Wrote results to {outdir}/")


if __name__ == "__main__":
    main()
