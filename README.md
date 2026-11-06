# covid-blast-analysis-
SARS-CoV-2 spike mutation profiling with BLAST and a validated Python caller
# SARS-CoV-2 Spike Mutation Profiling: BLAST + a Validated Mutation Caller

Identify what a SARS-CoV-2 genome is, find every amino acid change in its spike protein relative to the Wuhan-Hu-1 reference, and check the result against an independent method.

**Headline result:** on a Omicron BA.1 genome (GenBank [QB226460.1](https://www.ncbi.nlm.nih.gov/nuccore/QB226460.1)), the custom caller found **34 spike changes** (30 substitutions, 3 deletions, 1 insertion). All 34 match both a manual reading of the BLASTN alignment and an independent translation-based comparison (**34/34, 0 missing, 0 extra**).

---

## Why this project

BLAST tells you that two sequences are similar and where the nucleotides differ. It does not tell you which *protein* changes those differences cause, which is what matters for biology. This project closes that gap for the spike gene: from raw genome to a domain-annotated list of amino acid changes, with each step checked.

## Workflow

```
NCBI reference genome (NC_045512.2)         NCBI Virus: complete BA.1 genome (QB226460.1)
          |                                                |
          v                                                v
 extract spike CDS (21563-25384)                 Galaxy: Fasta Statistics (QC)
          |                                                |
          +----------------+-------------------------------+
                           v
            Galaxy: BLASTN (spike ref vs variant genome)
            -> confirms identity, shows nucleotide mismatches
                           |
                           v
            spike_mutation_caller.py
            local alignment -> translate -> protein alignment -> call variants
                           |
                           v
            results/*.tsv  (domain-annotated mutations)
                           |
                           v
            independent_spike_check.py  (stdlib-only validation)
```

## Data

| File | Description |
|---|---|
| `data/NC_045512.2.fasta` | SARS-CoV-2 reference genome, Wuhan-Hu-1 (29,903 nt) |
| `data/Spike_ref.fasta` | Reference spike CDS, NC_045512.2:21563-25384 (3,822 nt) |
| `data/sample_omicron.fasta` | QB226460.1, SARS-CoV-2/human/USA/NC-NCAT-Haley-20220111/2022 (29,851 nt) |

The sample was selected in NCBI Virus with filters for complete nucleotide sequence and Pango lineage BA.1. Its mutation profile is consistent with BA.1 (see results), but lineage was not independently confirmed with Pangolin or Nextclade (see Limitations).

## Methods

### 1. Quality control (Galaxy, Fasta Statistics)

FastQC is designed for sequencing reads, not assembled FASTA, so base-composition statistics were used instead. The reference spike:

| Metric | Value |
|---|---|
| Sequences | 1 |
| Length | 3,822 nt |
| GC content | 37.31% |
| A / T / C / G | 1125 / 1271 / 723 / 703 |
| Ambiguous bases (N) | 0 |
| Gaps | 0 |

The variant genome contains 51 ambiguous bases (N or IUPAC codes), which matters for downstream calling (see Notes).

### 2. BLASTN (Galaxy, NCBI BLAST+ 2.16.0+)

The reference spike was used as the query against the variant genome as a one-sequence subject, using the `blastn` task and an E-value cutoff of 0.001. The standard megablast task was avoided because it needs long exact seeds and tends to fragment alignments across dense mutation clusters.

| Metric | Value |
|---|---|
| Score | 6653 bits (7378) |
| E-value | 0.0 |
| Identities | 3772/3831 (98%) |
| Gaps | 27/3831 (1%) |
| Strand | Plus/Plus |
| Subject range | 21529-25341 |

The ~98% identity and 27 gap positions are what a heavily mutated Omicron spike looks like. The subject region is 9 nt shorter than the reference spike, which fits the net effect of the deletions and insertion reported below.

### 3. Spike mutation caller (`spike_mutation_caller.py`)

BLAST gives nucleotide differences. The caller converts them to protein-level changes:

1. Extract the spike CDS from the reference genome by coordinate (no hardcoded sequence).
2. Locate the spike in each query genome with a local nucleotide alignment on both strands, extended to the full reference span to preserve the reading frame.
3. Translate and globally align the protein to the 1273 aa reference (BLOSUM62, gap open -10, extend -0.5, free end gaps).
4. Walk the alignment to call substitutions (`D614G`), deletions (`del69-70`), insertions, and premature stops.
5. Annotate each change with its spike region (NTD, RBD, RBD/RBM, SD1/2, furin site, S2).
6. Flag possible frameshifts, reverse-complemented inputs, and high N content. Codons containing ambiguity codes translate to `X` and are **skipped rather than called**.

**Usage**

```bash
pip install biopython
python spike_mutation_caller.py -r data/NC_045512.2.fasta -q data/sample_omicron.fasta -o results/
```

The reference must be the **full genome**, not the spike-only file. The script asserts that the extracted spike translates to 1273 aa and stops otherwise. Multi-record query FASTA files are supported.

**Outputs:** `results/spike_mutations_long.tsv` (one row per mutation) and `results/spike_mutations_summary.tsv` (one row per sample).

### 4. Independent validation (`independent_spike_check.py`)

To avoid trusting the caller on its own output, a separate script re-derives the answer with a different approach and no third-party libraries: it locates the spike by exact end anchors, translates it with its own codon table, and diffs the protein with Python's `difflib`. It produced the same 34 events. The caller's output was then compared against the expected set in a Colab notebook: **34/34 matched, 0 missing, 0 extra.**

## Results

### Mutations found (QB226460.1 vs NC_045512.2)

| Region | Count | Mutations |
|---|---|---|
| NTD (14-305) | 8 | A67V, del69-70, T95I, G142D, del143-145, del211, L212I, ins214 |
| RBD (319-541, outside RBM) | 5 | G339D, S371L, S373P, S375F, K417N |
| RBD / RBM (437-508) | 10 | N440K, G446S, S477N, T478K, E484A, Q493R, G496S, Q498R, N501Y, Y505H |
| SD1/SD2 (542-685) | 4 | T547K, D614G, H655Y, N679K |
| Furin site (681-685) | 1 | P681H |
| S2 (686-1273) | 6 | N764K, D796Y, N856K, Q954H, N969K, L981F |
| **Total** | **34** | 30 substitutions, 3 deletions, 1 insertion |

The C-terminal ~350 nt of the spike are identical to the reference, which the BLAST alignment confirms as an unbroken run of matches through the `TAA` stop codon.

### Reading the alignment by hand

Three changes were traced from the BLAST alignment to the protein level as a cross-check on the caller:

- **D614G**: `GAT` to `GGT` at spike nucleotides 1840-1842, a single A-to-G change. Spike position 614 sits in the SD2 region, and this change became near-universal early in the pandemic and has been associated with higher infectivity and a more stable spike.
- **N501Y**: `AAT` to `TAT` at nucleotides ~1501-1503. Position 501 lies in the receptor-binding motif at the ACE2 interface, and this change has been linked to increased ACE2 binding affinity.
- **P681H**: `CCT` to `CAT`, adjacent to the S1/S2 furin cleavage site, a region linked to spike processing.

Other changes worth noting: `K417N`, `E484A` and `Q493R` fall in antibody-targeted parts of the RBD and are commonly discussed in the context of immune escape, and `del69-70` is the deletion behind S-gene target failure in some PCR assays. These functional statements summarize published literature and are not conclusions drawn from this dataset.

## Notes on the data

- **Insertion at 214.** The genome has masked bases (`N`) where the known BA.1 `ins214EPE` insertion sits, so the caller reports `ins214_215XXX`. The *position and length* are called correctly, but the amino acids cannot be resolved from this assembly.
- **Ambiguity codes.** Two codons contain IUPAC ambiguity codes: a silent one at codon 11 and one at codon 346 (A/G, which would encode K or R). These are skipped, not reported as mutations, since the base call is uncertain. A reported `R346K` here would be unsupported.

## Limitations

- **One validation genome.** Validated on a single BA.1 sample. Other lineages (Delta, XBB, recombinants) are untested.
- **No lineage calling.** Lineage assignment is outside the scope of this tool and BLAST. Results were not benchmarked against Pangolin or Nextclade; doing so is the first planned extension.
- **Ambiguity handling is conservative.** Codons with ambiguity codes are skipped, so mutations hidden behind Ns will be missed rather than guessed.
- **Equivalent indel representations.** Adjacent indels and substitutions (for example `del143-145` beside `G142D`) can be written more than one way. The caller's output matched the expected set here, but other genomes may need normalization before comparing to other tools.
- **Speed.** Whole-genome local alignment is fine for tens to hundreds of genomes but is not designed for millions. A seed-and-extend or reference-anchored approach would scale better.
- **Reverse-complement support is untested** on real data.

## Reproducing this

**Colab or local**

```bash
git clone https://github.com/<Oyedamola-bioinfomatics>/covid-blast-analysis.git
cd covid-blast-analysis
pip install biopython pandas
python spike_mutation_caller.py -r data/NC_045512.2.fasta -q data/sample_omicron.fasta -o results/
python independent_spike_check.py
```

**Galaxy steps** (usegalaxy.org): upload `Spike_ref.fasta` and `sample_omicron.fasta` as type `fasta`, run **Fasta Statistics** for QC, then **NCBI BLAST+ blastn** with the spike as the query and the genome as the subject (FASTA file from history), `blastn` task, Pairwise text output.

## Repository layout

```
covid-blast-analysis/
├── README.md
├── spike_mutation_caller.py
├── independent_spike_check.py
├── data/
│   ├── NC_045512.2.fasta
│   ├── Spike_ref.fasta
│   └── sample_omicron.fasta
├── results/
│   ├── spike_mutations_long.tsv
│   └── spike_mutations_summary.tsv
└── screenshots/
    ├── galaxy_fasta_stats.png
    └── blast_alignment.png
```

## Roadmap

- [ ] Benchmark against Pangolin and Nextclade on 20+ genomes across lineages
- [ ] Add a lineage-defining mutation lookup to infer the variant from the calls
- [ ] Mutation heatmap across samples, grouped by spike domain
- [ ] Wrap the pipeline in Snakemake with a conda `environment.yml`
- [ ] Normalize indel notation for direct comparison with other callers

## What I learned

- BLAST output is nucleotide-level; mapping it to amino acid changes needs translation and a protein alignment, and that step is where most of the care goes.
- Real genomes are messy: masked bases and IUPAC codes force a decision about whether to call, skip, or flag, and that decision should be explicit and documented.
- A caller should be validated against something that does not share its code. Checking the same logic twice proves nothing, so the independent script and the hand-read alignment matter.
- Choosing the right tool matters: FastQC is for reads, and megablast is the wrong BLAST task for dense Omicron-style mutation clusters.

## Tools

NCBI Nucleotide and NCBI Virus (data), Galaxy (Fasta Statistics, BLAST+ 2.16.0+), Python 3 with Biopython (PairwiseAligner), pandas, Google Colab.

## References and data sources

- Reference genome: NC_045512.2, SARS-CoV-2 isolate Wuhan-Hu-1 (NCBI RefSeq)
- Variant genome: QB226460.1 (GenBank)
- Altschul et al. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research* 25:3389-3402.

## Author

[Oyelakin Damola Ayomide] | Bioinformatics | [https://www.linkedin.com/in/oyelakin-damola-ayomide]
