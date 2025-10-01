# Functional Annotator (MAFFT + HMMER)

A lightweight, reproducible pipeline for **protein functional annotation** using **Multiple Sequence Alignment (MAFFT)** and **profile HMMs (HMMER)**.

---

## Features

- **Build per-family HMMs** from protein FASTA files (≥2 sequences per family), using `hmmbuild`.
- **Annotate unknown proteins** against a combined HMM database using `hmmscan`.
- **Search downloaded databases** (e.g., UniProt subsets) **with HMMs** using `hmmsearch`.
- Ready-to-run **Slurm scripts** for HPC clusters.

---

## Directory Layout

Recommended project layout (create at any path):

```
functional_annotator/
├─ scripts/
│  └─ functional_annotator.py         # main pipeline (MAFFT + HMMER)
├─ run/
   ├─ logs/                            # Slurm logs
│  ├─ run_build_hmm.sh                 # build per-family HMMs and combined DB
│  ├─ run_annotation.sh                # annotate unknown sequences with hmmscan
│  └─ run_searchdb.sh                  # search a downloaded DB with hmmsearch
├─ data/                               # raw inputs (UNaligned)
│  ├─ families/                        # one FASTA per family (≥2 sequences)
│  │  ├─ plastic_families
│  │  ├─ pfam_families
│  │  └─ clan_families
│  ├─ unknown_sequences                # sequences to annotate
   |   └─ proteins_converted.fasta
│  └─ database/                              # downloaded protein DBs for searchdb
│     └─ uniprot_swissprot_2025_08.fasta
├─ hmm_profiles/                       # outputs from build (auto)
├─ annotation_results/                 # outputs from annotate (auto)
└─ search_results/                     # outputs from searchdb (auto)
```


## Dependencies

Install with conda (recommended):

```bash
conda create -n hmm python=3.10 -y
conda activate hmm

# core tools
conda install -c bioconda hmmer mafft -y

# python libs for report
conda install -c conda-forge pandas openpyxl -y
```

Verify:
```bash
mafft --version
hmmbuild -h
hmmpress -h
hmmscan -h
hmmsearch -h
```

> **Input type:** All inputs are **protein FASTA** unless you have DNA-specific HMMs (not covered here).

---

## Usage (CLI)

> Change the path to your own path

### 1) Build per-family HMMs and combined DB

```bash
cd ../HMM_functional_annotator/run
sbatch run_build_hmm.sh

```

---

### 2) Annotate unknown proteins with the HMM DB

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation.sh

```

---

### 3) Search a downloaded protein database with your HMMs 

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation.sh

```

## Parameters & Tips

| Parameter | Meaning | Typical |
|---|---|---|
| `--threads` | CPU threads for MAFFT/HMMER | 8–32 |
| `--iE` | domain independent E-value cutoff | `1e-5` (strict) |
| `--bits` | domain bit-score cutoff | `20–30` |
| `--qcov` | query coverage threshold (0–1) | `0.30–0.50` |

- Ensure **protein** FASTA inputs for protein HMMs.
- Each family FASTA should have **≥2 sequences**; otherwise `hmmbuild` is not meaningful.

---
