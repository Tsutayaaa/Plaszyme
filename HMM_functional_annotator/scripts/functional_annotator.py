import argparse
import subprocess
import sys
import shutil
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd

# Utility functions

def run(cmd: str, cwd: Optional[Path] = None):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=str(cwd) if cwd else None)

def ensuredir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def check_bin(name: str, suggest: str):
    if shutil.which(name) is None:
        sys.exit(f"[ERROR] `{name}` not found in PATH. {suggest}")

def fasta_count_headers(fasta: Path) -> int:
    n = 0
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n

# MAFFT + HMMER Build

def align_fasta(in_fasta: Path, out_aln: Path, threads: int = 8, err_dir: Path = Path("logs")):
    ensuredir(err_dir)
    err_log = err_dir / f"mafft_{in_fasta.stem}.err"
    cmd = ["mafft", "--thread", str(threads), "--anysymbol", "--auto", str(in_fasta)]
    print("$", " ".join(cmd), f"> {out_aln}")
    with open(out_aln, "w") as fo:
        p = subprocess.run(cmd, stdout=fo, stderr=subprocess.PIPE, text=True)
    err_log.write_text(p.stderr or "")
    if p.returncode != 0:
        try:
            out_aln.unlink()
        except FileNotFoundError:
            pass
        raise RuntimeError(
            f"[MAFFT ERROR] returncode={p.returncode}\n"
            f"Input: {in_fasta}\nOutput: {out_aln}\nSee: {err_log}\n---- STDERR ----\n{p.stderr}"
        )
    
def hmmbuild_from_aln(aligned_fasta: Path, out_hmm: Path):
    ensuredir(out_hmm.parent)
    run(f"hmmbuild {out_hmm} {aligned_fasta}")

def combine_hmms(hmm_files: List[Path], out_combined: Path, run_press: bool = True):
    ensuredir(out_combined.parent)
    with out_combined.open("w") as out:
        for hmm in hmm_files:
            out.write(Path(hmm).read_text())
    if run_press:
        run(f"hmmpress {out_combined}")

# HMMER Scan & Search

def hmmscan_unknown(hmm_db: Path, unknown_fasta: Path, out_dir: Path, threads: int = 8):
    ensuredir(out_dir)
    tbl = out_dir / "hmmscan_tblout.tsv"
    dom = out_dir / "hmmscan_domtblout.tsv"
    log = out_dir / "hmmscan_alignment.txt"
    cmd = (
        f"hmmscan --cpu {threads} "
        f"--tblout {tbl} "
        f"--domtblout {dom} "
        f"{hmm_db} {unknown_fasta} > {log}"
    )
    run(cmd)
    return tbl, dom, log

def hmmsearch_one(hmm: Path, seq_db: Path, out_dir: Path, threads: int = 8):
    ensuredir(out_dir)
    base = hmm.stem
    tbl = out_dir / f"{base}.tblout.tsv"
    dom = out_dir / f"{base}.domtblout.tsv"
    log = out_dir / f"{base}.alignment.txt"
    cmd = (
        f"hmmsearch --cpu {threads} "
        f"--tblout {tbl} --domtblout {dom} "
        f"{hmm} {seq_db} > {log}"
    )
    run(cmd)
    return tbl, dom, log

# Parsing and filtering

def parse_domtblout(domtbl: Path) -> pd.DataFrame:
    rows = []
    with domtbl.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            target, tlen = parts[0], int(parts[2])
            query, qlen  = parts[3], int(parts[5])
            full_E, full_score = float(parts[6]), float(parts[7])
            dom_num, dom_of = int(parts[9]), int(parts[10])
            c_E, i_E, dom_score = float(parts[11]), float(parts[12]), float(parts[13])
            hmm_from, hmm_to = int(parts[15]), int(parts[16])
            ali_from, ali_to = int(parts[17]), int(parts[18])
            env_from, env_to = int(parts[19]), int(parts[20])
            acc = float(parts[21])
            rows.append({
                "hmm": target, "tlen": tlen, "query": query, "qlen": qlen,
                "full_E": full_E, "full_score": full_score,
                "c_E": c_E, "i_E": i_E, "dom_score": dom_score,
                "hmm_from": hmm_from, "hmm_to": hmm_to,
                "ali_from": ali_from, "ali_to": ali_to,
                "env_from": env_from, "env_to": env_to,
                "acc": acc, "dom_index": dom_num, "dom_of": dom_of
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["qcov"] = (df["ali_to"] - df["ali_from"] + 1) / df["qlen"]
    df["hcov"] = (df["hmm_to"] - df["hmm_from"] + 1) / df["tlen"]
    return df

def filter_hits(df: pd.DataFrame, iE: float = 1e-5, bits: float = 25.0, qcov: float = 0.3) -> pd.DataFrame:
    if df.empty:
        return df
    m = (df["i_E"] <= iE) & (df["dom_score"] >= bits) & (df["qcov"] >= qcov)
    return df.loc[m].copy()

def best_hit_per_query(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.sort_values(["query", "i_E", "dom_score"], ascending=[True, True, False])
    return df.groupby("query", as_index=False).first()

def discover_families(families_dir: Path) -> Dict[str, Path]:
    fams = {}
    for ext in ("*.fa", "*.faa", "*.fasta"):
        for p in families_dir.glob(ext):
            fams[p.stem] = p
    if not fams:
        sys.exit(f"[ERROR] No FASTA files found in {families_dir}")
    return fams

def build_family_hmms_separate_dirs(
    families: Dict[str, Path],
    aln_dir: Path,
    hmm_dir: Path,
    threads: int = 8
) -> List[Path]:
    ensuredir(aln_dir)
    ensuredir(hmm_dir)
    hmm_files = []
    for fam, fasta in families.items():
        if fasta_count_headers(fasta) < 2:
            print(f"[WARN] Family {fam}: <2 sequences, skip HMM build.")
            continue
        aln = aln_dir / f"{fam}.aln.fasta"
        hmm = hmm_dir / f"{fam}.hmm"
        print(f"[INFO] Building HMM for family: {fam}")
        align_fasta(fasta, aln, threads=threads)
        hmmbuild_from_aln(aln, hmm)
        hmm_files.append(hmm)
    if not hmm_files:
        sys.exit("[ERROR] No HMMs were created. Check family FASTAs.")
    return hmm_files

# Command-line interface

def cmd_build(args):
    families = discover_families(Path(args.families_dir))
    aln_dir = Path(args.aln_dir)
    hmm_dir = Path(args.hmm_dir)
    ensuredir(aln_dir)
    ensuredir(hmm_dir)
    hmm_files = build_family_hmms_separate_dirs(
        families,
        aln_dir=aln_dir,
        hmm_dir=hmm_dir,
        threads=args.threads
    )
    print(f"[OK] Built {len(hmm_files)} HMMs into {hmm_dir} (MSA in {aln_dir}). No combined.hmm generated.")

def cmd_annotate(args):
    unknown_fasta = Path(args.unknown_fasta)
    hmm_dir = Path(args.hmm_dir)
    out_dir = Path(args.out_dir)
    threads = args.threads
    ensuredir(out_dir)

    hmm_list = sorted(hmm_dir.glob("*.hmm")) 
    if not hmm_list: 
        sys.exit(f"[ERROR] No .hmm files found in {hmm_dir}")

    dom_paths = []
    for hmm in hmm_list:
        print(f"[INFO] hmmsearch: {hmm.name} vs {unknown_fasta.name}")
        sub_dir = out_dir / hmm.stem
        ensuredir(sub_dir)
        tbl = sub_dir / f"{hmm.stem}.tblout.tsv"
        dom = sub_dir / f"{hmm.stem}.domtblout.tsv"
        log = sub_dir / f"{hmm.stem}.alignment.txt"
        cmd = (
            f"hmmsearch --cpu {threads} "
            f"--tblout {tbl} --domtblout {dom} "
            f"{hmm} {unknown_fasta} > {log}"
            )
        run(cmd)
        dom_paths.append((hmm.name, dom))

    frames = []
    for hmm_name, p in dom_paths:
        df = parse_domtblout(p)
        if not df.empty:
            df["hmm_file"] = hmm_name
            frames.append(df)
            
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    
    filt = filter_hits(merged, iE=args.iE, bits=args.bits, qcov=args.qcov)
    if not filt.empty:
        for hmm_name, group in filt.groupby("hmm_file"):
            sub_dir = out_dir / Path(hmm_name).stem
            ensuredir(sub_dir)
            group.to_csv(sub_dir / f"{Path(hmm_name).stem}__filtered.tsv",
                        sep="\t", index=False)
            best_one = best_hit_per_query(group)
            best_one.to_csv(sub_dir / f"{Path(hmm_name).stem}__best.tsv",
                            sep="\t", index=False)
    best = best_hit_per_query(filt)
    
    filt.to_csv(out_dir / "annotation_report__dom_hits.tsv", sep="\t", index=False)
    best.to_csv(out_dir / "annotation_report__best.tsv", sep="\t", index=False)
    print(f"Annotation finished. TSV results saved in: {out_dir}")


def cmd_searchdb(args):
    seq_db = Path(args.seq_db)
    out_dir = Path(args.out_dir)
    ensuredir(out_dir)

    hmm_list: List[Path] = []
    if args.hmm:
        hmm_list.append(Path(args.hmm))
    if args.hmm_dir:
        hmm_list.extend(sorted(Path(args.hmm_dir).glob("*.hmm")))
    if not hmm_list:
        sys.exit("[ERROR] Please provide --hmm or --hmm_dir with .hmm files.")

    dom_paths: List[Tuple[str, Path, Path]] = [] 
    for hmm in hmm_list:
        hmm_name = hmm.stem                      
        sub_dir = out_dir / hmm_name
        ensuredir(sub_dir)
        print(f"[INFO] Searching DB with HMM: {hmm.name} -> {sub_dir}")

        _, dom, _ = hmmsearch_one(hmm, seq_db, sub_dir, threads=args.threads)
        dom_paths.append((hmm_name, dom, sub_dir))

    frames: List[pd.DataFrame] = []
    for hmm_name, dom_path, sub_dir in dom_paths:  
        df = parse_domtblout(dom_path)
        if not df.empty:
            df["hmm_file"] = hmm_name
            frames.append(df)
            df.to_csv(sub_dir / f"{hmm_name}__all_hits_enhanced.tsv",
                      sep="\t", index=False)

    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    filt = filter_hits(merged, iE=args.iE, bits=args.bits, qcov=args.qcov) if not merged.empty else merged
    best = best_hit_per_query(filt) if not filt.empty else filt

    if not filt.empty:
        for hmm_name, group in filt.groupby("hmm_file"):
            sub_dir = out_dir / hmm_name
            ensuredir(sub_dir)
            group.to_csv(sub_dir / f"{hmm_name}__filtered.tsv", sep="\t", index=False)
            best_one = best_hit_per_query(group)
            best_one.to_csv(sub_dir / f"{hmm_name}__best.tsv", sep="\t", index=False)

    filt.to_csv(out_dir / "searchdb_report__dom_hits.tsv", sep="\t", index=False)
    best.to_csv(out_dir / "searchdb_report__best.tsv", sep="\t", index=False)
    print(f"[OK] DB search finished. TSV results in: {out_dir}")


# Main entry point

def main():
    # HMMER + MAFFT
    check_bin("mafft", "Install MAFFT: conda install -c bioconda mafft")
    check_bin("hmmbuild", "Install HMMER: conda install -c bioconda hmmer")
    check_bin("hmmpress", "Install HMMER: conda install -c bioconda hmmer")
    check_bin("hmmscan", "Install HMMER: conda install -c bioconda hmmer")
    check_bin("hmmsearch", "Install HMMER: conda install -c bioconda hmmer")

    p = argparse.ArgumentParser(
        prog="functional_annotator",
        description="BLAST + HMMER functional annotation pipeline"
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # build
    p_build = sub.add_parser("build", help="Build per-family HMMs; save MSA and HMM to separate dirs (no combined.hmm)")
    p_build.add_argument("--families_dir", required=True)
    p_build.add_argument("--aln_dir", required=True, help="Where to write *.aln.fasta")
    p_build.add_argument("--hmm_dir", required=True, help="Where to write *.hmm")
    p_build.add_argument("--threads", type=int, default=8)
    p_build.set_defaults(func=cmd_build)

    # annotate 
    p_ann = sub.add_parser("annotate", help="Annotate unknown sequences by looping hmmsearch over all HMMs in a directory")
    p_ann.add_argument("--unknown_fasta", required=True)
    p_ann.add_argument("--hmm_dir", required=True, help="Directory containing *.hmm")
    p_ann.add_argument("--out_dir", required=True)
    p_ann.add_argument("--threads", type=int, default=8)
    p_ann.add_argument("--iE", type=float, default=1e-5)
    p_ann.add_argument("--bits", type=float, default=25.0)
    p_ann.add_argument("--qcov", type=float, default=0.30)
    p_ann.set_defaults(func=cmd_annotate)

    # searchdb
    p_sdb = sub.add_parser("searchdb", help="Search sequence DB with HMM(s) (hmmsearch)")
    p_sdb.add_argument("--seq_db", required=True, help="Downloaded protein DB in FASTA")
    p_sdb.add_argument("--out_dir", required=True, help="Output folder")
    p_sdb.add_argument("--hmm", help="A single HMM file")
    p_sdb.add_argument("--hmm_dir", help="Directory containing *.hmm (batch mode)")
    p_sdb.add_argument("--threads", type=int, default=8)
    p_sdb.add_argument("--iE", type=float, default=1e-5)
    p_sdb.add_argument("--bits", type=float, default=25.0)
    p_sdb.add_argument("--qcov", type=float, default=0.30)
    p_sdb.set_defaults(func=cmd_searchdb)

    args = p.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
