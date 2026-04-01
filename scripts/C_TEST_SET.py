#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Merlin ORF starts (Merlin coordinates; same as your app)
UL97_START = 141797
UL54_START = 78193

# --------------------- helpers ---------------------
def which_or_die(cmd: str):
    if shutil.which(cmd) is None:
        sys.exit(f"❌ Required tool not found in PATH: {cmd}")

def read_fasta_any(path: Path) -> tuple[str, str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    header = None
    seqs = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:].split()[0]
                else:
                    seqs.append("N" * 100)
            else:
                seqs.append(line.upper())
    if header is None:
        raise ValueError(f"Not a FASTA? {path}")
    return header, "".join(seqs)

def write_fasta(path: Path, header: str, seq: str):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

def fasta_len(path: Path) -> int:
    _, seq = read_fasta_any(path)
    return len(seq)

def run_art(ref_fa: Path, out_prefix: Path, pairs: int, read_len: int,
            mean_frag: int, sd_frag: int, seq_sys: str, seed: int):
    cmd = [
        "art_illumina",
        "-ss", seq_sys,
        "-p",
        "-l", str(read_len),
        "-m", str(mean_frag),
        "-s", str(sd_frag),
        "-c", str(pairs),
        "-na",
        "-rs", str(seed),
        "-i", str(ref_fa),
        "-o", str(out_prefix),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def append_fastq(src1: Path, src2: Path, dst1: Path, dst2: Path):
    with open(dst1, "ab") as o1, open(dst2, "ab") as o2:
        with open(src1, "rb") as i1, open(src2, "rb") as i2:
            shutil.copyfileobj(i1, o1)
            shutil.copyfileobj(i2, o2)

def mb(path: Path) -> float:
    return path.stat().st_size / (1024 * 1024)

# -------------------- mutation helpers (real resistance edits) --------------------
CODON_TABLE = {
    "I": {"ATA", "ATT", "ATC"},
    "M": {"ATG"},
    "V": {"GTT", "GTC", "GTA", "GTG"},
    "L": {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"},
    "A": {"GCT", "GCC", "GCA", "GCG"},
    "P": {"CCT", "CCC", "CCA", "CCG"},
    "S": {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"},
}

def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)

def codon_at(seq: str, pos1: int) -> str:
    i = pos1 - 1
    return seq[i:i+3]

def replace_codon_to_aa(seq: str, codon_start_1based: int, target_aa: str) -> tuple[str, str, str]:
    old = codon_at(seq, codon_start_1based)
    if len(old) != 3 or any(b not in "ACGT" for b in old):
        new = sorted(CODON_TABLE[target_aa])[0]
    else:
        candidates = sorted(CODON_TABLE[target_aa])
        new = min(candidates, key=lambda c: hamming(old, c))
    i = codon_start_1based - 1
    new_seq = seq[:i] + new + seq[i+3:]
    return new_seq, old, new

def delete_aa_block(seq: str, start_1based: int, aa_count: int) -> tuple[str, str]:
    i = start_1based - 1
    n = aa_count * 3
    deleted = seq[i:i+n]
    new_seq = seq[:i] + seq[i+n:]
    return new_seq, deleted

# -------------------- reference classification --------------------
def is_human_ref(path: Path) -> bool:
    name = path.name.lower()
    return bool(re.search(r"\bchr\d+\b|chr[xy]\b|chrm\b", name)) or "human" in name or "hg" in name

def is_cmv_ref(path: Path) -> bool:
    name = path.name.lower()
    return ("merlin" in name) or ("ad169" in name) or ("towne" in name) or ("cmv" in name) or ("hcmv" in name)

def discover_fastas(dir_path: Path) -> list[Path]:
    refs = []
    for ext in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
        refs.extend(dir_path.glob(ext))
    return sorted(set(refs))

def logspace_fracs(min_f: float, max_f: float, n_levels: int):
    if min_f <= 0 or max_f <= 0 or min_f >= max_f:
        raise ValueError("min/max fractions must be >0 and min < max.")
    if n_levels == 1:
        return [math.sqrt(min_f * max_f)]
    a = math.log10(min_f)
    b = math.log10(max_f)
    step = (b - a) / (n_levels - 1)
    return [10 ** (a + i * step) for i in range(n_levels)]

def sample_lognormal_pairs(rng: random.Random, median: int, spread: float, min_pairs: int, max_pairs: int) -> int:
    x = rng.lognormvariate(mu=math.log(median), sigma=spread)
    x = max(min_pairs, min(max_pairs, int(x)))
    return x

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate CMV diagnostic GUI test FASTQs: human-only, other-virus, CMV WT, contamination, co-infection, and antiviral resistance mutants."
    )
    ap.add_argument("--refs-dir", required=True,
                    help="Directory containing FASTAs for: CMV (e.g., Merlin), other viruses, and optional human refs. (You can also pass --human-ref explicitly).")
    ap.add_argument("--outdir", required=True, help="Output folder (fastq/ + truth/metadata.tsv).")

    ap.add_argument("--human-ref", default="", help="Explicit human FASTA (e.g. chr21.fasta). If not set, script tries to find one in --refs-dir.")

    # Filters so you don't accidentally classify tiny files as CMV
    ap.add_argument("--min-cmv-len", type=int, default=150000, help="Min length to treat a FASTA as CMV genome (default 150k).")
    ap.add_argument("--min-other-len", type=int, default=5000, help="Min length to treat a FASTA as other-virus (default 5k).")

    # Library / depth knobs
    ap.add_argument("--median-pairs", type=int, default=600000)
    ap.add_argument("--pair-spread", type=float, default=0.6)
    ap.add_argument("--min-pairs", type=int, default=120000)
    ap.add_argument("--max-pairs", type=int, default=1400000)

    ap.add_argument("--min-fastq-mb", type=float, default=100.0, help="If R1 < this, pad with extra human reads (if human-ref available).")

    ap.add_argument("--read-lens", default="100,150")
    ap.add_argument("--seq-sys-list", default="HS25")
    ap.add_argument("--mean-frag-list", default="350,500")
    ap.add_argument("--sd-frag-list", default="30,50,80")

    # Counts per group
    ap.add_argument("--n-neg-human", type=int, default=8)
    ap.add_argument("--n-neg-othervirus", type=int, default=8)
    ap.add_argument("--n-neg-contam", type=int, default=8)

    ap.add_argument("--n-pos-wt", type=int, default=8)
    ap.add_argument("--n-pos-coinf", type=int, default=8)          # CMV + other-virus
    ap.add_argument("--n-pos-resistance", type=int, default=10)     # resistance CMV
    ap.add_argument("--n-pos-res_co", type=int, default=6)          # resistance CMV + other-virus

    # Fraction knobs
    ap.add_argument("--cmv-frac-min", type=float, default=1e-4, help="Min CMV fraction for CMV-positive groups.")
    ap.add_argument("--cmv-frac-max", type=float, default=0.25, help="Max CMV fraction for CMV-positive groups.")
    ap.add_argument("--other-frac-min", type=float, default=0.02, help="Min other-virus fraction for other-virus groups.")
    ap.add_argument("--other-frac-max", type=float, default=0.30, help="Max other-virus fraction for other-virus groups.")
    ap.add_argument("--neg-contam-min-frac", type=float, default=1e-7)
    ap.add_argument("--neg-contam-max-frac", type=float, default=5e-5)

    # Resistance set (cycled)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    which_or_die("art_illumina")
    rng = random.Random(args.seed)

    refs_dir = Path(args.refs_dir).resolve()
    outdir = Path(args.outdir).resolve()
    out_fastq = outdir / "fastq"
    out_truth = outdir / "truth"
    out_fastq.mkdir(parents=True, exist_ok=True)
    out_truth.mkdir(parents=True, exist_ok=True)

    refs = discover_fastas(refs_dir)
    if not refs:
        sys.exit(f"❌ No FASTA found under: {refs_dir}")

    # Find candidates
    cmv_candidates: list[Path] = []
    other_candidates: list[Path] = []
    human_candidates: list[Path] = []

    for p in refs:
        try:
            L = fasta_len(p)
        except Exception:
            continue
        if is_human_ref(p):
            human_candidates.append(p)
        elif is_cmv_ref(p) and L >= args.min_cmv_len:
            cmv_candidates.append(p)
        elif (not is_cmv_ref(p)) and (not is_human_ref(p)) and L >= args.min_other_len:
            other_candidates.append(p)

    if not cmv_candidates:
        sys.exit("❌ No CMV references found. Put Merlin/AD169/Towne (full genome) FASTAs in --refs-dir (name must include merlin/ad169/towne/cmv/hcmv).")

    # Pick human ref
    human_ref: Path | None = None
    if args.human_ref:
        human_ref = Path(args.human_ref).resolve()
        if not human_ref.exists():
            sys.exit(f"❌ --human-ref not found: {human_ref}")
    else:
        human_ref = human_candidates[0] if human_candidates else None

    if human_ref is None:
        print("⚠️ No human reference found (or provided). Samples will be viral-only where applicable.", file=sys.stderr)

    if not other_candidates:
        print("⚠️ No OTHER-virus FASTAs detected in --refs-dir (after filtering). Your 'other-virus' groups will end up human-only (or empty if no human).", file=sys.stderr)

    # ART parameter choices
    read_lens = [int(x) for x in args.read_lens.split(",") if x.strip()]
    seq_sys_list = [x.strip() for x in args.seq_sys_list.split(",") if x.strip()]
    mean_frags = [int(x) for x in args.mean_frag_list.split(",") if x.strip()]
    sd_frags   = [int(x) for x in args.sd_frag_list.split(",") if x.strip()]

    # Resistance mutation “modes” (cycled)
    # (mode_name, gene, description)
    resistance_modes = [
        ("UL97_M460I", "UL97", "UL97 M460I"),
        ("UL97_A594V", "UL97", "UL97 A594V"),
        ("UL97_L595S", "UL97", "UL97 L595S"),
        ("UL54_V812L", "UL54", "UL54 V812L"),
        ("UL54_A834P", "UL54", "UL54 A834P"),
        ("UL54_DEL981_982", "UL54", "UL54 Δ981–982"),
    ]

    # Plan samples
    plan: list[tuple[str, str, str]] = []
    for i in range(1, args.n_neg_human + 1):
        plan.append(("NEG_HUMAN_ONLY", f"NEG_H{i:03d}", ""))

    for i in range(1, args.n_neg_othervirus + 1):
        plan.append(("NEG_OTHER_VIRUS", f"NEG_O{i:03d}", ""))

    for i in range(1, args.n_neg_contam + 1):
        plan.append(("NEG_CONTAM", f"NEG_C{i:03d}", ""))

    for i in range(1, args.n_pos_wt + 1):
        plan.append(("POS_CMV_WT", f"POS_WT{i:03d}", ""))

    for i in range(1, args.n_pos_coinf + 1):
        plan.append(("POS_CMV_PLUS_OTHER", f"COINF{i:03d}", ""))

    for i in range(1, args.n_pos_resistance + 1):
        mode = resistance_modes[(i - 1) % len(resistance_modes)][0]
        plan.append(("POS_CMV_RESISTANCE", f"RES{i:03d}", mode))

    for i in range(1, args.n_pos_res_co + 1):
        mode = resistance_modes[(i - 1) % len(resistance_modes)][0]
        plan.append(("POS_RES_PLUS_OTHER", f"RESCO{i:03d}", mode))

    # Prepare output truth file
    meta_path = out_truth / "metadata.tsv"
    with open(meta_path, "w", encoding="utf-8") as meta:
        meta.write("\t".join([
            "sample","group","truth_cmv_present",
            "cmv_ref","cmv_pairs","cmv_fraction_target",
            "other_ref","other_pairs","other_fraction_target",
            "human_ref","human_pairs",
            "expected_resistance",
            "total_pairs","read_len","seq_sys","mean_frag","sd_frag"
        ]) + "\n")

        tmp_root = Path(tempfile.mkdtemp(prefix="cmv_gui_sim_"))
        try:
            for si, (group, sample, res_mode) in enumerate(plan, start=1):
                total_pairs = sample_lognormal_pairs(
                    rng,
                    median=args.median_pairs,
                    spread=args.pair_spread,
                    min_pairs=args.min_pairs,
                    max_pairs=args.max_pairs,
                )
                read_len = rng.choice(read_lens)
                seq_sys = rng.choice(seq_sys_list)
                mean_frag = rng.choice(mean_frags)
                sd_frag = rng.choice(sd_frags)

                r1 = out_fastq / f"{sample}_R1.fastq"
                r2 = out_fastq / f"{sample}_R2.fastq"
                if r1.exists(): r1.unlink()
                if r2.exists(): r2.unlink()

                cmv_present = 0
                expected_res = ""
                cmv_pairs = 0
                other_pairs = 0
                human_pairs = 0

                cmv_ref_path: Path | None = None
                other_ref_path: Path | None = None

                # Choose refs
                if group in {"POS_CMV_WT","POS_CMV_PLUS_OTHER","POS_CMV_RESISTANCE","POS_RES_PLUS_OTHER","NEG_CONTAM"}:
                    cmv_ref_path = rng.choice(cmv_candidates)

                if group in {"NEG_OTHER_VIRUS","POS_CMV_PLUS_OTHER","POS_RES_PLUS_OTHER","NEG_CONTAM"} and other_candidates:
                    other_ref_path = rng.choice(other_candidates)

                # Fractions
                if group == "NEG_HUMAN_ONLY":
                    cmv_present = 0

                elif group == "NEG_OTHER_VIRUS":
                    cmv_present = 0
                    other_frac = rng.uniform(args.other_frac_min, args.other_frac_max) if other_ref_path else 0.0
                    other_pairs = max(1, int(total_pairs * other_frac)) if other_ref_path else 0

                elif group == "NEG_CONTAM":
                    cmv_present = 0
                    contam_frac = 10 ** rng.uniform(math.log10(args.neg_contam_min_frac), math.log10(args.neg_contam_max_frac))
                    cmv_pairs = max(1, int(total_pairs * contam_frac))
                    # optionally add other virus too (if available)
                    other_frac = rng.uniform(0.05, 0.25) if other_ref_path else 0.0
                    other_pairs = max(1, int(total_pairs * other_frac)) if other_ref_path else 0

                elif group == "POS_CMV_WT":
                    cmv_present = 1
                    cmv_frac = 10 ** rng.uniform(math.log10(args.cmv_frac_min), math.log10(args.cmv_frac_max))
                    cmv_pairs = max(1, int(total_pairs * cmv_frac))

                elif group == "POS_CMV_PLUS_OTHER":
                    cmv_present = 1
                    cmv_frac = 10 ** rng.uniform(math.log10(args.cmv_frac_min), math.log10(args.cmv_frac_max))
                    cmv_pairs = max(1, int(total_pairs * cmv_frac))
                    other_frac = rng.uniform(args.other_frac_min, args.other_frac_max) if other_ref_path else 0.0
                    other_pairs = max(1, int(total_pairs * other_frac)) if other_ref_path else 0

                elif group == "POS_CMV_RESISTANCE":
                    cmv_present = 1
                    cmv_frac = 10 ** rng.uniform(math.log10(args.cmv_frac_min), math.log10(args.cmv_frac_max))
                    cmv_pairs = max(1, int(total_pairs * cmv_frac))
                    expected_res = dict((m[0], m[2]) for m in resistance_modes).get(res_mode, "")

                elif group == "POS_RES_PLUS_OTHER":
                    cmv_present = 1
                    cmv_frac = 10 ** rng.uniform(math.log10(args.cmv_frac_min), math.log10(args.cmv_frac_max))
                    cmv_pairs = max(1, int(total_pairs * cmv_frac))
                    other_frac = rng.uniform(args.other_frac_min, args.other_frac_max) if other_ref_path else 0.0
                    other_pairs = max(1, int(total_pairs * other_frac)) if other_ref_path else 0
                    expected_res = dict((m[0], m[2]) for m in resistance_modes).get(res_mode, "")

                else:
                    raise RuntimeError(f"Unknown group: {group}")

                viral_pairs = cmv_pairs + other_pairs
                human_pairs = max(0, total_pairs - viral_pairs) if human_ref is not None else 0

                tmp = tmp_root / sample
                tmp.mkdir(parents=True, exist_ok=True)

                # ---- Build CMV mutant ref if needed ----
                cmv_ref_used_name = ""
                if cmv_ref_path is not None:
                    hdr, merlin_seq = read_fasta_any(cmv_ref_path)
                    cmv_seq = merlin_seq

                    if res_mode:
                        # Apply specific edits on Merlin-like coordinate system
                        if res_mode == "UL97_M460I":
                            codon_start = UL97_START + (460 - 1) * 3
                            cmv_seq, _, _ = replace_codon_to_aa(cmv_seq, codon_start, "I")
                        elif res_mode == "UL97_A594V":
                            codon_start = UL97_START + (594 - 1) * 3
                            cmv_seq, _, _ = replace_codon_to_aa(cmv_seq, codon_start, "V")
                        elif res_mode == "UL97_L595S":
                            codon_start = UL97_START + (595 - 1) * 3
                            cmv_seq, _, _ = replace_codon_to_aa(cmv_seq, codon_start, "S")
                        elif res_mode == "UL54_V812L":
                            codon_start = UL54_START + (812 - 1) * 3
                            cmv_seq, _, _ = replace_codon_to_aa(cmv_seq, codon_start, "L")
                        elif res_mode == "UL54_A834P":
                            codon_start = UL54_START + (834 - 1) * 3
                            cmv_seq, _, _ = replace_codon_to_aa(cmv_seq, codon_start, "P")
                        elif res_mode == "UL54_DEL981_982":
                            del_start = UL54_START + (981 - 1) * 3
                            cmv_seq, _ = delete_aa_block(cmv_seq, del_start, aa_count=2)
                        else:
                            raise RuntimeError(f"Unknown resistance mode: {res_mode}")

                    cmv_out = tmp / ("cmv_mut.fa" if res_mode else "cmv_wt.fa")
                    write_fasta(cmv_out, hdr, cmv_seq)
                    cmv_ref_used_name = cmv_ref_path.name + (f":{res_mode}" if res_mode else ":WT")

                    if cmv_pairs > 0:
                        prefix = tmp / "cmv_"
                        run_art(cmv_out, prefix, cmv_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 10_000 + si)
                        append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # ---- Other virus ----
                other_ref_used_name = ""
                if other_ref_path is not None and other_pairs > 0:
                    hdr_o, seq_o = read_fasta_any(other_ref_path)
                    oth_out = tmp / "other.fa"
                    write_fasta(oth_out, hdr_o, seq_o)
                    other_ref_used_name = other_ref_path.name

                    prefix = tmp / "oth_"
                    run_art(oth_out, prefix, other_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 20_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # ---- Human ----
                if human_ref is not None and human_pairs > 0:
                    hdr_h, seq_h = read_fasta_any(human_ref)
                    hum_out = tmp / "human.fa"
                    write_fasta(hum_out, hdr_h, seq_h)

                    prefix = tmp / "hum_"
                    run_art(hum_out, prefix, human_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 30_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # ---- Pad (adds more human) ----
                if human_ref is not None and mb(r1) < args.min_fastq_mb:
                    while mb(r1) < args.min_fastq_mb:
                        add_pairs = 50_000
                        hdr_h, seq_h = read_fasta_any(human_ref)
                        hum_pad = tmp / "human_pad.fa"
                        write_fasta(hum_pad, hdr_h, seq_h)

                        prefix = tmp / f"pad_{int(mb(r1))}_"
                        run_art(hum_pad, prefix, add_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 40_000 + si)
                        append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)
                        human_pairs += add_pairs
                        total_pairs += add_pairs

                cmv_frac_target = (cmv_pairs / total_pairs) if total_pairs else 0.0
                other_frac_target = (other_pairs / total_pairs) if total_pairs else 0.0

                meta.write("\t".join([
                    sample,
                    group,
                    "1" if cmv_present else "0",
                    cmv_ref_used_name,
                    str(cmv_pairs),
                    f"{cmv_frac_target:.6g}",
                    other_ref_used_name,
                    str(other_pairs),
                    f"{other_frac_target:.6g}",
                    str(human_ref) if human_ref is not None else "",
                    str(human_pairs),
                    expected_res,
                    str(total_pairs),
                    str(read_len),
                    seq_sys,
                    str(mean_frag),
                    str(sd_frag),
                ]) + "\n")

                print(f"✅ {sample:8s} {group:18s} cmv_pairs={cmv_pairs:<8d} other_pairs={other_pairs:<8d} human_pairs={human_pairs:<8d} R1={mb(r1):.1f}MB")

        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    print("\n🎉 Done.")
    print(f"FASTQ folder: {out_fastq}")
    print(f"Truth file:   {meta_path}")

if __name__ == "__main__":
    main()
