#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

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

def parse_mutations_tsv(tsv_path: Path):
    muts = []
    with open(tsv_path, "r") as fh:
        header = fh.readline().strip().split("\t")
        cols = {c.lower(): i for i, c in enumerate(header)}
        if {"pos", "ref", "alt"}.issubset(cols.keys()):
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                try:
                    pos = int(parts[cols["pos"]])
                    ref = parts[cols["ref"]].upper()
                    alt = parts[cols["alt"]].upper()
                    muts.append((pos, ref, alt))
                except Exception:
                    continue
        else:
            fh.seek(0)
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    pos = int(parts[0])
                    ref = parts[1].upper()
                    alt = parts[2].upper()
                    muts.append((pos, ref, alt))
                except Exception:
                    continue
    if not muts:
        raise ValueError(f"No mutations parsed from {tsv_path}")
    return muts

def apply_point_mutations(seq: str, muts, fail_if_mismatch=False):
    s = list(seq)
    applied = []
    for pos, ref, alt in muts:
        i = pos - 1
        if i < 0 or i >= len(s):
            continue
        if s[i] != ref and fail_if_mismatch:
            raise ValueError(f"REF mismatch at {pos}: expected {ref}, found {s[i]}")
        s[i] = alt
        applied.append((pos, ref, alt))
    return "".join(s), applied

def apply_random_snvs(seq: str, rate: float, rng: random.Random):
    if rate <= 0:
        return seq, 0
    bases = ["A", "C", "G", "T"]
    s = list(seq)
    n = 0
    for i, b in enumerate(s):
        if b not in "ACGT":
            continue
        if rng.random() < rate:
            s[i] = rng.choice([x for x in bases if x != b])
            n += 1
    return "".join(s), n

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

def classify_refs(refs):
    cmv = []
    human = []
    other = []
    for p in refs:
        name = p.name.lower()
        # human heuristics
        if re.search(r"\b(chr\d+|chr[xy]|chrm)\b", name) or "human" in name or "hg" in name:
            human.append(p)
        # cmv heuristics
        elif ("merlin" in name) or ("ad169" in name) or ("towne" in name) or ("cmv" in name) or ("hcmv" in name):
            cmv.append(p)
        else:
            other.append(p)
    return cmv, human, other

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
    # spread ~ 0.6 gives wide depth variation
    x = rng.lognormvariate(mu=math.log(median), sigma=spread)
    x = max(min_pairs, min(max_pairs, int(x)))
    return x

# --------------------- main ---------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate a ROC-ready CMV threshold stress-test sample set (includes hard negatives + contamination)."
    )
    ap.add_argument("--refs-dir", required=True, help="Directory with FASTA/FASTA.GZ refs (CMV strains + other viruses + optional chr*.fa).")
    ap.add_argument("--outdir", required=True, help="Output folder (fastq/ + truth/metadata.tsv).")
    ap.add_argument("--mutations", required=True, help="TSV with POS REF ALT (used for mutation-spike samples).")

    ap.add_argument("--human-ref", default="", help="Optional explicit human FASTA(.gz). If not provided, tries chr*.fa in refs-dir.")

    # IMPORTANT: separate min lengths
    ap.add_argument("--min-cmv-len", type=int, default=150000, help="Minimum length to consider a reference as CMV candidate (default 150k).")
    ap.add_argument("--min-other-len", type=int, default=5000, help="Minimum length to consider as other-virus candidate (default 5k).")

    # size/heterogeneity knobs
    ap.add_argument("--median-pairs", type=int, default=600000, help="Median read pairs per sample (depth varies lognormally).")
    ap.add_argument("--pair-spread", type=float, default=0.6, help="Lognormal sigma controlling depth heterogeneity.")
    ap.add_argument("--min-pairs", type=int, default=120000)
    ap.add_argument("--max-pairs", type=int, default=1400000)

    ap.add_argument("--min-fastq-mb", type=float, default=100.0, help="If R1 < this, pad with extra human reads.")
    ap.add_argument("--read-lens", default="75,100,150", help="Comma list of read lengths to sample from.")
    ap.add_argument("--seq-sys-list", default="HS20,HS25,MSv3", help="Comma list of ART Illumina profiles to sample from.")
    ap.add_argument("--mean-frag-list", default="250,350,500")
    ap.add_argument("--sd-frag-list", default="30,50,80")

    # ROC ladder
    ap.add_argument("--roc-levels", type=int, default=14, help="Number of CMV fraction levels (log-spaced).")
    ap.add_argument("--roc-reps", type=int, default=4, help="Replicates per CMV fraction level.")
    ap.add_argument("--roc-min-frac", type=float, default=1e-6, help="Min CMV fraction for ladder positives.")
    ap.add_argument("--roc-max-frac", type=float, default=0.40, help="Max CMV fraction for ladder positives.")

    # add HARD negatives (truth=0 but tiny CMV spike)
    ap.add_argument("--n-neg-human", type=int, default=12)
    ap.add_argument("--n-neg-othervirus", type=int, default=12)
    ap.add_argument("--n-neg-contam", type=int, default=16, help="Truth-negative samples with tiny CMV spike-in (index hopping / lab contam).")
    ap.add_argument("--neg-contam-min-frac", type=float, default=1e-7)
    ap.add_argument("--neg-contam-max-frac", type=float, default=5e-5)

    # positives
    ap.add_argument("--n-coinf", type=int, default=10)
    ap.add_argument("--n-mixed", type=int, default=10)
    ap.add_argument("--n-mut", type=int, default=10)

    # diversity knobs
    ap.add_argument("--snv-rate-min", type=float, default=1e-5)
    ap.add_argument("--snv-rate-max", type=float, default=3e-3)

    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    which_or_die("art_illumina")
    rng = random.Random(args.seed)

    refs_dir = Path(args.refs_dir)
    outdir = Path(args.outdir)
    out_fastq = outdir / "fastq"
    out_truth = outdir / "truth"
    out_fastq.mkdir(parents=True, exist_ok=True)
    out_truth.mkdir(parents=True, exist_ok=True)

    # discover refs
    refs = []
    for ext in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
        refs.extend(refs_dir.glob(ext))
    refs = sorted(set(refs))
    if not refs:
        sys.exit(f"❌ No FASTA found under: {refs_dir}")

    # classify by name, but enforce min lengths separately
    cmv_candidates = []
    other_candidates = []
    human_candidates = []

    cmv_guess, human_guess, other_guess = classify_refs(refs)

    for p in cmv_guess:
        try:
            if fasta_len(p) >= args.min_cmv_len:
                cmv_candidates.append(p)
        except Exception:
            pass

    for p in other_guess:
        try:
            if fasta_len(p) >= args.min_other_len:
                other_candidates.append(p)
        except Exception:
            pass

    for p in human_guess:
        human_candidates.append(p)

    if not cmv_candidates:
        sys.exit("❌ No CMV references detected. Put merlin/ad169/towne/cmv FASTAs (full genomes) in refs-dir.")

    # pick human ref
    human_ref = None
    if args.human_ref:
        human_ref = Path(args.human_ref)
        if not human_ref.exists():
            sys.exit(f"❌ --human-ref not found: {human_ref}")
    else:
        # prefer a chr* file if present
        chr_candidates = []
        for p in refs:
            if re.search(r"\bchr\d+\b|chr[xy]|chrm", p.name.lower()):
                chr_candidates.append(p)
        if chr_candidates:
            human_ref = sorted(chr_candidates)[0]
        elif human_candidates:
            human_ref = human_candidates[0]

    if human_ref is None:
        print("⚠️ No human reference found. Samples will be viral-only unless you pass --human-ref.", file=sys.stderr)

    if not other_candidates:
        print("⚠️ No other-virus refs detected. 'other-virus' groups will become human-only.", file=sys.stderr)

    muts_all = parse_mutations_tsv(Path(args.mutations))
    muts_subset_size = min(6, len(muts_all))

    read_lens = [int(x) for x in args.read_lens.split(",") if x.strip()]
    seq_sys_list = [x.strip() for x in args.seq_sys_list.split(",") if x.strip()]
    mean_frags = [int(x) for x in args.mean_frag_list.split(",") if x.strip()]
    sd_frags   = [int(x) for x in args.sd_frag_list.split(",") if x.strip()]

    # plan samples
    plan = []

    for i in range(1, args.n_neg_human + 1):
        plan.append(("NEG_HUMAN_ONLY", f"NEG_H{i:03d}", None))

    for i in range(1, args.n_neg_othervirus + 1):
        plan.append(("NEG_OTHER_VIRUS", f"NEG_O{i:03d}", None))

    for i in range(1, args.n_neg_contam + 1):
        plan.append(("NEG_CONTAM", f"NEG_C{i:03d}", None))

    # ROC ladder positives
    levels = logspace_fracs(args.roc_min_frac, args.roc_max_frac, args.roc_levels)
    idx = 1
    for frac in levels:
        for r in range(args.roc_reps):
            plan.append(("POS_ROC_LADDER", f"ROC{idx:03d}", frac))
            idx += 1

    for i in range(1, args.n_coinf + 1):
        plan.append(("POS_CMV_PLUS_OTHER", f"COINF{i:03d}", None))

    for i in range(1, args.n_mixed + 1):
        plan.append(("POS_MIXED_STRAINS", f"MIX{i:03d}", None))

    mut_af_targets = [0.02, 0.05, 0.10, 0.30, 0.60]
    for i in range(1, args.n_mut + 1):
        plan.append(("POS_MUT_SPIKE", f"MUT{i:03d}", mut_af_targets[(i - 1) % len(mut_af_targets)]))

    meta_path = out_truth / "metadata.tsv"
    with open(meta_path, "w") as meta:
        meta.write("\t".join([
            "sample","group","truth_cmv_present",
            "cmv_fraction_target","total_pairs",
            "cmv_refs","cmv_pairs",
            "other_refs","other_pairs",
            "human_ref","human_pairs",
            "read_len","seq_sys","mean_frag","sd_frag",
            "snv_rate","random_snv_count",
            "mutation_spiked","mutation_af_target","mutations_applied"
        ]) + "\n")

        tmp_root = Path(tempfile.mkdtemp(prefix="roc_sim_v2_"))
        try:
            for si, (group, sample, param) in enumerate(plan, start=1):

                # heterogeneity per-sample
                total_pairs = sample_lognormal_pairs(
                    rng, median=args.median_pairs, spread=args.pair_spread,
                    min_pairs=args.min_pairs, max_pairs=args.max_pairs
                )
                read_len = rng.choice(read_lens)
                seq_sys = rng.choice(seq_sys_list)
                mean_frag = rng.choice(mean_frags)
                sd_frag = rng.choice(sd_frags)
                snv_rate = rng.uniform(args.snv_rate_min, args.snv_rate_max)

                r1 = out_fastq / f"{sample}_R1.fastq"
                r2 = out_fastq / f"{sample}_R2.fastq"
                if r1.exists(): r1.unlink()
                if r2.exists(): r2.unlink()

                cmv_components = []   # (ref, pairs, do_res_muts, label)
                other_components = [] # (ref, pairs)
                cmv_present = 0
                mutation_spiked = 0
                mutation_af_target = ""
                mutations_applied = []
                random_snv_count = 0

                other_pairs = 0
                cmv_pairs_total = 0

                if group == "NEG_HUMAN_ONLY":
                    cmv_present = 0

                elif group == "NEG_OTHER_VIRUS":
                    cmv_present = 0
                    other_pairs = int(total_pairs * rng.uniform(0.10, 0.45)) if other_candidates else 0
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))

                elif group == "NEG_CONTAM":
                    # truth-negative, but tiny CMV contamination → creates grey-zone false positives
                    cmv_present = 0
                    contam_frac = 10 ** rng.uniform(math.log10(args.neg_contam_min_frac), math.log10(args.neg_contam_max_frac))
                    cmv_pairs_total = max(1, int(total_pairs * contam_frac))
                    cmv_ref = rng.choice(cmv_candidates)
                    cmv_components.append((cmv_ref, cmv_pairs_total, False, cmv_ref.name + ":CONTAM"))
                    other_pairs = int(total_pairs * rng.uniform(0.05, 0.30)) if other_candidates else 0
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))

                elif group == "POS_ROC_LADDER":
                    cmv_present = 1
                    frac = float(param)
                    cmv_pairs_total = max(1, int(total_pairs * frac))
                    cmv_ref = rng.choice(cmv_candidates)
                    cmv_components.append((cmv_ref, cmv_pairs_total, False, cmv_ref.name))
                    other_pairs = int(total_pairs * rng.uniform(0.00, 0.15)) if other_candidates else 0
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))

                elif group == "POS_CMV_PLUS_OTHER":
                    cmv_present = 1
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.01, 0.25))
                    other_pairs = int(total_pairs * rng.uniform(0.05, 0.30)) if other_candidates else 0
                    cmv_ref = rng.choice(cmv_candidates)
                    cmv_components.append((cmv_ref, cmv_pairs_total, False, cmv_ref.name))
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))

                elif group == "POS_MIXED_STRAINS":
                    cmv_present = 1
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.03, 0.50))
                    n_strains = 2 if len(cmv_candidates) == 2 else rng.choice([2, 3])
                    picks = rng.sample(cmv_candidates, k=min(n_strains, len(cmv_candidates)))
                    weights = [rng.random() ** 0.5 for _ in picks]
                    s = sum(weights)
                    props = [w / s for w in weights]
                    for ref, p in zip(picks, props):
                        cmv_components.append((ref, max(1, int(cmv_pairs_total * p)), False, ref.name))
                    other_pairs = int(total_pairs * rng.uniform(0.00, 0.20)) if other_candidates else 0
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))

                elif group == "POS_MUT_SPIKE":
                    cmv_present = 1
                    mutation_spiked = 1
                    af_target = float(param)
                    mutation_af_target = str(af_target)
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.05, 0.50))
                    mut_pairs = max(1, int(cmv_pairs_total * af_target))
                    wt_pairs  = max(1, cmv_pairs_total - mut_pairs)
                    cmv_ref = rng.choice(cmv_candidates)
                    cmv_components.append((cmv_ref, wt_pairs, False, cmv_ref.name + ":WT"))
                    cmv_components.append((cmv_ref, mut_pairs, True,  cmv_ref.name + f":MUT_AF{af_target}"))
                    other_pairs = int(total_pairs * rng.uniform(0.00, 0.15)) if other_candidates else 0
                    if other_pairs > 0:
                        other_components.append((rng.choice(other_candidates), other_pairs))
                else:
                    raise RuntimeError(f"Unknown group: {group}")

                viral_pairs = cmv_pairs_total + other_pairs
                human_pairs = max(0, total_pairs - viral_pairs) if human_ref is not None else 0

                tmp = tmp_root / sample
                tmp.mkdir(parents=True, exist_ok=True)

                def prep_ref(original: Path, do_random_snvs: bool, do_resistance_muts: bool):
                    hdr, seq = read_fasta_any(original)
                    n_snv = 0
                    if do_random_snvs:
                        seq, n_snv = apply_random_snvs(seq, snv_rate, rng)
                    applied = []
                    if do_resistance_muts:
                        subset = rng.sample(muts_all, k=muts_subset_size)
                        seq, applied = apply_point_mutations(seq, subset, fail_if_mismatch=False)
                    out_ref = tmp / (original.stem + (".mut" if do_resistance_muts else ".snv") + ".fa")
                    write_fasta(out_ref, hdr, seq)
                    return out_ref, n_snv, applied

                cmv_refs_used, cmv_pairs_list = [], []
                for ci, (ref, pairs, do_muts, label) in enumerate(cmv_components, start=1):
                    cmv_refs_used.append(label)
                    cmv_pairs_list.append(str(pairs))
                    ref_prepped, n_snv, applied = prep_ref(ref, do_random_snvs=True, do_resistance_muts=do_muts)
                    random_snv_count += n_snv
                    if applied:
                        mutations_applied.extend([f"{pos}{refb}>{altb}" for pos, refb, altb in applied])

                    prefix = tmp / f"cmv_{ci}_"
                    run_art(ref_prepped, prefix, pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 10_000 + si + ci)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                other_refs_used, other_pairs_list = [], []
                for ci, (ref, pairs) in enumerate(other_components, start=1):
                    other_refs_used.append(ref.name)
                    other_pairs_list.append(str(pairs))
                    ref_prepped, n_snv, _ = prep_ref(ref, do_random_snvs=True, do_resistance_muts=False)
                    random_snv_count += n_snv

                    prefix = tmp / f"oth_{ci}_"
                    run_art(ref_prepped, prefix, pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 20_000 + si + ci)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                if human_pairs > 0 and human_ref is not None:
                    hdr, seq = read_fasta_any(human_ref)
                    ref_h = tmp / "human.fa"
                    write_fasta(ref_h, hdr, seq)
                    prefix = tmp / "hum_"
                    run_art(ref_h, prefix, human_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 30_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # pad to minimum size (adds more human only)
                if human_ref is not None and mb(r1) < args.min_fastq_mb:
                    while mb(r1) < args.min_fastq_mb:
                        add_pairs = 50_000
                        hdr, seq = read_fasta_any(human_ref)
                        ref_h = tmp / "human_pad.fa"
                        write_fasta(ref_h, hdr, seq)
                        prefix = tmp / f"pad_{int(mb(r1))}_"
                        run_art(ref_h, prefix, add_pairs, read_len, mean_frag, sd_frag, seq_sys, args.seed + 40_000 + si)
                        append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)
                        human_pairs += add_pairs
                        total_pairs += add_pairs

                cmv_fraction_target = (cmv_pairs_total / total_pairs) if total_pairs else 0.0

                meta.write("\t".join([
                    sample,
                    group,
                    "1" if cmv_present else "0",
                    f"{cmv_fraction_target:.6g}" if cmv_pairs_total > 0 else "0",
                    str(total_pairs),
                    ",".join(cmv_refs_used),
                    ",".join(cmv_pairs_list) if cmv_pairs_list else "0",
                    ",".join(other_refs_used),
                    ",".join(other_pairs_list) if other_pairs_list else "0",
                    str(human_ref) if human_ref is not None else "",
                    str(human_pairs),
                    str(read_len),
                    seq_sys,
                    str(mean_frag),
                    str(sd_frag),
                    f"{snv_rate:.6g}",
                    str(random_snv_count),
                    "1" if mutation_spiked else "0",
                    mutation_af_target,
                    ",".join(mutations_applied) if mutations_applied else ""
                ]) + "\n")

                print(f"✅ {sample:8s} {group:16s} truth={cmv_present} cmv_frac={cmv_fraction_target:.3g} R1={mb(r1):.1f}MB")

        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    print("\n🎉 Done.")
    print(f"FASTQ: {out_fastq}")
    print(f"TRUTH: {meta_path}")

if __name__ == "__main__":
    main()
