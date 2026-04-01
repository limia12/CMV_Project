#!/usr/bin/env python3
import argparse
import gzip
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------- helpers ----------
def which_or_die(cmd: str):
    if shutil.which(cmd) is None:
        sys.exit(f"❌ Required tool not found in PATH: {cmd}\n   Install it or load your module, then retry.")

def read_fasta_any(path: Path) -> tuple[str, str]:
    """Return (header, sequence) for single-contig fasta. If multi-contig, concatenates with Ns."""
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
                    # new contig
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
    """
    Accepts TSV with columns like:
      POS REF ALT
    or
      pos ref alt
    or
      CHROM POS REF ALT
    We only use POS/REF/ALT.
    """
    muts = []
    with open(tsv_path, "r") as fh:
        header = fh.readline().strip().split("\t")
        cols = {c.lower(): i for i, c in enumerate(header)}
        need = {"pos", "ref", "alt"}
        if not need.issubset(cols.keys()):
            # allow 2nd line no header format: POS REF ALT
            # rewind and try parse raw
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
                raise ValueError(f"Could not parse mutations TSV columns in {tsv_path}")
            return muts

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
        if s[i] != ref:
            if fail_if_mismatch:
                raise ValueError(f"REF mismatch at {pos}: expected {ref}, found {s[i]}")
            # still apply (but record mismatch)
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
            choices = [x for x in bases if x != b]
            s[i] = rng.choice(choices)
            n += 1
    return "".join(s), n

def run_art(ref_fa: Path, out_prefix: Path, pairs: int, read_len: int, mean_frag: int, sd_frag: int,
            seq_sys: str, seed: int):
    """
    ART outputs:
      <prefix>1.fq and <prefix>2.fq (paired)
    We use -c for number of read pairs.
    """
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

# ---------- classification ----------
def classify_refs(refs):
    """
    Heuristic:
      - CMV: filename contains merlin/ad169/towne/cmv/hcmv
      - HUMAN: filename starts with chr or contains hg/human
      - OTHER: EBV/HHV and anything not CMV/HUMAN but long enough
    """
    cmv = []
    human = []
    other = []
    for p in refs:
        name = p.name.lower()
        if re.search(r"\b(chr\d+|chr[xy]|chrm)\b", name) or "human" in name or "hg" in name:
            human.append(p)
        elif ("merlin" in name) or ("ad169" in name) or ("towne" in name) or ("cmv" in name) or ("hcmv" in name):
            cmv.append(p)
        else:
            other.append(p)
    return cmv, human, other

# ---------- plan ----------
def make_plan():
    """
    30 samples:
      10 CMV-neg: 5 human-only, 5 human+other-virus
      5 CMV+other-virus
      5 mixed CMV (2-3 strains)
      5 2-strain CMV with varying proportions
      5 mutation-spiked CMV
    """
    plan = []
    # 10 negatives
    for i in range(1, 6):
        plan.append(("NEG_HUMAN_ONLY", f"NEG_H{i:02d}"))
    for i in range(1, 6):
        plan.append(("NEG_OTHER_VIRUS", f"NEG_O{i:02d}"))
    # 5 CMV + other viruses
    for i in range(1, 6):
        plan.append(("CMV_PLUS_OTHER", f"COINF{i:02d}"))
    # 5 mixed CMV (2-3)
    for i in range(1, 6):
        plan.append(("CMV_MIXED_2_3", f"MIX{i:02d}"))
    # 5 varying CMV strain proportions (2 strains)
    for i in range(1, 6):
        plan.append(("CMV_PROP_SERIES", f"PROP{i:02d}"))
    # 5 mutation spiked
    for i in range(1, 6):
        plan.append(("CMV_MUT_SPIKE", f"MUT{i:02d}"))
    return plan

def main():
    ap = argparse.ArgumentParser(
        description="Generate a diverse 30-sample ART Illumina stress-test set with truth metadata."
    )
    ap.add_argument("--refs-dir", required=True, help="Directory containing FASTA/FASTA.GZ references (CMV + other viruses + possibly chr*.fasta).")
    ap.add_argument("--outdir", required=True, help="Output folder.")
    ap.add_argument("--mutations", required=True, help="TSV of resistance mutations (POS REF ALT...). Used for mutation-spike group.")
    ap.add_argument("--human-ref", default="", help="Optional: explicit human FASTA (or .gz) to use instead of chr*.fasta in refs-dir.")
    ap.add_argument("--min-ref-len", type=int, default=150000, help="Minimum reference length to include as viral genome (default 150k).")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seq-sys", default="HS25")
    ap.add_argument("--seed", type=int, default=1)

    # size control
    ap.add_argument("--total-pairs", type=int, default=600000, help="Total read pairs per sample (default 600k => usually >100MB/sample uncompressed).")
    ap.add_argument("--min-fastq-mb", type=float, default=100.0, help="If final R1 < this, pad with extra human reads.")

    # diversity knobs
    ap.add_argument("--snv-rate-min", type=float, default=1e-5)
    ap.add_argument("--snv-rate-max", type=float, default=2e-3)
    ap.add_argument("--mean-frag-list", default="250,350,500")
    ap.add_argument("--sd-frag-list", default="30,50,80")

    args = ap.parse_args()

    which_or_die("art_illumina")

    rng = random.Random(args.seed)

    refs_dir = Path(args.refs_dir)
    outdir = Path(args.outdir)
    out_fastq = outdir / "fastq"
    outdir.mkdir(parents=True, exist_ok=True)
    out_fastq.mkdir(parents=True, exist_ok=True)

    # discover refs (fa/fasta/fna + gz)
    refs = []
    for ext in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
        refs.extend(refs_dir.glob(ext))
    refs = sorted(set(refs))

    if not refs:
        sys.exit(f"❌ No FASTA found under: {refs_dir}")

    # filter to long genomes for virus sets
    long_refs = []
    for p in refs:
        try:
            L = fasta_len(p)
            if L >= args.min_ref_len:
                long_refs.append(p)
        except Exception:
            continue

    cmv_refs, human_refs_guess, other_refs = classify_refs(long_refs)

    # human ref selection: explicit or from ANY chr*.fasta even if short
    human_ref = None
    if args.human_ref:
        human_ref = Path(args.human_ref)
        if not human_ref.exists():
            sys.exit(f"❌ --human-ref not found: {human_ref}")
    else:
        # allow chr*.fasta even if short: search original refs list
        chr_candidates = []
        for p in refs:
            if re.search(r"\bchr\d+\b|chr[xy]|chrm", p.name.lower()):
                chr_candidates.append(p)
        if chr_candidates:
            human_ref = sorted(chr_candidates)[0]

    if not cmv_refs:
        sys.exit("❌ No CMV references detected. Put merlin/ad169/towne/cmv fasta files in refs-dir (full genomes).")

    if not other_refs:
        print("⚠️ No 'other virus' references detected (EBV/HHV etc). 'other virus' groups will be human-only.", file=sys.stderr)

    if human_ref is None:
        print("⚠️ No human reference found (chr*.fasta) and --human-ref not provided. Samples will NOT include human reads.", file=sys.stderr)

    muts_all = parse_mutations_tsv(Path(args.mutations))
    # we'll spike a small random subset per sample to avoid overdoing it
    # (you can increase if you want)
    muts_subset_size = min(6, len(muts_all))

    mean_frags = [int(x) for x in args.mean_frag_list.split(",") if x.strip()]
    sd_frags   = [int(x) for x in args.sd_frag_list.split(",") if x.strip()]

    # proportion series for 2-strain mixes
    prop_series = [(0.90, 0.10), (0.70, 0.30), (0.50, 0.50), (0.30, 0.70), (0.10, 0.90)]
    # mutation AF targets
    mut_af_targets = [0.02, 0.05, 0.10, 0.30, 0.60]

    plan = make_plan()

    meta_path = outdir / "metadata.tsv"
    with open(meta_path, "w") as meta:
        meta.write("\t".join([
            "sample","group","read_len","mean_frag","sd_frag","total_pairs",
            "cmv_present","cmv_refs","cmv_pairs",
            "other_refs","other_pairs",
            "human_ref","human_pairs",
            "snv_rate","random_snv_count",
            "mutation_spiked","mutation_af_target","mutations_applied"
        ]) + "\n")

        tmp_root = Path(tempfile.mkdtemp(prefix="art_sim_"))
        try:
            for idx, (group, sample) in enumerate(plan, start=1):
                read_len = args.read_len
                mean_frag = rng.choice(mean_frags)
                sd_frag = rng.choice(sd_frags)
                snv_rate = rng.uniform(args.snv_rate_min, args.snv_rate_max)

                # output names
                r1 = out_fastq / f"{sample}_R1.fastq"
                r2 = out_fastq / f"{sample}_R2.fastq"
                # reset
                if r1.exists(): r1.unlink()
                if r2.exists(): r2.unlink()

                # decide composition (pairs per component)
                total_pairs = args.total_pairs

                cmv_present = False
                cmv_components = []  # list of (ref_path, pairs, mutated?, label)
                other_components = []  # list of (ref_path, pairs)
                human_pairs = 0
                other_pairs = 0

                # baseline: always reserve a lot for human unless no human_ref
                # We'll compute viral pairs first then give the rest to human.
                if group == "NEG_HUMAN_ONLY":
                    cmv_present = False
                    other_pairs = 0
                    viral_pairs = 0

                elif group == "NEG_OTHER_VIRUS":
                    cmv_present = False
                    # 10–30% other virus
                    other_pairs = int(total_pairs * rng.uniform(0.10, 0.30)) if other_refs else 0
                    viral_pairs = other_pairs
                    if other_pairs > 0:
                        v = rng.choice(other_refs)
                        other_components.append((v, other_pairs))

                elif group == "CMV_PLUS_OTHER":
                    cmv_present = True
                    # CMV 10–40%, other virus 5–20%
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.10, 0.40))
                    other_pairs = int(total_pairs * rng.uniform(0.05, 0.20)) if other_refs else 0
                    viral_pairs = cmv_pairs_total + other_pairs
                    # pick 1 CMV strain (single) + 1 other virus
                    cmv_ref = rng.choice(cmv_refs)
                    cmv_components.append((cmv_ref, cmv_pairs_total, False, cmv_ref.name))
                    if other_pairs > 0:
                        other_ref = rng.choice(other_refs)
                        other_components.append((other_ref, other_pairs))

                elif group == "CMV_MIXED_2_3":
                    cmv_present = True
                    # CMV 20–60% split across 2–3 strains
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.20, 0.60))
                    n_strains = 2 if len(cmv_refs) == 2 else rng.choice([2,3])
                    picks = rng.sample(cmv_refs, k=min(n_strains, len(cmv_refs)))
                    # Dirichlet-like proportions
                    weights = [rng.random() ** 0.5 for _ in picks]
                    s = sum(weights)
                    props = [w/s for w in weights]
                    for ref, p in zip(picks, props):
                        cmv_components.append((ref, max(1, int(cmv_pairs_total * p)), False, ref.name))
                    viral_pairs = cmv_pairs_total

                elif group == "CMV_PROP_SERIES":
                    cmv_present = True
                    # fixed CMV fraction 40%, but vary strain proportions 90/10 ... 10/90
                    cmv_pairs_total = int(total_pairs * 0.40)
                    if len(cmv_refs) < 2:
                        # fallback: single strain with varying viral load
                        cmv_pairs_total = int(total_pairs * rng.uniform(0.05, 0.60))
                        cmv_ref = rng.choice(cmv_refs)
                        cmv_components.append((cmv_ref, cmv_pairs_total, False, cmv_ref.name))
                    else:
                        a, b = rng.sample(cmv_refs, 2)
                        pA, pB = prop_series[(idx-1) % 5]
                        cmv_components.append((a, int(cmv_pairs_total * pA), False, a.name))
                        cmv_components.append((b, int(cmv_pairs_total * pB), False, b.name))
                    viral_pairs = cmv_pairs_total

                elif group == "CMV_MUT_SPIKE":
                    cmv_present = True
                    # CMV 30–60%, spike mutations by mixing WT and MUT reads to hit target AF
                    cmv_pairs_total = int(total_pairs * rng.uniform(0.30, 0.60))
                    af_target = mut_af_targets[(idx-1) % 5]
                    # choose 1 CMV ref as the base
                    cmv_ref = rng.choice(cmv_refs)
                    mut_pairs = max(1, int(cmv_pairs_total * af_target))
                    wt_pairs  = max(1, cmv_pairs_total - mut_pairs)
                    cmv_components.append((cmv_ref, wt_pairs, False, cmv_ref.name + ":WT"))
                    cmv_components.append((cmv_ref, mut_pairs, True,  cmv_ref.name + f":MUT_AF{af_target}"))
                    viral_pairs = cmv_pairs_total
                else:
                    # should not happen
                    viral_pairs = 0

                # human gets the rest
                if human_ref is not None:
                    human_pairs = max(0, total_pairs - viral_pairs)
                else:
                    # no human ref available
                    human_pairs = 0

                # simulate components into temp and append
                tmp = tmp_root / sample
                tmp.mkdir(parents=True, exist_ok=True)

                random_snv_count = 0
                mutations_applied = []
                mutation_spiked = (group == "CMV_MUT_SPIKE")
                mutation_af_target = ""
                if mutation_spiked:
                    mutation_af_target = str(mut_af_targets[(idx-1) % 5])

                def prep_ref(original: Path, do_random_snvs: bool, do_resistance_muts: bool):
                    hdr, seq = read_fasta_any(original)
                    if do_random_snvs:
                        seq, n = apply_random_snvs(seq, snv_rate, rng)
                        nonlocal_random = n
                    else:
                        nonlocal_random = 0

                    applied = []
                    if do_resistance_muts:
                        # choose a subset so it’s not insane
                        subset = rng.sample(muts_all, k=muts_subset_size)
                        seq2, applied = apply_point_mutations(seq, subset, fail_if_mismatch=False)
                        seq = seq2
                    out_ref = tmp / (original.stem + (".mut" if do_resistance_muts else ".snv") + ".fa")
                    write_fasta(out_ref, hdr, seq)
                    return out_ref, nonlocal_random, applied

                # CMV + other virus refs: we always apply random SNVs to viral refs for diversity
                cmv_refs_used = []
                cmv_pairs_list = []

                for comp_i, (ref, pairs, is_mut, label) in enumerate(cmv_components, start=1):
                    cmv_refs_used.append(label)
                    cmv_pairs_list.append(str(pairs))

                    ref_prepped, n_snv, applied = prep_ref(ref, do_random_snvs=True, do_resistance_muts=is_mut)
                    random_snv_count += n_snv
                    if applied:
                        mutations_applied.extend([f"{pos}{refb}>{altb}" for pos, refb, altb in applied])

                    prefix = tmp / f"cmv_{comp_i}_"
                    run_art(ref_prepped, prefix, pairs, read_len, mean_frag, sd_frag, args.seq_sys, args.seed + idx + comp_i)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                other_refs_used = []
                other_pairs_list = []
                for comp_i, (ref, pairs) in enumerate(other_components, start=1):
                    other_refs_used.append(ref.name)
                    other_pairs_list.append(str(pairs))
                    ref_prepped, n_snv, _ = prep_ref(ref, do_random_snvs=True, do_resistance_muts=False)
                    random_snv_count += n_snv
                    prefix = tmp / f"oth_{comp_i}_"
                    run_art(ref_prepped, prefix, pairs, read_len, mean_frag, sd_frag, args.seq_sys, args.seed + 1000 + idx + comp_i)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                if human_pairs > 0 and human_ref is not None:
                    # For human, we do NOT random-SNV the reference (keeps mapping stable)
                    hdr, seq = read_fasta_any(human_ref)
                    ref_h = tmp / "human.fa"
                    write_fasta(ref_h, hdr, seq)
                    prefix = tmp / "hum_"
                    run_art(ref_h, prefix, human_pairs, read_len, mean_frag, sd_frag, args.seq_sys, args.seed + 2000 + idx)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # pad if too small
                if human_ref is not None and mb(r1) < args.min_fastq_mb:
                    # add extra human reads until we cross threshold (in chunks)
                    while mb(r1) < args.min_fastq_mb:
                        add_pairs = 50000
                        hdr, seq = read_fasta_any(human_ref)
                        ref_h = tmp / "human_pad.fa"
                        write_fasta(ref_h, hdr, seq)
                        prefix = tmp / f"pad_{int(mb(r1))}_"
                        run_art(ref_h, prefix, add_pairs, read_len, mean_frag, sd_frag, args.seq_sys, args.seed + 3000 + idx)
                        append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)
                        human_pairs += add_pairs
                        total_pairs += add_pairs

                # write metadata row
                meta.write("\t".join([
                    sample,
                    group,
                    str(read_len),
                    str(mean_frag),
                    str(sd_frag),
                    str(total_pairs),
                    "1" if cmv_present else "0",
                    ",".join(cmv_refs_used) if cmv_refs_used else "",
                    ",".join(cmv_pairs_list) if cmv_pairs_list else "",
                    ",".join(other_refs_used) if other_refs_used else "",
                    ",".join(other_pairs_list) if other_pairs_list else "",
                    str(human_ref) if human_ref is not None else "",
                    str(human_pairs),
                    f"{snv_rate:.6g}",
                    str(random_snv_count),
                    "1" if mutation_spiked else "0",
                    mutation_af_target,
                    ",".join(mutations_applied) if mutations_applied else ""
                ]) + "\n")

                print(f"✅ {sample:8s} {group:16s} R1={mb(r1):.1f}MB  R2={mb(r2):.1f}MB")

        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    print(f"\n🎉 Done.\nFASTQ:   {out_fastq}\nTRUTH:   {meta_path}")

if __name__ == "__main__":
    main()
