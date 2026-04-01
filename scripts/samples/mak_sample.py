#!/usr/bin/env python3
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
                    # separate contigs with Ns so ART can still run
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
        # if mismatch and not failing, still replace (keeps generator robust across refs)
        s[i] = alt
        applied.append((pos, ref, alt))
    return "".join(s), applied


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


def classify_refs(refs):
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


def pick_human_ref(args, refs, human_guess):
    if args.human_ref:
        human_ref = Path(args.human_ref)
        if not human_ref.exists():
            sys.exit(f"❌ --human-ref not found: {human_ref}")
        return human_ref

    # prefer chr* if present
    chr_candidates = []
    for p in refs:
        if re.search(r"\bchr\d+\b|chr[xy]|chrm", p.name.lower()):
            chr_candidates.append(p)
    if chr_candidates:
        return sorted(chr_candidates)[0]
    if human_guess:
        return human_guess[0]
    return None


# --------------------- main ---------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate a SMALL 'definite' truth sample set (1 sample per scenario) for report-led Results."
    )
    ap.add_argument("--refs-dir", required=True, help="Directory with FASTA/FASTA.GZ refs (CMV strains + other viruses + optional chr*.fa).")
    ap.add_argument("--outdir", required=True, help="Output folder (fastq/ + truth/metadata.tsv).")
    ap.add_argument("--mutations", required=True, help="TSV with POS REF ALT for resistance spike-in sample.")
    ap.add_argument("--human-ref", default="", help="Optional explicit human FASTA(.gz). If not provided, tries chr*.fa in refs-dir.")
    ap.add_argument("--min-cmv-len", type=int, default=150000, help="Minimum length to consider a ref as CMV (default 150k).")
    ap.add_argument("--min-other-len", type=int, default=5000, help="Minimum length to consider a ref as other virus (default 5k).")

    # fixed sequencing settings (keeps samples comparable + smaller)
    ap.add_argument("--total-pairs", type=int, default=500000, help="Total read pairs per sample (default 500k).")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seq-sys", default="HS25")
    ap.add_argument("--mean-frag", type=int, default=350)
    ap.add_argument("--sd-frag", type=int, default=50)

    # composition settings (definite, not grey-zone)
    ap.add_argument("--pos-cmv-frac", type=float, default=0.10, help="CMV fraction for DEF_POS (default 0.10).")
    ap.add_argument("--coinf-cmv-frac", type=float, default=0.08, help="CMV fraction for DEF_COINF (default 0.08).")
    ap.add_argument("--coinf-other-frac", type=float, default=0.12, help="Other virus fraction for DEF_COINF (default 0.12).")
    ap.add_argument("--mixed-cmv-frac", type=float, default=0.12, help="Total CMV fraction for DEF_MIXED (default 0.12).")
    ap.add_argument("--mixed-strains", type=int, default=3, help="Number of CMV strains for DEF_MIXED (2 or 3; default 3).")
    ap.add_argument("--resist-cmv-frac", type=float, default=0.12, help="Total CMV fraction for DEF_RESIST (default 0.12).")
    ap.add_argument("--resist-af", type=float, default=0.30, help="Mutant allele fraction target within CMV reads (default 0.30).")
    ap.add_argument("--mutations-n", type=int, default=6, help="How many mutations to apply from TSV (default 6).")

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

    cmv_guess, human_guess, other_guess = classify_refs(refs)

    cmv_candidates = [p for p in cmv_guess if fasta_len(p) >= args.min_cmv_len]
    other_candidates = [p for p in other_guess if fasta_len(p) >= args.min_other_len]
    human_ref = pick_human_ref(args, refs, human_guess)

    if not cmv_candidates:
        sys.exit("❌ No CMV references detected (need full CMV genomes named merlin/ad169/towne/cmv/hcmv).")
    if not other_candidates:
        print("⚠️ No other-virus refs detected: DEF_COINF will become CMV+human only.", file=sys.stderr)
    if human_ref is None:
        print("⚠️ No human reference found: samples will be viral-only unless you pass --human-ref.", file=sys.stderr)

    muts_all = parse_mutations_tsv(Path(args.mutations))
    muts_n = min(args.mutations_n, len(muts_all))

    # Define the 5 "definite" samples
    plan = [
        ("DEF_NEG",    "DEF_NEG",    0),
        ("DEF_POS",    "DEF_POS",    1),
        ("DEF_COINF",  "DEF_COINF",  1),
        ("DEF_MIXED",  "DEF_MIXED",  1),
        ("DEF_RESIST", "DEF_RESIST", 1),
    ]

    meta_path = out_truth / "metadata.tsv"
    with open(meta_path, "w") as meta:
        meta.write("\t".join([
            "sample", "group", "truth_cmv_present",
            "total_pairs",
            "cmv_refs", "cmv_pairs",
            "other_refs", "other_pairs",
            "human_ref", "human_pairs",
            "read_len", "seq_sys", "mean_frag", "sd_frag",
            "mutation_spiked", "mutation_af_target", "mutations_applied"
        ]) + "\n")

        tmp_root = Path(tempfile.mkdtemp(prefix="def_truth_"))
        try:
            for si, (group, sample, cmv_present) in enumerate(plan, start=1):
                total_pairs = int(args.total_pairs)

                r1 = out_fastq / f"{sample}_R1.fastq"
                r2 = out_fastq / f"{sample}_R2.fastq"
                if r1.exists(): r1.unlink()
                if r2.exists(): r2.unlink()

                cmv_refs_used, cmv_pairs_list = [], []
                other_refs_used, other_pairs_list = [], []
                human_pairs = 0
                mutation_spiked = 0
                mutation_af_target = ""
                mutations_applied = []

                # Decide composition
                cmv_pairs_total = 0
                other_pairs = 0

                if group == "DEF_NEG":
                    cmv_pairs_total = 0
                    other_pairs = 0

                elif group == "DEF_POS":
                    cmv_pairs_total = max(1, int(total_pairs * args.pos_cmv_frac))

                elif group == "DEF_COINF":
                    cmv_pairs_total = max(1, int(total_pairs * args.coinf_cmv_frac))
                    if other_candidates:
                        other_pairs = max(1, int(total_pairs * args.coinf_other_frac))

                elif group == "DEF_MIXED":
                    cmv_pairs_total = max(1, int(total_pairs * args.mixed_cmv_frac))

                elif group == "DEF_RESIST":
                    cmv_pairs_total = max(1, int(total_pairs * args.resist_cmv_frac))
                    mutation_spiked = 1
                    mutation_af_target = str(args.resist_af)

                viral_pairs = cmv_pairs_total + other_pairs
                human_pairs = max(0, total_pairs - viral_pairs) if human_ref is not None else 0

                tmp = tmp_root / sample
                tmp.mkdir(parents=True, exist_ok=True)

                def prep_ref(original: Path, do_resistance_muts: bool):
                    hdr, seq = read_fasta_any(original)
                    applied = []
                    if do_resistance_muts:
                        subset = rng.sample(muts_all, k=muts_n)
                        seq, applied = apply_point_mutations(seq, subset, fail_if_mismatch=False)
                    out_ref = tmp / (original.stem + (".resmut" if do_resistance_muts else ".wt") + ".fa")
                    write_fasta(out_ref, hdr, seq)
                    return out_ref, applied

                # --- CMV components ---
                if group in ("DEF_POS", "DEF_COINF"):
                    cmv_ref = rng.choice(cmv_candidates)
                    cmv_refs_used.append(cmv_ref.name)
                    cmv_pairs_list.append(str(cmv_pairs_total))

                    ref_prepped, applied = prep_ref(cmv_ref, do_resistance_muts=False)
                    prefix = tmp / "cmv_1_"
                    run_art(ref_prepped, prefix, cmv_pairs_total, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 10_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                elif group == "DEF_MIXED":
                    n = 2 if args.mixed_strains <= 2 else 3
                    n = min(n, len(cmv_candidates))
                    picks = rng.sample(cmv_candidates, k=n)

                    # split CMV pairs evenly (definite mixed)
                    per = max(1, cmv_pairs_total // n)
                    remainder = cmv_pairs_total - per * n

                    for idx, ref in enumerate(picks, start=1):
                        pairs = per + (1 if idx <= remainder else 0)
                        cmv_refs_used.append(ref.name)
                        cmv_pairs_list.append(str(pairs))

                        ref_prepped, applied = prep_ref(ref, do_resistance_muts=False)
                        prefix = tmp / f"cmv_{idx}_"
                        run_art(ref_prepped, prefix, pairs, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 10_000 + si + idx)
                        append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                elif group == "DEF_RESIST":
                    cmv_ref = rng.choice(cmv_candidates)

                    # Mix WT and mutant reads at specified AF *within CMV reads*
                    mut_pairs = max(1, int(cmv_pairs_total * args.resist_af))
                    wt_pairs = max(1, cmv_pairs_total - mut_pairs)

                    # WT
                    cmv_refs_used.append(cmv_ref.name + ":WT")
                    cmv_pairs_list.append(str(wt_pairs))
                    ref_wt, _ = prep_ref(cmv_ref, do_resistance_muts=False)
                    prefix_wt = tmp / "cmv_wt_"
                    run_art(ref_wt, prefix_wt, wt_pairs, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 11_000 + si)
                    append_fastq(Path(str(prefix_wt) + "1.fq"), Path(str(prefix_wt) + "2.fq"), r1, r2)

                    # MUT
                    cmv_refs_used.append(cmv_ref.name + f":MUT_AF{args.resist_af}")
                    cmv_pairs_list.append(str(mut_pairs))
                    ref_mut, applied = prep_ref(cmv_ref, do_resistance_muts=True)
                    if applied:
                        mutations_applied.extend([f"{pos}{refb}>{altb}" for pos, refb, altb in applied])

                    prefix_mut = tmp / "cmv_mut_"
                    run_art(ref_mut, prefix_mut, mut_pairs, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 12_000 + si)
                    append_fastq(Path(str(prefix_mut) + "1.fq"), Path(str(prefix_mut) + "2.fq"), r1, r2)

                # --- Other virus component ---
                if group == "DEF_COINF" and other_pairs > 0 and other_candidates:
                    oth_ref = rng.choice(other_candidates)
                    other_refs_used.append(oth_ref.name)
                    other_pairs_list.append(str(other_pairs))

                    ref_oth, _ = prep_ref(oth_ref, do_resistance_muts=False)
                    prefix = tmp / "oth_1_"
                    run_art(ref_oth, prefix, other_pairs, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 20_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                # --- Human component ---
                if human_pairs > 0 and human_ref is not None:
                    hdr, seq = read_fasta_any(human_ref)
                    ref_h = tmp / "human.fa"
                    write_fasta(ref_h, hdr, seq)
                    prefix = tmp / "hum_"
                    run_art(ref_h, prefix, human_pairs, args.read_len, args.mean_frag, args.sd_frag, args.seq_sys, args.seed + 30_000 + si)
                    append_fastq(Path(str(prefix) + "1.fq"), Path(str(prefix) + "2.fq"), r1, r2)

                meta.write("\t".join([
                    sample,
                    group,
                    str(int(cmv_present)),
                    str(total_pairs),
                    ",".join(cmv_refs_used) if cmv_refs_used else "",
                    ",".join(cmv_pairs_list) if cmv_pairs_list else "0",
                    ",".join(other_refs_used) if other_refs_used else "",
                    ",".join(other_pairs_list) if other_pairs_list else "0",
                    str(human_ref) if human_ref is not None else "",
                    str(human_pairs),
                    str(args.read_len),
                    args.seq_sys,
                    str(args.mean_frag),
                    str(args.sd_frag),
                    str(int(mutation_spiked)),
                    mutation_af_target,
                    ",".join(mutations_applied) if mutations_applied else ""
                ]) + "\n")

                print(f"✅ {sample:10s} {group:10s} truth={cmv_present} pairs={total_pairs}")

        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    print("\n🎉 Done.")
    print(f"FASTQ: {out_fastq}")
    print(f"TRUTH: {meta_path}")


if __name__ == "__main__":
    main()