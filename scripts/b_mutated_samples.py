#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Merlin ORF starts (same as your app’s CONSERVE_REGIONS starts)
UL97_START = 141797
UL54_START = 78193

# -------------------- FASTA helpers --------------------
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
                    # separate contigs (rare for Merlin; fine)
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

# -------------------- FASTQ concatenation --------------------
def append_fastq(src1: Path, src2: Path, dst1: Path, dst2: Path):
    with open(dst1, "ab") as o1, open(dst2, "ab") as o2:
        with open(src1, "rb") as i1, open(src2, "rb") as i2:
            shutil.copyfileobj(i1, o1)
            shutil.copyfileobj(i2, o2)

# -------------------- ART wrapper --------------------
def run_art(ref_fa: Path, out_prefix: Path, pairs: int, read_len: int, mean_frag: int, sd_frag: int, seq_sys: str, seed: int):
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

def which_or_die(cmd: str):
    if shutil.which(cmd) is None:
        sys.exit(f"❌ Required tool not found in PATH: {cmd}")

# -------------------- Mutation helpers --------------------
CODON_TABLE = {
    # minimal needed AA sets
    "I": {"ATA","ATT","ATC"},
    "M": {"ATG"},
    "V": {"GTT","GTC","GTA","GTG"},
    "L": {"TTA","TTG","CTT","CTC","CTA","CTG"},
}

def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)

def codon_at(seq: str, pos1: int) -> str:
    # pos1 is 1-based, points to codon first base
    i = pos1 - 1
    return seq[i:i+3]

def replace_codon_to_aa(seq: str, codon_start_1based: int, target_aa: str) -> tuple[str, str, str]:
    """
    Replace the existing codon with a codon encoding target_aa, choosing minimal edits.
    Returns: (new_seq, old_codon, new_codon)
    """
    old = codon_at(seq, codon_start_1based)
    if len(old) != 3 or any(b not in "ACGT" for b in old):
        # just force a default codon for target AA
        new = sorted(CODON_TABLE[target_aa])[0]
    else:
        candidates = sorted(CODON_TABLE[target_aa])
        new = min(candidates, key=lambda c: hamming(old, c))
    i = codon_start_1based - 1
    new_seq = seq[:i] + new + seq[i+3:]
    return new_seq, old, new

def delete_aa_block(seq: str, start_1based: int, aa_count: int) -> tuple[str, str]:
    """
    Delete aa_count amino acids starting at start_1based nucleotide (codon boundary assumed).
    Deletes aa_count*3 bases. Returns (new_seq, deleted_bases).
    """
    i = start_1based - 1
    n = aa_count * 3
    deleted = seq[i:i+n]
    new_seq = seq[:i] + seq[i+n:]
    return new_seq, deleted

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser(
        description="Create 1-of-each CMV resistance mixed FASTQ set (Merlin + human chr21) for GUI pipeline testing."
    )
    ap.add_argument("--merlin", required=True, help="Merlin FASTA (NC_006273.2)")
    ap.add_argument("--chr21", required=True, help="Human chr21 FASTA")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--total-pairs", type=int, default=400000, help="Total read pairs per sample (default 400k)")
    ap.add_argument("--cmv-frac", type=float, default=0.05, help="Fraction of pairs from CMV (default 0.05)")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seq-sys", default="HS25")
    ap.add_argument("--mean-frag", type=int, default=350)
    ap.add_argument("--sd-frag", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    which_or_die("art_illumina")

    outdir = Path(args.outdir).resolve()
    fastq_dir = outdir / "fastq"
    truth_dir = outdir / "truth"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    truth_dir.mkdir(parents=True, exist_ok=True)

    merlin_path = Path(args.merlin).resolve()
    chr21_path = Path(args.chr21).resolve()

    merlin_hdr, merlin_seq = read_fasta_any(merlin_path)
    chr21_hdr, chr21_seq = read_fasta_any(chr21_path)

    # sanity: Merlin should be ~235646
    if len(merlin_seq) < 200000:
        print(f"⚠️ Merlin sequence length looks short: {len(merlin_seq)}", file=sys.stderr)

    rng = random.Random(args.seed)

    # Plan: one of each
    samples = [
        ("SUSC001", "none"),
        ("RES_UL97_001", "ul97_m460i"),
        ("RES_UL54_001", "ul54_v812l"),
        ("RES_UL54_DEL_001", "ul54_del981_982"),
    ]

    # Make mutated Merlin references in a temp folder
    tmp_root = Path(tempfile.mkdtemp(prefix="cmv_oneofeach_"))
    try:
        # write chr21 temp fasta (ART wants a file)
        chr21_fa = tmp_root / "chr21.fa"
        write_fasta(chr21_fa, chr21_hdr, chr21_seq)

        # write base merlin
        merlin_base_fa = tmp_root / "merlin.base.fa"
        write_fasta(merlin_base_fa, merlin_hdr, merlin_seq)

        # build mutant refs
        mutant_refs: dict[str, Path] = {"none": merlin_base_fa}

        # UL97 M460I
        seq_ul97 = merlin_seq
        ul97_codon_start = UL97_START + (460 - 1) * 3
        seq_ul97, old_c, new_c = replace_codon_to_aa(seq_ul97, ul97_codon_start, "I")
        ul97_fa = tmp_root / "merlin.UL97_M460I.fa"
        write_fasta(ul97_fa, merlin_hdr, seq_ul97)
        mutant_refs["ul97_m460i"] = ul97_fa

        # UL54 V812L
        seq_ul54 = merlin_seq
        ul54_codon_start = UL54_START + (812 - 1) * 3
        seq_ul54, old_c2, new_c2 = replace_codon_to_aa(seq_ul54, ul54_codon_start, "L")
        ul54_fa = tmp_root / "merlin.UL54_V812L.fa"
        write_fasta(ul54_fa, merlin_hdr, seq_ul54)
        mutant_refs["ul54_v812l"] = ul54_fa

        # UL54 del981_982 (delete 2 amino acids => 6 bases)
        seq_del = merlin_seq
        del_start = UL54_START + (981 - 1) * 3
        seq_del, deleted = delete_aa_block(seq_del, del_start, aa_count=2)
        del_fa = tmp_root / "merlin.UL54_del981_982.fa"
        write_fasta(del_fa, merlin_hdr, seq_del)
        mutant_refs["ul54_del981_982"] = del_fa

        # pairs split
        total_pairs = int(args.total_pairs)
        cmv_pairs = max(1, int(total_pairs * float(args.cmv_frac)))
        human_pairs = max(1, total_pairs - cmv_pairs)

        truth_path = truth_dir / "expected.tsv"
        with truth_path.open("w", encoding="utf-8") as t:
            t.write("sample\tcmv_pairs\thuman_pairs\texpected_flag\n")

            for i, (sample, mode) in enumerate(samples, start=1):
                r1 = fastq_dir / f"{sample}_R1.fastq"
                r2 = fastq_dir / f"{sample}_R2.fastq"
                if r1.exists(): r1.unlink()
                if r2.exists(): r2.unlink()

                # --- CMV component ---
                cmv_ref = mutant_refs[mode]
                cmv_prefix = tmp_root / f"{sample}.cmv_"
                run_art(
                    cmv_ref, cmv_prefix, cmv_pairs,
                    args.read_len, args.mean_frag, args.sd_frag, args.seq_sys,
                    seed=args.seed + 1000 + i
                )
                append_fastq(Path(str(cmv_prefix) + "1.fq"), Path(str(cmv_prefix) + "2.fq"), r1, r2)

                # --- Human chr21 component ---
                hum_prefix = tmp_root / f"{sample}.hum_"
                run_art(
                    chr21_fa, hum_prefix, human_pairs,
                    args.read_len, args.mean_frag, args.sd_frag, args.seq_sys,
                    seed=args.seed + 2000 + i
                )
                append_fastq(Path(str(hum_prefix) + "1.fq"), Path(str(hum_prefix) + "2.fq"), r1, r2)

                if mode == "none":
                    exp = "None detected"
                elif mode == "ul97_m460i":
                    exp = "UL97 M460I"
                elif mode == "ul54_v812l":
                    exp = "UL54 V812L"
                else:
                    exp = "UL54 Δ981–982"

                t.write(f"{sample}\t{cmv_pairs}\t{human_pairs}\t{exp}\n")
                print(f"✅ {sample}: CMV pairs={cmv_pairs} + human(chr21) pairs={human_pairs}  expected={exp}")

        print("\n🎉 Done.")
        print(f"FASTQs: {fastq_dir}")
        print(f"Truth:  {truth_path}")
        print("\nNotes:")
        print(" - UL97 M460I and UL54 V812L are encoded by editing the Merlin codon at the expected ORF position.")
        print(" - UL54 del981_982 deletes 6 bases at the expected codon boundary.")
        print(" - Run these FASTQs through your GUI pipeline; it should produce annotated VCFs that your app flags.")

    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

if __name__ == "__main__":
    main()
