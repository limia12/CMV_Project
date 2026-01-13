# #!/usr/bin/env bash
# set -euo pipefail

# # ========== CONFIG ========== #
# REF_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/References/Simulated_Reads_Ref"
# OUT_DIR="${REF_DIR}/out_fastq/CMV_Resistant_Mix"
# mkdir -p "$OUT_DIR/tmp"

# READ_LENGTH=150
# FRAG_MEAN=350
# FRAG_SD=50
# PROFILE="HS25"
# PAIRED="-p"
# SAMFLAG="-na"

# CMV=("AD169.fasta" "Towne.fasta" "hcmv_genome.fasta")  # Merlin = hcmv_genome.fasta

# # ========== FUNCTIONS ========== #

# insert_merlin_resistance() {
#     local inref="$1"
#     local outref="$2"

#     python3 - <<EOF
# from Bio import SeqIO
# from Bio.Seq import Seq

# record = SeqIO.read("$inref", "fasta")
# seq = list(str(record.seq))

# # UL97 edits
# replacements = {
#     143173: "ATA",  # M460I
#     143353: "CAA",  # H520Q
#     143569: "GGT",  # C592G
#     143575: "GTT",  # A594V
#     143578: "TCA",  # L595S
#     143602: "TGG",  # C603W
# }

# # UL54 edits
# replacements.update({
#     138692: "CTT",  # V812L
#     138758: "CCC",  # A834P
# })

# # UL54 deletion: Δ981–982 (6 nt)
# del_start = 139294
# del_end = del_start + 6

# for pos, codon in replacements.items():
#     seq[pos:pos+3] = list(codon)

# del seq[del_start:del_end]

# record.seq = Seq("".join(seq))
# SeqIO.write(record, "$outref", "fasta")
# EOF
# }

# simulate_reads() {
#     local fasta="$1"
#     local outpref="$2"

#     art_illumina -ss "$PROFILE" $SAMFLAG $PAIRED \
#         -i "$fasta" \
#         -l "$READ_LENGTH" \
#         -f 10 \
#         -m "$FRAG_MEAN" \
#         -s "$FRAG_SD" \
#         -o "$OUT_DIR/tmp/${outpref}." \
#         -rs 42 >/dev/null
# }

# combine_fastqs() {
#     local out="$1"
#     local r1="$OUT_DIR/${out}_R1.fq"
#     local r2="$OUT_DIR/${out}_R2.fq"

#     > "$r1"
#     > "$r2"
#     for fq in $OUT_DIR/tmp/*.1.fq; do cat "$fq" >> "$r1"; done
#     for fq in $OUT_DIR/tmp/*.2.fq; do cat "$fq" >> "$r2"; done

#     gzip -f "$r1" "$r2"
# }

# # ========== MAIN ========== #
# echo "[*] Generating CMV strains with hardcoded resistance mutations..."

# for cmv in "${CMV[@]}"; do
#     strain=$(basename "$cmv" .fasta)
#     inref="$REF_DIR/$cmv"
#     mutref="$OUT_DIR/tmp/${strain}_resistance.fa"

#     if [[ "$cmv" == "hcmv_genome.fasta" ]]; then
#         echo "  - Injecting resistance mutations into Merlin..."
#         insert_merlin_resistance "$inref" "$mutref"
#     else
#         echo "  - Copying $strain unchanged (no mapped UL97/UL54 positions)..."
#         cp "$inref" "$mutref"
#     fi

#     echo "  - Simulating reads for $strain..."
#     simulate_reads "$mutref" "$strain"
# done

# echo "[*] Combining reads into one sample..."
# combine_fastqs "CMV_ALL_RESISTANT"

# echo "[\u2713] Done. Output FASTQs:"
# echo "  $OUT_DIR/CMV_ALL_RESISTANT_R1.fq.gz"
# echo "  $OUT_DIR/CMV_ALL_RESISTANT_R2.fq.gz"

#!/usr/bin/env bash
set -euo pipefail

# ========== CONFIG ========== #
REF_DIR="/home/stp/Bioinformatics/Limia/CMV_Project/References/Simulated_Reads_Ref"
OUT_DIR="${REF_DIR}/out_fastq/CMV_UL54_Resistant"
mkdir -p "$OUT_DIR/tmp"

READ_LENGTH=150
FRAG_MEAN=350
FRAG_SD=50
PROFILE="HS25"
PAIRED="-p"
SAMFLAG="-na"

CMV=("hcmv_genome.fasta")  # Only Merlin for UL54

# ========== FUNCTIONS ========== #

insert_ul54_resistance() {
    local inref="$1"
    local outref="$2"

    python3 - <<EOF
from Bio import SeqIO
from Bio.Seq import Seq

record = SeqIO.read("$inref", "fasta")
seq = list(str(record.seq))

# === UL54 resistance mutations ===
# Codon 812 (nt 79487–79489): GTT (V) → CTT (L)
seq[79486:79489] = list("CTT")

# Codon 834 (nt 79419–79421): GCC (A) → CCC (P)
seq[79418:79421] = list("CCC")

# Δ981–982 (nt 78977–78982): delete 6 bp
del seq[78976:78982]

record.seq = Seq("".join(seq))
SeqIO.write(record, "$outref", "fasta")
EOF
}

simulate_reads() {
    local fasta="$1"
    local outpref="$2"

    art_illumina -ss "$PROFILE" $SAMFLAG $PAIRED \
        -i "$fasta" \
        -l "$READ_LENGTH" \
        -f 10 \
        -m "$FRAG_MEAN" \
        -s "$FRAG_SD" \
        -o "$OUT_DIR/tmp/${outpref}." \
        -rs 54 >/dev/null
}

combine_fastqs() {
    local out="$1"
    local r1="$OUT_DIR/${out}_R1.fq"
    local r2="$OUT_DIR/${out}_R2.fq"

    > "$r1"
    > "$r2"
    for fq in $OUT_DIR/tmp/*.1.fq; do cat "$fq" >> "$r1"; done
    for fq in $OUT_DIR/tmp/*.2.fq; do cat "$fq" >> "$r2"; done

    gzip -f "$r1" "$r2"
}

# ========== MAIN ========== #
echo "[*] Generating Merlin CMV with UL54 resistance mutations..."

for cmv in "${CMV[@]}"; do
    strain=$(basename "$cmv" .fasta)
    inref="$REF_DIR/$cmv"
    mutref="$OUT_DIR/tmp/${strain}_UL54_resistance.fa"

    echo "  - Injecting UL54 mutations into $strain..."
    insert_ul54_resistance "$inref" "$mutref"

    echo "  - Simulating reads for $strain with UL54 mutations..."
    simulate_reads "$mutref" "$strain"
done

echo "[*] Combining reads into one UL54-resistant sample..."
combine_fastqs "CMV_UL54_RESISTANT"

echo "[\u2713] Done. Output FASTQs:"
echo "  $OUT_DIR/CMV_UL54_RESISTANT_R1.fq.gz"
echo "  $OUT_DIR/CMV_UL54_RESISTANT_R2.fq.gz"
