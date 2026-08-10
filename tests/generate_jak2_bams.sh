#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
output_dir="$script_dir/data/jak2"
mkdir -p "$output_dir"

generate_bam() {
    sample=$1
    read_count=$2
    alt_count=$3
    seed=$4
    contig=$5
    output="$output_dir/$sample.bam"

    awk -v sample="$sample" -v read_count="$read_count" -v alt_count="$alt_count" -v seed="$seed" -v contig="$contig" '
        function reference_base(position, base_index) {
            if (position == 5073770) return "G"
            base_index = ((position * 7 + 3) % 4) + 1
            return substr("ACGT", base_index, 1)
        }
        BEGIN {
            srand(seed)
            print "@HD\tVN:1.6\tSO:coordinate"
            print "@SQ\tSN:" contig "\tLN:141213431"
            print "@RG\tID:" sample "\tSM:" sample "\tPL:ILLUMINA"
            target = 5073770
            read_length = 151
            quality = ""
            for (quality_index = 1; quality_index <= read_length; ++quality_index) quality = quality "I"
            for (read_index = 1; read_index <= read_count; ++read_index) {
                start = target - 120 + int(rand() * 80)
                sequence = ""
                for (offset = 0; offset < read_length; ++offset) {
                    position = start + offset
                    base = reference_base(position)
                    if (position == target && read_index <= alt_count) base = "T"
                    sequence = sequence base
                }
                flag = read_index % 2 == 0 ? 16 : 0
                mapq = 45 + int(rand() * 16)
                printf "%s_read_%03d\t%d\t%s\t%d\t%d\t151M\t*\t0\t0\t%s\t%s\tRG:Z:%s\tMI:Z:%s_molecule_%03d\n", \
                    sample, read_index, flag, contig, start, mapq, sequence, quality, sample, sample, read_index
            }
        }
    ' | samtools view -b - | samtools sort -o "$output" -
    samtools index "$output"
    samtools quickcheck "$output"
}

# Negative control plus two V617F-positive samples in the requested 20-30% range.
generate_bam jak2_negative_120x 120 0 6170 chr9
generate_bam jak2_v617f_25pct_150x 150 38 6171 chr9
generate_bam jak2_v617f_28pct_180x 180 50 6172 9

printf '%s\n' "Generated indexed JAK2 test BAMs in $output_dir"
