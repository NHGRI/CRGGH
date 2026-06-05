#!/bin/bash
# Generate a swarm file for cellranger count
# Usage: bash cellranger_count_swarm_generator.sh
# Then submit with:
#   swarm -f cellranger_count.swarm -g 64 -t 12 --time=48:00:00 \
#         --merge-output --module cellranger \
#         --sbatch "--mail-type=BEGIN,END,FAIL"

FASTQ_PATH="./fastqfiles"
TRANSCRIPTOME="$CELLRANGER_REF/refdata-gex-GRCh38-2024-A"
SWARM_OUT="cellranger_count.swarm"

# Discover samples from FASTQ directory (looks for unique sample-name prefixes)
SAMPLES=($(ls "$FASTQ_PATH" | grep -oP '^[^_]+' | sort -u))

> "$SWARM_OUT"  # clear/create output file

for SAMPLE in "${SAMPLES[@]}"; do
    cat >> "$SWARM_OUT" <<EOF
FASTQ_PATH=$FASTQ_PATH; \\
ulimit -u 10240 -n 16384; \\
cellranger count --id=$SAMPLE \\
     --transcriptome=$TRANSCRIPTOME \\
     --fastqs="\$FASTQ_PATH" \\
     --create-bam=true \\
     --sample=$SAMPLE \\
     --localcores=\$SLURM_CPUS_PER_TASK \\
     --localmem=34 \\
     --jobmode=slurm \\
     --maxjobs=10

EOF
done

echo "Swarm file written to: $SWARM_OUT (${#SAMPLES[@]} jobs)"