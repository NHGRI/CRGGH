#!/bin/bash
# Generate a swarm file for cellbender remove-background + ptrepack
# Usage: bash cellbender_swarm_generator.sh
# Then submit with:
#   swarm -f cellbender.swarm -g 64 -t 4 --time=48:00:00 \
#         --merge-output --module cellbender \
#         --sbatch "--mail-type=BEGIN,END,FAIL"

BASE_DIR="."
SWARM_OUT="cellbender.swarm"

# Discover samples by looking for */outs directories
SAMPLES=($(ls -d "$BASE_DIR"/*/outs 2>/dev/null | awk -F'/' '{print $(NF-1)}' | sort -u))

> "$SWARM_OUT"  # clear/create output file

for SAMPLE in "${SAMPLES[@]}"; do
    OUTS_DIR="$BASE_DIR/$SAMPLE/outs"
    cat >> "$SWARM_OUT" <<EOF
cd $OUTS_DIR; \\
cellbender remove-background --cpu-threads 4 \\
    --input $OUTS_DIR/raw_feature_bc_matrix.h5 \\
    --output $OUTS_DIR/cb_feature_bc_matrix.h5; \\
ptrepack --complevel 5 cb_feature_bc_matrix_filtered.h5:/matrix cb-seurat_feature_bc_matrix_filtered.h5:/matrix

EOF
done

echo "Swarm file written to: $SWARM_OUT (${#SAMPLES[@]} jobs)"