#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
       echo "Usage: $0 <path-to-data-directory> [--cpus N] [--memory SIZE]"
     exit 1
     fi
     
     # Convert the provided path to an absolute path
     DATA_DIR="$(cd "$1" && pwd)"
     shift
     
     DOCKER_OPTS=""
     while [ $# -gt 0 ]; do
               case "$1" in
             --cpus) DOCKER_OPTS="$DOCKER_OPTS --cpus=$2"; shift 2 ;;
             --memory) DOCKER_OPTS="$DOCKER_OPTS --memory=$2"; shift 2 ;;
             *) echo "Unknown option: $1"; exit 1 ;;
esac
done

echo "Setting up output directories in $DATA_DIR..."

# Create output directories inside the local data directory
# This ensures they are owned by your user, not the Docker root user
mkdir -p "$DATA_DIR/processed_data"
mkdir -p "$DATA_DIR/processed_data_without_artifact_correction"
mkdir -p "$DATA_DIR/biomarker_results"
mkdir -p "$DATA_DIR/biomarker_results_for_TK_RTK_vignette"
mkdir -p "$DATA_DIR/release_files"




# Note: The single mount -v "$DATA_DIR:/app/data" covers everything, 
# because your R scripts use relative paths starting with "data/..."

echo "Running Pipeline..."

# Script 1 — data processing
docker run --rm $DOCKER_OPTS \
-v "$DATA_DIR:/app/data" \
oncology-reference "scripts/1 - DATA_PROCESSING.R"

# Script 2 — uncorrected data processing for benchmarking
docker run --rm $DOCKER_OPTS \
-v "$DATA_DIR:/app/data" \
oncology-reference "scripts/2 - DATA_PROCESSING - UNCORRECTED FILES FOR COMPARISON.R"

# Script 3 — biomarker table generation
docker run --rm $DOCKER_OPTS \
-v "$DATA_DIR:/app/data" \
oncology-reference "scripts/3 - GENERATE BIOMARKER TABLES.R"

# Script 4 — depmap release files
docker run --rm $DOCKER_OPTS \
-v "$DATA_DIR:/app/data" \
oncology-reference "scripts/4 - DEPMAP RELEASE FILES.R"

echo "Pipeline completed successfully!"