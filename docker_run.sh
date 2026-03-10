#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <path-to-data-directory>"
  exit 1
fi

DATA_DIR="$(cd "$1" && pwd)"

# Create output directories inside the data directory
mkdir -p "$DATA_DIR/processed data"
mkdir -p "$DATA_DIR/processed data without artifact correction (for benchmarking)"
mkdir -p "$DATA_DIR/results/biomarker results"
mkdir -p "$DATA_DIR/results/biomarker results for TK:RTK vignette"

# Script 1 — data processing
docker run --rm \
  -v "$DATA_DIR:/app/data" \
  -v "$DATA_DIR/results:/app/results" \
  oncology-reference "scripts/1 - DATA_PROCESSING.R"

# Script 2 — uncorrected data processing for benchmarking
docker run --rm \
  -v "$DATA_DIR:/app/data" \
  -v "$DATA_DIR/results:/app/results" \
  oncology-reference "scripts/2 - DATA_PROCESSING - UNCORRECTED FILES FOR COMPARISON.R"

# Script 3 — biomarker table generation
docker run --rm \
  -v "$DATA_DIR:/app/data" \
  -v "$DATA_DIR/results:/app/results" \
  oncology-reference "scripts/3 - GENERATE BIOMARKER TABLES.R"
