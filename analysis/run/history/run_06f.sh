#!/usr/bin/env bash
# Historical: how stage 6 (D18) was run (moved from results/; analysis/run/run_all.sh replaces it).
cd "$(dirname "$0")/../../.." || exit 1
"${RSCRIPT:-Rscript}" analysis/06f_hdbscan_pca.R > results/run_06f.log 2>&1 || echo "FAILED" >> results/run_06f.log
echo "06f done" >> results/run_06f.log
