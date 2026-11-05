#!/usr/bin/env bash
# Historical: how stage 6 (D18 and the protocol tests) was run (moved from results/; analysis/run/run_all.sh replaces it).
cd "$(dirname "$0")/../../.." || exit 1
R="${RSCRIPT:-Rscript}"
"$R" analysis/06f_hdbscan_pca.R > results/run_06f.log 2>&1 || echo "FAILED 06f" >> results/run_06f.log
"$R" analysis/06g_protocol_tests.R > results/run_06g.log 2>&1 || echo "FAILED 06g" >> results/run_06f.log
echo "06fg done" >> results/run_06f.log
