#!/usr/bin/env bash
# Historical: how stage 7 (amendment protocol-v1.1) was run (moved from results/; analysis/run/run_all.sh replaces it).
cd "$(dirname "$0")/../../.." || exit 1
"${RSCRIPT:-Rscript}" analysis/07c_amendment.R > results/run_07c.log 2>&1 || echo "FAILED" >> results/run_07c.log
echo "07c done" >> results/run_07c.log
