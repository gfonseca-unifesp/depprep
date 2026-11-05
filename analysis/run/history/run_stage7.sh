#!/usr/bin/env bash
# Historical: how stage 7 was run (moved from results/; analysis/run/run_all.sh replaces it).
# Stage 7 of protocol-v1: the full grid, then H2a on new seeds.
R="${RSCRIPT:-Rscript}"
cd "$(dirname "$0")/../../.." || exit 1
for s in 07_full_grid; do
  echo "$(date '+%F %T') start $s" >> results/run_stage7.log
  "$R" analysis/$s.R > results/run_stage7_$s.log 2>&1 || echo "$(date '+%F %T') FAILED $s" >> results/run_stage7.log
  echo "$(date '+%F %T') end $s" >> results/run_stage7.log
done
echo "$(date '+%F %T') stage 7 done" >> results/run_stage7.log
