#!/usr/bin/env bash
# Historical: how stages 5 and 6 after the fix a9a609e were run (moved from results/; analysis/run/run_all.sh replaces it).
# D14 rerun after the fix of the variable union (a9a609e): stages 5b and 6.
R="${RSCRIPT:-Rscript}"
cd "$(dirname "$0")/../../.." || exit 1
for s in 05b_correction 06_pilot; do
  echo "$(date '+%F %T') start $s" >> results/run_g.log
  "$R" analysis/$s.R > results/run_g_$s.log 2>&1 || echo "$(date '+%F %T') FAILED $s" >> results/run_g.log
  echo "$(date '+%F %T') end $s" >> results/run_g.log
done
echo "$(date '+%F %T') fix run done" >> results/run_g.log
