#!/usr/bin/env bash
# Historical: how stages 5 and 6 after D14 were run (moved from results/; analysis/run/run_all.sh replaces it).
# D14 rerun, after stages 3 and 4: 5a, 5b, 5c, 6 and the stage 6 extras, in order.
R="${RSCRIPT:-Rscript}"
cd "$(dirname "$0")/../../.." || exit 1
for s in 05a_h2_mechanism 05b_correction 05c_whitening_st 06_pilot 06b_pilot_absorption 06c_q3_warning 06d_whitening_filters; do
  echo "$(date '+%F %T') start $s" >> results/run_g.log
  "$R" analysis/$s.R > results/run_g_$s.log 2>&1 || echo "$(date '+%F %T') FAILED $s" >> results/run_g.log
  echo "$(date '+%F %T') end $s" >> results/run_g.log
done
