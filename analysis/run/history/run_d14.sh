#!/usr/bin/env bash
# Historical: how stages 4 to 6 after D14 were run (moved from results/; analysis/run/run_all.sh replaces it).
# D14 rerun: finish stage 4 (resumable), then stages 5a to 6d (results/run_g.sh).
R="${RSCRIPT:-Rscript}"
cd "$(dirname "$0")/../../.." || exit 1
echo "$(date '+%F %T') start 04_grouping" >> results/run_g.log
"$R" analysis/04_grouping.R > results/grouping_run_d14c.log 2>&1 || echo "$(date '+%F %T') FAILED 04_grouping" >> results/run_g.log
echo "$(date '+%F %T') end 04_grouping" >> results/run_g.log
bash analysis/run/history/run_g.sh
echo "$(date '+%F %T') all done" >> results/run_g.log
