#!/usr/bin/env bash
# Run the analysis scripts in their order (analysis/run/steps.tsv: the step, the script,
# its arguments (- for none), the steps it needs, the network or login it needs and the
# time it took).
# Each script runs through analysis/run/run_script.R, which records its provenance in
# results/provenance/; its output goes to results/logs/<step>.log and the progress to
# results/logs/run_all.log. A failed step stops the run.
#
# Usage, from anywhere:
#   bash analysis/run/run_all.sh              every step
#   bash analysis/run/run_all.sh 05a 06g      the steps from 05a to 06g, in the order of the table
#   bash analysis/run/run_all.sh 09d          one step
#   SKIP_NETWORK=1 bash analysis/run/run_all.sh   leave out the steps that download (they
#                                             need the raw data already in data/raw/)
#   DRY_RUN=1 bash analysis/run/run_all.sh    list the steps without running them
# Rscript is taken from $RSCRIPT, else from the PATH.

set -u
cd "$(dirname "$0")/../.." || exit 1
RSCRIPT="${RSCRIPT:-Rscript}"
mkdir -p results/logs
from="${1:-}"
to="${2:-$from}"
[ -z "$from" ] && to=""
on=0; [ -z "$from" ] && on=1
log() { echo "$(date '+%F %T') $*" | tee -a results/logs/run_all.log; }

while IFS=$'\t' read -r step script arguments needs access expected note; do
  [ "$step" = "step" ] && continue
  [ "$step" = "$from" ] && on=1
  if [ "$on" = 1 ]; then
    if [ "${SKIP_NETWORK:-0}" = 1 ] && [ "${access#network}" != "$access" ]; then
      log "skip $step ($script): $access"
    elif [ "${DRY_RUN:-0}" = 1 ]; then
      echo "$step  $script $arguments  (needs: $needs; $access; expected: $expected)"
    else
      log "start $step ($script $arguments; needs: $needs; expected: $expected)"
      # shellcheck disable=SC2086
      [ "$arguments" = "-" ] && arguments=""
      if ! "$RSCRIPT" analysis/run/run_script.R "analysis/$script" $arguments > "results/logs/$step.log" 2>&1; then
        log "FAILED $step: see results/logs/$step.log"
        exit 1
      fi
      log "end $step"
    fi
  fi
  [ -n "$to" ] && [ "$step" = "$to" ] && break
done < analysis/run/steps.tsv
[ "${DRY_RUN:-0}" = 1 ] || log "done"
