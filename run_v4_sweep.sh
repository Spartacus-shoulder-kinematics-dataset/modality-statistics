#!/usr/bin/env bash
# =============================================================================
# GENERALIZED SWEEP v4 — run the pipeline, resumably
#
#   ./run_v4_sweep.sh              # step 1 -> 2 -> 3 -> figures, resuming caches
#   ./run_v4_sweep.sh step1        # step 1 + figures only. Stop here to look at
#                                  # cells/*/01_choosing_k and edit selected_K in
#                                  # figures/generalized_v4/K_selection.csv
#   ./run_v4_sweep.sh figures      # redraw only
#
# SAFE TO RE-RUN AND SAFE TO INTERRUPT. Every step caches one RDS per cell as
# that cell finishes and skips cells already cached, so a Ctrl-C, a kill, or a
# reboot costs only the cells in flight. Editing selected_K and re-running the
# full pipeline refits exactly the cells whose K changed.
#
# Step 1 (K = 1..10) is sharded across processes. hlme(nproc=) already forks
# internally, so each shard uses ONE thread.
#
# Tunables:
#   V4_SHARDS   parallel fitting processes (default: cores - 2)
#   V4_FIGCORES parallel drawing processes (default: 10)
#   V4_REFIT=1  ignore the caches and refit everything
# =============================================================================
set -u
cd "$(dirname "$0")"

SHARDS=${V4_SHARDS:-$(( $(nproc) - 2 ))}
[ "$SHARDS" -lt 1 ] && SHARDS=1
FIGCORES=${V4_FIGCORES:-10}

cached () { find "figures/generalized_v4/cache/step$1" -name '*.rds' 2>/dev/null | wc -l; }
stamp  () { date +%H:%M; }

run_sharded () {   # $1 = step number
  local step=$1 s
  echo "[$(stamp)] step $step — $SHARDS shards ($(cached "$step")/72 already cached)"
  for s in $(seq 0 $((SHARDS - 1))); do
    V4_NSHARD=$SHARDS V4_SHARD=$s V4_NPROC=1 \
      Rscript "analyze_all_v4_step$step.R" > ".v4_s${step}_$s.log" 2>&1 &
  done
  wait
  echo "[$(stamp)] step $step done — $(cached "$step")/72 cached"
  # A shard that died leaves its reason in its own log; surface it rather than
  # letting the next step run on a half-built cache.
  if grep -lE "Erreur|Error|Execution halted" ".v4_s${step}_"*.log >/dev/null 2>&1; then
    echo "  !! errors in step $step:"
    grep -hE "Erreur|Error|Execution halted" ".v4_s${step}_"*.log | sort -u | head -5
  fi
}

# Shards finish in any order and each writes the step-1 tables from whatever is
# cached at that moment; one final pass after all of them guarantees the files
# describe every cell.
step1_tables () {
  Rscript -e 'source("analyze_all_v4_common.R"); write_step1_tables()' > /dev/null
  echo "[$(stamp)] K_selection.csv written — figures/generalized_v4/K_selection.csv"
}

step3 () {
  echo "[$(stamp)] step 3 — LRT + Bonferroni across the plate (single process)"
  Rscript analyze_all_v4_step3.R 2>&1 | tee .v4_s3.log | tail -8
}

figures () {
  echo "[$(stamp)] figures — $FIGCORES cores"
  V4_FIGCORES=$FIGCORES Rscript analyze_all_v4_figures.R 2>&1 | tee .v4_figs.log | tail -6
  echo "[$(stamp)] DONE — figures/generalized_v4/"
}

case "${1:-all}" in
  step1)
    run_sharded 1; step1_tables; figures
    echo "Now inspect cells/*/01_choosing_k, edit selected_K, then run: ./run_v4_sweep.sh" ;;
  figures)
    figures ;;
  all)
    run_sharded 1; step1_tables; run_sharded 2; step3; figures ;;
  *)
    echo "usage: $0 [all|step1|figures]"; exit 2 ;;
esac
