#!/usr/bin/env bash
# =============================================================================
# GENERALIZED SWEEP v3 — run the whole pipeline, resumably
#
#   ./run_v3_sweep.sh              # run/resume the sweep, then draw everything
#   ./run_v3_sweep.sh figures      # skip the fitting, just redraw
#
# SAFE TO RE-RUN AND SAFE TO INTERRUPT. Every step caches one RDS per cell as
# that cell finishes and skips cells already cached, so a Ctrl-C, a kill, or a
# reboot costs only the cells that were in flight. Re-running after one picks up
# where it stopped rather than starting over.
#
# Step 1 (choose K over K = 1..20) is the long pole — twenty fits per cell
# against steps 2-3's two — so it is sharded across processes. hlme(nproc=)
# already forks internally, so the shards each use ONE thread: N shards x 1
# thread beats N/4 shards x 4 threads, which measured only ~60% CPU per process.
#
# Tunables:
#   V3_SHARDS   how many parallel fitting processes (default: cores - 2)
#   V3_FIGCORES how many parallel drawing processes (default: 10)
#   V3_REFIT=1  ignore the caches and refit everything
# =============================================================================
set -u
cd "$(dirname "$0")"

SHARDS=${V3_SHARDS:-$(( $(nproc) - 2 ))}
[ "$SHARDS" -lt 1 ] && SHARDS=1
FIGCORES=${V3_FIGCORES:-10}

cached () { find "figures/generalized_v3/cache/step$1" -name '*.rds' 2>/dev/null | wc -l; }
stamp  () { date +%H:%M; }

run_sharded () {   # $1 = step number
  local step=$1 s
  echo "[$(stamp)] step $step — $SHARDS shards ($(cached "$step")/72 already cached)"
  for s in $(seq 0 $((SHARDS - 1))); do
    V3_NSHARD=$SHARDS V3_SHARD=$s V3_NPROC=1 \
      Rscript "analyze_all_v3_step$step.R" > ".v3_s${step}_$s.log" 2>&1 &
  done
  wait
  echo "[$(stamp)] step $step done — $(cached "$step")/72 cached"
  # A shard that died leaves its reason in its own log; surface it rather than
  # letting the next step run on a half-built cache.
  if grep -lE "Erreur|Error|Execution halted" ".v3_s${step}_"*.log >/dev/null 2>&1; then
    echo "  !! errors in step $step:"
    grep -hE "Erreur|Error|Execution halted" ".v3_s${step}_"*.log | sort -u | head -5
  fi
}

if [ "${1:-}" != "figures" ]; then
  run_sharded 1
  run_sharded 2
  echo "[$(stamp)] step 3 — LRT + Bonferroni across the plate (single process: the"
  echo "           adjustment spans every compared cell, so it cannot be sharded)"
  Rscript analyze_all_v3_step3.R 2>&1 | tee .v3_s3.log | tail -6
fi

echo "[$(stamp)] figures — $FIGCORES cores"
V3_FIGCORES=$FIGCORES Rscript analyze_all_v3_figures.R 2>&1 | tee .v3_figs.log | tail -6
echo "[$(stamp)] DONE — figures/generalized_v3/"
