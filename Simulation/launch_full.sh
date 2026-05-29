#!/bin/zsh
## Launch the simulation as 3 parallel shards (M4 Pro 14-core:
## 3 shards x 4 chains = 12 cores, 2 spare).
## Run from the Analysis/ directory.
## Args: $1 = scenario list (default all 4)   $2 = SIM_TAG (default full)
SCEN="${1:-S1_jan,S2_jul,S4_miss,S3_nogp}"
TAG="${2:-full}"
NSHARD=3
for s in 0 1 2; do
  SIM_R=200 SIM_SCEN="$SCEN" \
  SIM_CHAINS=4 SIM_ITER=1500 SIM_WARMUP=750 SIM_TAG="$TAG" \
  SIM_NSHARD=$NSHARD SIM_SHARD=$s \
    nohup Rscript Simulation/run_sim.R > "Simulation/log_${TAG}_shard$s.txt" 2>&1 &
  echo "shard $s ($TAG) launched PID $!"
  sleep 3
done
