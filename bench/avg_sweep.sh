#!/bin/sh
# Averaging sensitivity: decode probability of N averaged minutes vs S/N.
#   sh bench/avg_sweep.sh [N=4] [trials=20] [snr list]
N=${1:-4}; T=${2:-20}; shift 2 2>/dev/null
SNRS=${*:--24 -25 -26 -27 -28 -29 -30}
D=$(mktemp -d)
echo "averaging $N minutes, $T trials per step"
echo "  S/N dB   decoded"
seed=5000
for snr in $SNRS; do
  ok=0
  for k in $(seq 1 $T); do
    args=""
    for m in $(seq 1 $N); do
      seed=$((seed+1))
      ./pi4gen --text OZ7IGY --snr $snr --f0 688.5 --lead 3.2 --dur 36 --seed $seed $D/m$m.wav
      args="$args $D/m$m.wav 3.0"
    done
    if ./pi4decoder --variants PI4 --avg $args | grep -q '"text": "OZ7IGY"'; then ok=$((ok+1)); fi
  done
  echo "  $snr   $ok/$T"
done
rm -rf $D
