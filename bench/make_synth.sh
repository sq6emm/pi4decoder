#!/bin/sh
# Generate the synthetic benchmark set with pi4gen (deterministic seeds).
# Run inside the dev container from the repo root:  sh bench/make_synth.sh
# Layout (testdata/syn): OZ7IGY, PI4, 12 kHz, 63 s, minute at 3.0 s + dt.
#   s-10_dfX.wav, s-16_dfX.wav : 7 files each, df = -15,-9,-5,0,5,9,15 Hz, dt 0.5 s
#   syn_<snr>_<k>.wav          : snr -20..-32 step 2, k = 0..2, random df/dt
#   spread_<hz>_<k>.wav        : -16 dB, Doppler spread 10..100 Hz
#   drift_<k>.wav              : -16 dB, 2 Hz drift over the message
#   noise_<k>.wav              : 10 pure-noise minutes
set -e
GEN=${GEN:-./pi4gen}
OUT=${OUT:-testdata/syn}
mkdir -p "$OUT"
seed=100
for snr in -10 -16; do
  for df in -15 -9 -5 0 5 9 15; do
    seed=$((seed+1))
    $GEN --text OZ7IGY --snr $snr --f0 $(awk "BEGIN{print 682.8125+$df}") --lead 3.5 --seed $seed "$OUT/s${snr}_df${df}.wav"
  done
done
for snr in -20 -22 -24 -26 -28 -30 -32; do
  for k in 0 1 2; do
    seed=$((seed+1))
    # pseudo-random df in -20..20 Hz and dt in -0.5..1.0 s derived from the seed
    df=$(awk "BEGIN{srand($seed); print int(rand()*41)-20}")
    dt=$(awk "BEGIN{srand($seed+7); printf \"%.2f\", rand()*1.5-0.5}")
    $GEN --text OZ7IGY --snr $snr --f0 $(awk "BEGIN{print 682.8125+$df}") --lead $(awk "BEGIN{print 3.0+$dt}") --seed $seed "$OUT/syn_${snr}_${k}.wav"
  done
done
for sp in 10 20 50 100; do
  for k in 0 1 2; do
    seed=$((seed+1))
    $GEN --text OZ7IGY --snr -16 --spread $sp --lead 3.2 --seed $seed "$OUT/spread_${sp}_${k}.wav"
  done
done
for k in 0 1 2; do
  seed=$((seed+1))
  $GEN --text OZ7IGY --snr -16 --drift 2 --lead 3.1 --seed $seed "$OUT/drift_${k}.wav"
done
for k in 0 1 2 3 4 5 6 7 8 9; do
  seed=$((seed+1))
  $GEN --noise --seed $seed "$OUT/noise_${k}.wav"
done
echo "synthetic set written to $OUT"
