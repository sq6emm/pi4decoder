#!/bin/sh
# Reproduce the pi4decoder benchmark tables. Runs everything in Docker.
#
#   sh bench/benchmark.sh            # pi4decoder on all available test sets
#   REFS=1 sh bench/benchmark.sh     # also the reference decoders (if present)
#
# Test sets (set the variables to your copies; missing sets are skipped):
#   synthetic set  : generated into testdata/syn by bench/make_synth.sh
#   OZ2M_DIR       : OZ2M's off-air PI4 recordings, 6 kHz WAV
#                    (https://rudius.net/oz2m/ngnb/pi4_cwid_car_files.zip, unzipped)
#   REC_DIR        : your own 12 kHz minute files (minute at 3.0 s), e.g.
#                    13 noise-only minutes and their *_inj16.wav versions
# Reference decoders (REFS=1): PI4DEC_IMAGE (image with /usr/local/bin/pi4dec),
#   PI4D_BIN, MSHV_BIN (host paths, run in BUILD_IMAGE).
set -e
cd "$(dirname "$0")/.."
IMG=${IMG:-pi4decoder-dev}
OUT=bench/results
mkdir -p "$OUT"
RUN="docker run --rm --user $(id -u):$(id -g) --label owner=claude -v $(pwd):/work -w /work"
[ -n "$OZ2M_DIR" ] && RUN="$RUN -v $OZ2M_DIR:/oz2m:ro"
[ -n "$REC_DIR" ] && RUN="$RUN -v $REC_DIR:/rec:ro"

$RUN --name pi4decoder-bench-build $IMG sh -c 'make -s pi4decoder pi4gen tests/test_noise >/dev/null && gcc -O2 -o bench/rusage bench/rusage.c && gcc -O2 -o bench/wavcut bench/wavcut.c src/wav.c'
[ -d testdata/syn ] || $RUN --name pi4decoder-bench-synth $IMG sh bench/make_synth.sh

echo "== pi4decoder: synthetic set" | tee "$OUT/pi4decoder_syn.txt"
$RUN --name pi4decoder-bench-syn $IMG sh -c 'for f in testdata/syn/*.wav; do echo "== $(basename $f)"; ./bench/rusage ./pi4decoder --variants PI4 $f 3.0 2>&1; done' >> "$OUT/pi4decoder_syn.txt"

if [ -n "$OZ2M_DIR" ]; then
  echo "== pi4decoder: OZ2M scan" | tee "$OUT/pi4decoder_oz2m.txt"
  $RUN --name pi4decoder-bench-oz2m $IMG sh -c 'for f in /oz2m/*.wav; do echo "== $(basename $f)"; ./bench/rusage ./pi4decoder --scan $f 2>&1; done' >> "$OUT/pi4decoder_oz2m.txt"
fi
if [ -n "$REC_DIR" ]; then
  echo "== pi4decoder: own minutes (all variants, +-150 Hz)" | tee "$OUT/pi4decoder_rec.txt"
  $RUN --name pi4decoder-bench-rec $IMG sh -c 'for f in /rec/*.wav; do echo "== $(basename $f)"; ./bench/rusage ./pi4decoder $f 3.0 2>&1; done' >> "$OUT/pi4decoder_rec.txt"
fi

echo "== false decodes on pure noise (200 minutes)" | tee "$OUT/pi4decoder_noise.txt"
$RUN --name pi4decoder-bench-noise $IMG ./tests/test_noise 200 42 >> "$OUT/pi4decoder_noise.txt" || true

if [ "$REFS" = 1 ]; then
  PI4DEC_IMAGE=${PI4DEC_IMAGE:-rsanttracker-beacon}
  BUILD_IMAGE=${BUILD_IMAGE:-rsanttracker-mshvbuild}
  echo "== pi4dec: synthetic set" | tee "$OUT/pi4dec_syn.txt"
  docker run -i --rm --user $(id -u):$(id -g) --label owner=claude --name pi4decoder-bench-pi4dec -v $(pwd):/work -w /work --entrypoint python3 $PI4DEC_IMAGE - <<'PY' >> "$OUT/pi4dec_syn.txt"
import glob, resource, subprocess
for f in sorted(glob.glob("testdata/syn/*.wav")):
    print("==", f.split("/")[-1])
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    out = subprocess.run(["/usr/local/bin/pi4dec", "--variants", "PI4", f, "3.0"], capture_output=True, text=True).stdout
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    print(out, end="")
    print(f"RUSAGE user {r1.ru_utime - r0.ru_utime:.2f} s maxrss {r1.ru_maxrss} kB")
PY
  if [ -n "$PI4D_BIN" ]; then
    echo "== pi4d: synthetic set (file re-cut to start at the minute)" | tee "$OUT/pi4d_syn.txt"
    $RUN --name pi4decoder-bench-cut $IMG sh -c 'mkdir -p testdata/syn_m0; for f in testdata/syn/*.wav; do ./bench/wavcut $f testdata/syn_m0/$(basename $f) 3.0; done'
    docker run --rm --user $(id -u):$(id -g) --label owner=claude --name pi4decoder-bench-pi4d -v $(pwd):/work -v $PI4D_BIN:/pi4d:ro -w /work $BUILD_IMAGE bash -c '
      for f in testdata/syn_m0/*.wav; do echo "== $(basename $f)"; ./bench/rusage /pi4d $f 2>&1; done' >> "$OUT/pi4d_syn.txt"
  fi
  if [ -n "$MSHV_BIN" ]; then
    echo "== mshvpi4 --tol 110: synthetic set" | tee "$OUT/mshv_syn.txt"
    docker run --rm --user $(id -u):$(id -g) --label owner=claude --name pi4decoder-bench-mshv --ulimit stack=1073741824:1073741824 -v $(pwd):/work -v $MSHV_BIN:/mshvpi4:ro -w /work $BUILD_IMAGE bash -c '
      for f in testdata/syn/*.wav; do echo "== $(basename $f)"; ./bench/rusage /mshvpi4 --tol 110 $f 2>&1; done' >> "$OUT/mshv_syn.txt"
  fi
fi
$RUN --name pi4decoder-bench-table $IMG python3 bench/table.py
