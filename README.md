# pi4decoder

A small, dependency-free decoder for **PI4**, the 4-FSK beacon mode of the
IARU Region 1 "Next Generation Beacon" platform (OZ2M / DJ5HG,
<https://rudius.net/oz2m/ngnb/pi4.htm>).

* plain C99, one library (`libpi4.a`, four source files, ~1700 lines) and a
  command-line tool; needs only `libm` (no FFTW, no Qt, no Fortran)
* 35 kB static binary, ~15 MB RAM and ~0.2-0.5 s of CPU per minute of audio
* reads 6 / 8 / 12 / 24 / 48 kHz mono WAV (PCM 8/16/24/32 bit or float)
* all four tone-spacing variants (PI4, PI4-80, PI4-96, PI4-120)
* Doppler-spread / aircraft-scatter / rain-scatter signals (spread-matched
  energy windows up to ~100 Hz), slow drift, AGC holes and dropouts
* zero false decodes measured on 3000 minutes of pure noise plus 13 real
  noise-only minutes, without sacrificing sensitivity (see Benchmarks)
* blind decoding to about -23 dB S/N (2500 Hz), a known-message mode
  (`--calls`) to about -27 dB, and multi-minute averaging (`--avg`)
* JSON-lines output, CLI compatible with the `pi4dec` call convention

License: GPL-3.0-or-later.

## Why C

The goal was "outperform the others and be the lightest possible". Plain C99
gives the smallest binary and the fewest moving parts: it links into
WSJT-X/MSHV-style C++ hosts directly, builds with every microcontroller
toolchain (the core is allocation-light and could be given a static arena),
and has no runtime. Rust with a no_std core would have been a fine second
choice, but it adds a toolchain and a larger binary for no decoding benefit.

## Building

```sh
make                 # pi4decoder, pi4gen, libpi4.a  (gcc or clang, -lm)
make test            # unit tests: spec round trip, decoder, noise
make small           # -Os static-ish build for size measurement
```

Everything in this repository was built and tested inside Docker; the
`Dockerfile` is the image used (`debian:bookworm-slim` + gcc + python3 for the
benchmark scripts):

```sh
docker build -t pi4decoder-dev .
./run.sh make test
```

## Usage

```
pi4decoder [options] file.wav [boundary_s] [file2.wav boundary2_s ...]
  boundary_s      where the even UTC minute (PI4 start) falls in the file (default 0)
  --radius HZ     tone-0 search +-HZ around the convention (default 150)
  --center HZ     CW carrier convention: tone0 = center - spacing/2 (default 800)
  --band LO HI    search tone 0 in LO..HI Hz instead of center/radius
  --variants L    comma list of PI4,PI4-80,PI4-96,PI4-120 (default all)
  --dt S          message start searched +-S s around the boundary (default 2.5)
  --scan          search the whole file for any message start (boundary ignored;
                  default band 300..2700 Hz unless --band/--radius given)
  --spread N      widest tone energy window: 0=2.9 Hz .. 5=97 Hz (default 5)
  --calls LIST    also match these known messages, reported with listed=true
  --avg           average all given files (same alignment) before decoding
  --gate Z        blind-decode acceptance gate on fit_z (default 6.0)
  --relax         accept marginal decodes one unit below the gate
  --fs HZ         override the WAV sample rate (clock-error experiments)
  --verbose       diagnostics on stderr (twice: per-symbol dump)
```

Typical monitor call (as the `pi4dec` convention: 12 kHz file, minute at
3.0 s, OZ7IGY on the 800 Hz CW convention):

```sh
pi4decoder --variants PI4 --radius 50 minute.wav 3.0
{"text": "OZ7IGY", "variant": "PI4", "tone0_hz": 682.81, "dt_s": -0.001, "t_s": 2.999,
 "snr_db": -16.2, "fit": 1.000, "fit_z": 92.1, "sync_z": 50.8, "hard_errors": 1,
 "spread_hz": 2.9, "plausible": true, "listed": false, "edge": false}
```

Fields: `text` message (trailing spaces trimmed); `variant`; `tone0_hz` audio
frequency of tone 0; `dt_s` start relative to `boundary_s`; `t_s` absolute
start in the file; `snr_db` S/N in 2500 Hz; `fit` soft agreement of the
re-encoded message with the received symbols (-1..1); `fit_z` the same as a
z-score (the acceptance statistic); `sync_z` sync-correlation z-score of the
cell; `hard_errors` symbols whose data bit disagrees; `spread_hz` energy
window used; `plausible` text looks like a callsign or `/status`; `listed`
found by the known-message matcher; `edge` tone 0 of this variant lies outside
the requested band (`--band`, or `--center` - spacing/2 ± `--radius`) by more
than half a bin. A signal just outside the band is still caught by the wide
spread windows at the band edge; it is then re-centred to its true tone 0
(so `tone0_hz`, `spread_hz` and `snr_db` are right) and kept, with
`edge: true`.

Scan a long recording of unknown timing and frequency (e.g. OZ2M's 120 s
6 kHz files):

```sh
pi4decoder --scan 432M_130815_190746_261.wav
```

Known-message mode and averaging (a monitor that knows which beacon it is
listening to):

```sh
pi4decoder --variants PI4 --calls OZ7IGY,SR6LEG minute.wav 3.0
pi4decoder --variants PI4 --avg m1.wav 3.0 m2.wav 3.0 m3.wav 3.0 m4.wav 3.0
```

Synthetic test signals:

```sh
pi4gen --text OZ7IGY --snr -22 --f0 690 --lead 3.4 --spread 30 --drift 1 out.wav
pi4gen --noise --seed 7 noise.wav
```

### Library

`src/pi4.h` is the whole API: `pi4_encode` / `pi4_pack` / `pi4_unpack`
(spec), `pi4_fano` (sequential decoder), `pi4_decode` (audio in, results out),
`pi4_spec_build` / `pi4_spec_add` / `pi4_decode_spec` (spectrogram level, for
averaging), `pi4_resample12k`, `pi4_text_plausible`.

## Method

1. **Resample** to 12 kHz (windowed-sinc polyphase; 6/8/24/48 kHz in).
2. **Spectrogram**: one 2000-sample (one symbol, rectangular = matched
   filter) window every 250 samples (1/8 symbol), zero-padded to a 4096-point
   real FFT (2.93 Hz bins, so a tone spacing is exactly 2K bins). Only the
   band of interest is kept. Each bin is **whitened** by a robust noise floor
   (30th percentile over time, median over 13 bins) and capped at +40 dB.
3. **Sync search**: PI4's 146-bit sync vector is the low bit of every symbol,
   so for each (time, tone-0 bin, variant) the statistic
   sum_n sign_n * ((P0+P2) - (P1+P3)) is computed on three spread windows
   (1, 5 and 17 bins). One lookup per symbol per cell, z-normalised by the
   median/MAD of all cells, local maxima above z = 4 become candidates. A
   quota per hypothesis and a "same frequency within 30 s" rule keep one
   interferer from crowding out a real signal. Cost: ~50 ms for a +-2.5 s,
   +-150 Hz, four-variant minute search.
4. **Refinement**: per candidate the best of six windows (1..33 bins) by a
   self-normalised sync statistic; narrow winners are refined on the audio
   with Goertzel DFTs at exact frequency (+-1/2 bin in 1/4 bins, +-1/2 hop in
   1/4 hops).
5. **Soft metrics**: per symbol only the two tones the sync bit allows can
   carry the message; the other two measure the noise. LLR models:
   Rayleigh-prior (`(P1-P0)*rho/(1+rho)`, amplitude unknown per symbol -
   robust to fading, AGC and dropouts), exact noncoherent Bessel
   (`ln I0(2 sqrt(rho P1)) - ln I0(2 sqrt(rho P0))`) and Gaussian for wide
   windows / averaged spectrograms. Fano metrics are the standard
   `log2(2/(1+e^-L)) - 1/2` in 1/32 bit.
6. **Sequential decoding**: a Fano decoder over the rate-1/2 K=32 code
   (73 nodes, zero tail), run forward and over the time-reversed code
   (bit-reversed polynomials; bidirectional attempts fail on different noise
   bursts), with LLR scalings 1, 0.5, 2, 0.25. Effort is scaled by the
   candidate's sync z (weak cells get a quarter budget).
7. **Gates**: the decoded message is re-encoded and correlated with the soft
   symbols (`fit_z`, a z-score against the measured noise); accepted at
   `fit_z >= 6` for callsign-like texts, 9 otherwise; the all-zero codeword is
   refused. Overlapping decodes (same time, overlapping tone spans) are
   reduced to the one explaining the most signal energy per sqrt(window).
8. **Known-message mode** correlates each listed message's 146 tones directly
   with the spectrogram (a matched filter; gate z >= 6.5). **Averaging** sums
   whitened spectrograms of consecutive minutes before step 3.

PI4 has no CRC, so points 3, 7 and the text plausibility check are what keep
noise from turning into "callsigns".

## Limits

* Blind decoding stops at the sequential-decoding cutoff rate: for
  noncoherent 2-of-4 FSK at 6 baud, R0 = 1/2 falls at about -23.8 dB S/N
  (2500 Hz). Measured 50 % point: about -23.3 dB. Beyond that only the
  known-message matcher (-27 dB) or averaging (+3 dB per doubling) help.
* Spread signals cost sensitivity in proportion to the window (see table).
* Scan mode of a 120 s file uses ~120 MB (the whole spectrogram is kept);
  minute mode ~15 MB.
* OZ2M's `pi4_cwid_car.mp4` audio (tones ~238 Hz apart and, as it turns out,
  symbol clock scaled too) is not decoded, nor with a sample-rate override.

## Benchmarks

All measurements were made in this repository's Docker image on one core
of a heavily loaded host (CPU seconds are user time). Reference decoders:
`pi4dec` (Rust, sdroxide/rsanttracker), `pi4d` (C port, WSJT-X PR #81),
MSHV 2.x `mshvpi4 --tol 110` (LZ2HV, JT4-derived), PI-RX (OZ1CKG, Windows;
figures from its published results on the same recordings, documented 50 %
limit -22/-23 dB).

### Files with a correct decode

Test material: OZ2M's off-air recordings
(<https://rudius.net/oz2m/ngnb/pi4_cwid_car_files.zip>, 6 kHz, 120 s, timing
and frequency unknown: decoded with `--scan`); 13 of our own 23 cm minutes
containing no PI4 at all, and the same minutes with a -16 dB OZ7IGY injected
at 682.8 Hz (`pi4decoder file.wav 3.0`, all variants, +-150 Hz, the monitor's
configuration); the synthetic sets of the original comparison (reference
columns measured by the parent project with the same files).

| Test | Files | pi4dec | pi4d | MSHV | PI-RX | **pi4decoder** |
|---|---|---|---|---|---|---|
| Synthetic -10 dB | 7 | 7 | 7 | 7 | 7 | **7** |
| Synthetic -16 dB | 7 | 7 | 7 | 7 | 7 | **7** |
| Synthetic -20 dB | 3 | 3 | 3 | 0 | 3 | **3** |
| Synthetic -24..-32 dB | 15 | 0 | 0 | 0 | 0 | **1** (one -24 dB file) |
| OZ2M 10 GHz ED5YAE/B (tones spread 30-80 Hz) | 2 | 0 | 0 | 2 | 2 | **2** |
| OZ2M 28 MHz G3ZJO | 1 | 1 | 1 | 1 | 0 | **1** |
| OZ2M 432 MHz OZ7IGY 2013 | 16 | 14 | 14 | 15 | 15 | **16** |
| OZ2M 70 MHz IW9GDC/B | 2 | 2 | 2 | 2 | 0 | **2** |
| OZ2M 70 MHz OZ7IGY 2014 | 8 | 8 | 8 | 7 | 8 | **8** |
| OZ2M 70 MHz OZ7IGY 2015 (weak, drifting) | 4 | 4 | 4 | 4 | 1 | **4** |
| OZ2M total | 33 | 29 | 29 | 31 | 26 | **33** |
| OZ2M pi4_cwid_car.mp3 audio | 1 | 1 | 1 | 1 | 1 | **1** |
| OZ2M pi4_cwid_car.mp4 audio (scaled) | 1 | 0 | 0 | 0 | 0 | 0 |
| Our minutes (no signal) | 13 | 0 | 0 | 0 | 0 | **0** |
| Our minutes + -16 dB OZ7IGY | 13 | 13 | 13 | 13 | 13 | **13** |
| Garbage decodes, all of the above | - | 1 | 1 | 1 | 0 | **0** |

The one OZ2M file every other decoder misses (432M_130815_185710_897) has
AGC/dropout holes in a quarter of its symbols; the Rayleigh metric treats
them as erasures. The two 10 GHz files decode with 9 and 15 Hz windows, the
2015 70 MHz files (receiver drifting tens of Hz within a message) with
50-100 Hz windows.

### Sensitivity (synthetic OZ7IGY, AWGN, random df +-20 Hz and dt, 50 trials per step)

`bench/sensitivity.py --trials 50 --snr=-18:-28`:

| S/N (2500 Hz) | -18 | -19 | -20 | -21 | -22 | -23 | -24 | -25 | -26 |
|---|---|---|---|---|---|---|---|---|---|
| decoded | 50/50 | 50/50 | 50/50 | 50/50 | 48/50 | 40/50 | 7/50 | 0/50 | 0/50 |
| wrong texts | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

50 % point: **-23.3 dB** (PI-RX documents -22/-23 dB; pi4dec/pi4d decode
the same -22 dB files and the same single -24 dB file as this decoder in the
regenerated set below, MSHV stops at about -19 dB). This is the
sequential-decoding cutoff-rate limit for noncoherent 2-of-4 FSK at 6 baud
(R0 = 1/2 at -23.8 dB); a bigger cycle budget, the reversed-code attempts or
other LLR scalings do not move it. Reported `snr_db` tracks the true value
within 0.3 dB down to -23 dB.

Other channels (30 trials per step, 50 % points):

| channel | 50 % point |
|---|---|
| Doppler spread 10 Hz | -20.7 dB |
| Doppler spread 20 Hz | -19.7 dB |
| Doppler spread 50 Hz | -18.7 dB |
| Doppler spread 100 Hz | -17.3 dB |
| linear drift 2 Hz over the message | -23.3 dB (no loss) |
| linear drift 5 Hz over the message | -22.9 dB |
| known-message mode (`--calls`), AWGN | -27.5 dB |
| 4-minute averaging (`--avg`), AWGN | -27 dB |

### False decodes

* 4200 minutes of pure Gaussian noise (`tests/test_noise`, all four
  variants, +-150 Hz, gate disabled): **0 Fano convergences**, hence 0
  decodes at any gate. Noise cells never reach the sequential decoder's
  threshold; the `fit_z >= 6` gate is a second line of defence that was
  never needed.
* 13 real noise-only minutes, 33 OZ2M scans (120 s each, 300-2700 Hz), 10
  synthetic noise minutes: 0 false decodes.
* Known-message mode: gate z >= 6.5 on a ~25000-cell matched-filter search
  per message and window (Gaussian tail ~1e-5 per minute per message); no
  false "SR6LEG" in 240 list-mode trials.

### Size, speed, memory (per 63 s minute file, one core)

| decoder | binary | dependencies | CPU s / minute | peak RSS |
|---|---|---|---|---|
| **pi4decoder** | 35 kB (static -Os), 56 kB (-O3) | libm | 0.10-0.15 (PI4 only), 0.25-0.45 (all variants, real noise) | 10 MB |
| pi4dec (Rust) | 2.2 MB | rustfft, rayon, mfsk-core | 0.6-0.8 | 11 MB |
| pi4d (C) | 28 kB | libfftw3 | 0.6-1.0 | 7 MB |
| MSHV mshvpi4 | 1.26 MB | Qt5, fftw3 | 1.3-1.6 | 243 MB |

A 120 s scan (`--scan`, 2400 Hz band, 24 candidates) takes 3-7 s and
60-95 MB.

### Regenerated synthetic set, all decoders run here (`bench/benchmark.sh REFS=1`)

`bench/make_synth.sh` regenerates the set with `pi4gen` (OZ7IGY, PI4, random
df/dt, fixed seeds), so the reference columns below were measured by this
repository's scripts rather than copied:

| Test | Files | pi4decoder | pi4dec | pi4d | MSHV |
|---|---|---|---|---|---|
| Synthetic -10 dB | 7 | 7 | 7 | 7 | 7 |
| Synthetic -16 dB | 7 | 7 | 7 | 7 | 7 |
| Synthetic -20/-22 dB | 6 | 6 | 6 | 6 | 0 |
| Synthetic -24..-32 dB | 15 | 1 | 1 | 1 | 0 |
| Synthetic drift 2 Hz, -16 dB | 3 | 3 | 3 | 3 | 3 |
| Synthetic spread 10-100 Hz, -16 dB | 12 | 12 | 7 | 7 | 12 |
| Synthetic pure noise (false decodes) | 10 | 0 | 0 | 0 | 0 |
| mean CPU s / file | | 0.12 | 0.70 | 0.87 | 1.38 |
| peak RSS | | 7.6 MB | 11 MB | 7.7 MB | 244 MB |

`bench/results/` holds the raw outputs; `python3 bench/table.py` rebuilds
the tables.

## Repository

```
src/pi4.h          public API
src/pi4_spec.c     alphabet, packing, convolutional code, interleaver, sync, text plausibility
src/pi4_fano.c     sequential decoder (forward / reversed code)
src/pi4_fft.c      radix-2 FFT and real-FFT wrapper
src/pi4_decode.c   resampler, spectrogram, search, refinement, metrics, gates
src/main.c         CLI        src/pi4gen.c  signal generator      src/wav.c  WAV I/O
tests/             test_spec (round trip vs OZ2M encoder, Fano, plausibility),
                   test_decode (FFT, resampler, end-to-end at 6/8/12/24/48 kHz,
                   noise, list mode, averaging), test_noise (false-decode corpus)
bench/             sensitivity.py, avg_sweep.sh, make_synth.sh, benchmark.sh, table.py
```

Test WAVs are not committed (`.gitignore`); the synthetic ones are
regenerated by `bench/make_synth.sh`, OZ2M's recordings are third-party
(<https://rudius.net/oz2m/ngnb/pi4_cwid_car_files.zip>).

## Acknowledgements

PI4 is by Bo, OZ2M, and Klaus, DJ5HG (sync vector). The reference encoder is
OZ2M's PI4Ino. The K=32 code and the Fano algorithm go back to Layland,
Lushbaugh, Fano and Phil Karn, KA9Q; the sync-bit observation and the
"re-encode and compare" gate follow the pi4dec design.
