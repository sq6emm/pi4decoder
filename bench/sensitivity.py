#!/usr/bin/env python3
"""Synthetic sensitivity sweep: decode probability vs S/N (2500 Hz) with
random frequency offset, start time and noise seed per trial.

Runs inside the pi4decoder-dev container (needs ./pi4gen and ./pi4decoder).
  python3 bench/sensitivity.py [--trials N] [--snr -18:-28] [--jobs 4] [--spread HZ]
                               [--drift HZ] [--text MSG] [--args "extra decoder args"]
"""
import argparse, json, os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser()
ap.add_argument("--trials", type=int, default=50)
ap.add_argument("--snr", default="-18:-28")
ap.add_argument("--step", type=float, default=1.0)
ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--spread", type=float, default=0.0)
ap.add_argument("--drift", type=float, default=0.0)
ap.add_argument("--text", default="OZ7IGY")
ap.add_argument("--args", default="--variants PI4")
ap.add_argument("--gen", default="./pi4gen")
ap.add_argument("--dec", default="./pi4decoder")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--json", default=None, help="write per-trial records here")
a = ap.parse_args()
lo, hi = [float(x) for x in a.snr.split(":")]
snrs = []
s = lo
while s >= hi - 1e-9:
    snrs.append(round(s, 2)); s -= a.step
rng = random.Random(a.seed)
tmp = tempfile.mkdtemp(prefix="pi4sens_")
jobs = []
for snr in snrs:
    for k in range(a.trials):
        jobs.append((snr, k, rng.uniform(-20, 20), rng.uniform(2.5, 4.0), rng.randrange(1 << 30)))

def run(job):
    snr, k, df, lead, seed = job
    wav = os.path.join(tmp, f"t_{snr}_{k}.wav")
    gen = [a.gen, "--text", a.text, "--snr", str(snr), "--f0", str(682.8125 + df), "--lead", f"{lead:.3f}",
           "--dur", "40", "--seed", str(seed)]
    if a.spread: gen += ["--spread", str(a.spread)]
    if a.drift: gen += ["--drift", str(a.drift)]
    subprocess.run(gen + [wav], check=True)
    out = subprocess.run([a.dec] + a.args.split() + [wav, "3.0"], capture_output=True, text=True).stdout
    os.unlink(wav)
    decs = [json.loads(l) for l in out.splitlines() if l.startswith("{")]
    ok = any(d["text"] == a.text for d in decs)
    wrong = [d["text"] for d in decs if d["text"] != a.text]
    return {"snr": snr, "k": k, "df": df, "lead": lead, "ok": ok, "wrong": wrong,
            "snr_est": next((d["snr_db"] for d in decs if d["text"] == a.text), None)}

with ThreadPoolExecutor(max_workers=a.jobs) as ex:
    recs = list(ex.map(run, jobs))
os.rmdir(tmp)
if a.json:
    with open(a.json, "w") as f: json.dump(recs, f)
print(f"text={a.text} spread={a.spread} drift={a.drift} args='{a.args}' trials/step={a.trials}")
print("  S/N dB   decoded   wrong   mean est")
for snr in snrs:
    r = [x for x in recs if x["snr"] == snr]
    ok = sum(x["ok"] for x in r); wrong = sum(len(x["wrong"]) for x in r)
    est = [x["snr_est"] for x in r if x["snr_est"] is not None]
    print(f"  {snr:6.1f}   {ok:3d}/{len(r):<3d}   {wrong:3d}    {sum(est)/len(est) if est else float('nan'):6.1f}")
