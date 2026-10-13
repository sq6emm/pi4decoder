#!/usr/bin/env python3
"""Compile benchmark tables from the raw outputs in bench/results.

Reads pi4decoder_{syn,oz2m,rec}.txt (and the reference decoders' files when
present) and prints Markdown tables: files with a correct decode per set,
false/garbage decodes, CPU time and peak memory. Ground truth: synthetic files
carry OZ7IGY (noise_*.wav carry nothing); OZ2M files carry the beacon named
in GROUND_TRUTH; own minutes carry OZ7IGY only when named *_inj16.wav.
"""
import json, os, re, sys
from collections import defaultdict

R = os.path.join(os.path.dirname(__file__), "results")
GROUND_TRUTH = {"10G": "ED5YAE/B", "28M": "G3ZJO", "432M": "OZ7IGY", "70M": "OZ7IGY", "70M_1308": "IW9GDC/B"}

def truth_oz2m(name):
    if name.startswith("70M_1308"): return "IW9GDC/B"
    if name.startswith("pi4_cwid_car"): return "OZ7IGY"
    return GROUND_TRUTH.get(name.split("_")[0])

def parse(path):
    """-> {file: {"decodes": [dict], "user": s, "rss": kB}}"""
    out = {}
    if not os.path.exists(path): return out
    cur = None
    for ln in open(path, errors="replace"):
        ln = ln.rstrip("\n")
        if ln.startswith("== ") and ln.endswith(".wav"):
            cur = ln[3:].strip(); out[cur] = {"decodes": [], "user": None, "rss": None}
        elif cur is None: continue
        elif ln.startswith("{"):
            try: out[cur]["decodes"].append(json.loads(ln))
            except ValueError: pass
        elif ln.startswith("RUSAGE"):
            m = re.search(r"user ([\d.]+) s.*maxrss (\d+) kB", ln)
            if m: out[cur]["user"] = float(m.group(1)); out[cur]["rss"] = int(m.group(2))
        else:
            # reference decoders' plain-text lines: pi4d "SNR DT Hz Variant Message", MSHV tab separated
            parts = ln.split("\t")
            if len(parts) >= 8 and parts[7].strip() and parts[7].strip() != "?":
                out[cur]["decodes"].append({"text": parts[7].strip(), "ref": True})
            else:
                m = re.match(r"\s*-?\d+\s+-?[\d.]+\s+\d+\s+(PI4\S*)\s+(\S+)", ln)
                if m: out[cur]["decodes"].append({"text": m.group(2), "variant": m.group(1), "ref": True})
    return out

def summarize(name, data, truth):
    rows = defaultdict(lambda: [0, 0, 0])   # group -> [files, correct, garbage]
    times, rss = [], []
    for f, d in sorted(data.items()):
        want = truth(f)
        texts = [x["text"] for x in d["decodes"]]
        g = group_of(f)
        rows[g][0] += 1
        if want and want in texts: rows[g][1] += 1
        rows[g][2] += sum(1 for t in texts if t != want)
        if d["user"] is not None: times.append(d["user"]); rss.append(d["rss"])
    return rows, times, rss

def group_of(f):
    if f.startswith("s-10"): return "Synthetic -10 dB"
    if f.startswith("s-16"): return "Synthetic -16 dB"
    if f.startswith("syn_-20") or f.startswith("syn_-22"): return "Synthetic -20/-22 dB"
    if f.startswith("syn_"): return "Synthetic -24..-32 dB"
    if f.startswith("spread_"): return "Synthetic spread 10-100 Hz, -16 dB"
    if f.startswith("drift_"): return "Synthetic drift 2 Hz, -16 dB"
    if f.startswith("noise_"): return "Synthetic pure noise"
    if f.startswith("pi4_cwid_car"): return "OZ2M pi4_cwid_car (mp3/mp4 audio)"
    if f.startswith("10G"): return "OZ2M 10 GHz ED5YAE/B (spread)"
    if f.startswith("28M"): return "OZ2M 28 MHz G3ZJO"
    if f.startswith("432M"): return "OZ2M 432 MHz OZ7IGY 2013"
    if f.startswith("70M_1308"): return "OZ2M 70 MHz IW9GDC/B"
    if f.startswith("70M_1409"): return "OZ2M 70 MHz OZ7IGY 2014"
    if f.startswith("70M_1503"): return "OZ2M 70 MHz OZ7IGY 2015 (weak)"
    if f.endswith("_inj16.wav"): return "Own minutes + -16 dB OZ7IGY"
    if re.match(r"\d{8}_\d{4}\.wav", f): return "Own minutes (no signal)"
    return "other"

def truth_syn(f): return None if f.startswith("noise_") else "OZ7IGY"
def truth_rec(f): return "OZ7IGY" if f.endswith("_inj16.wav") else None

sets = [("syn", truth_syn), ("oz2m", truth_oz2m), ("rec", truth_rec)]
decoders = ["pi4decoder", "pi4dec", "pi4d", "mshv"]
table = defaultdict(dict); perf = {}
for dec in decoders:
    for s, truth in sets:
        data = parse(os.path.join(R, f"{dec}_{s}.txt"))
        if not data: continue
        rows, times, rss = summarize(dec, data, truth)
        for g, (n, ok, bad) in rows.items():
            table[g][dec] = (n, ok, bad)
        if times:
            perf.setdefault(dec, {})[s] = (sum(times) / len(times), max(rss))
present = [d for d in decoders if any(d in v for v in table.values())]
print("| Test | Files | " + " | ".join(present) + " |")
print("|---|---|" + "---|" * len(present))
for g in sorted(table, key=lambda g: (not g.startswith("Synthetic"), g)):
    n = max(v[0] for v in table[g].values())
    cells = []
    for d in present:
        v = table[g].get(d)
        cells.append("-" if v is None else (f"{v[1]}" + (f" (+{v[2]} wrong)" if v[2] else "")))
    print(f"| {g} | {n} | " + " | ".join(cells) + " |")
print()
print("| Decoder | set | mean CPU s/file | peak RSS kB |")
print("|---|---|---|---|")
for d in present:
    for s, (t, m) in perf.get(d, {}).items():
        print(f"| {d} | {s} | {t:.2f} | {m} |")
noise = os.path.join(R, "pi4decoder_noise.txt")
if os.path.exists(noise):
    for ln in open(noise):
        if ln.startswith("test_noise:"): print("\n" + ln.strip())
