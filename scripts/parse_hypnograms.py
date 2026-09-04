#!/usr/bin/env python3
"""Parse Sleep-EDF Expanded EDF+ hypnogram files into per-record sleep-stage summaries.

Two windows are computed for every recording:
  * full   -- the entire annotated recording (Sleep-Cassette spans ~20 h of ambulatory
              recording, so this includes daytime wake and naps)
  * noct   -- the nocturnal sleep period only, anchored at the lights-off clock time
              from SC-subjects.xls and closed at the last sleep epoch within 12 h.
The nocturnal window is the one comparable to in-lab PSG used in the amyloid literature.
"""
import os, re, glob, csv, sys, datetime as dt
import pandas as pd

STAGE_MAP = {
    "Sleep stage W": "W", "Sleep stage 1": "N1", "Sleep stage 2": "N2",
    "Sleep stage 3": "N3", "Sleep stage 4": "N3",   # R&K S3+S4 -> AASM N3
    "Sleep stage R": "REM", "Sleep stage ?": "UNK", "Movement time": "MT",
}
SLEEP = ("N1", "N2", "N3", "REM")

def read_edf(path):
    """Return (start_datetime, [(onset_s, dur_s, label), ...])."""
    with open(path, "rb") as f:
        hdr = f.read(256)
        sd, st_ = hdr[168:176].decode(), hdr[176:184].decode()
        n_records = int(hdr[236:244].decode().strip())
        n_sig = int(hdr[252:256].decode().strip())
        sig = f.read(256 * n_sig)
        off = n_sig * (16 + 80 + 8 + 8 + 8 + 8 + 8 + 80)
        nsamp = [int(sig[off + 8 * i: off + 8 * (i + 1)].decode().strip()) for i in range(n_sig)]
        data = f.read()
    dd, mm, yy = (int(v) for v in sd.split("."))
    yy += 1900 if yy >= 85 else 2000
    hh, mi, ss = (int(v) for v in st_.replace(":", ".").split("."))
    start = dt.datetime(yy, mm, dd, hh, mi, ss)
    bpr = sum(n * 2 for n in nsamp)
    if n_records < 0: n_records = len(data) // bpr
    anns, pos = [], 0
    for _ in range(n_records):
        chunk = data[pos: pos + bpr]; pos += bpr
        for tal in chunk.split(b"\x00"):
            if not tal.strip(b"\x00 "): continue
            parts = tal.split(b"\x14")
            if len(parts) < 2: continue
            tm = parts[0].decode("ascii", "ignore").split("\x15")
            try:
                onset = float(tm[0]); dur = float(tm[1]) if len(tm) > 1 and tm[1] else 0.0
            except ValueError:
                continue
            for lab in parts[1:]:
                lab = lab.decode("ascii", "ignore").strip("\x00 ")
                if lab: anns.append((onset, dur, lab))
    return start, anns

def summarise(epochs, prefix):
    """epochs: list of (onset_s, dur_s, stage). Returns dict of metrics."""
    mins = {}
    for _, d, s in epochs: mins[s] = mins.get(s, 0.0) + d / 60.0
    tst = sum(mins.get(s, 0.0) for s in SLEEP)
    idx = [i for i, (_, _, s) in enumerate(epochs) if s in SLEEP]
    out = {}
    if not idx or tst <= 0:
        return {f"{prefix}_{k}": None for k in
                ["tst_min","spt_min","sleep_eff_pct","waso_min","n1_pct_tst","n2_pct_tst",
                 "n3_pct_tst","rem_pct_tst","n3_min"]}
    lo, hi = idx[0], idx[-1]
    spt = (epochs[hi][0] + epochs[hi][1] - epochs[lo][0]) / 60.0
    waso = sum(d for o, d, s in epochs[lo:hi + 1] if s == "W") / 60.0
    out[f"{prefix}_tst_min"] = round(tst, 2)
    out[f"{prefix}_spt_min"] = round(spt, 2)
    out[f"{prefix}_sleep_eff_pct"] = round(100 * tst / spt, 3) if spt > 0 else None
    out[f"{prefix}_waso_min"] = round(waso, 2)
    out[f"{prefix}_n3_min"] = round(mins.get("N3", 0.0), 2)
    for s in SLEEP:
        out[f"{prefix}_{s.lower()}_pct_tst"] = round(100 * mins.get(s, 0.0) / tst, 3)
    return out

# lights-off lookup for Sleep-Cassette
sc_meta = pd.read_excel("datasets/sleep_edfx/SC-subjects.xls")
lights = {(int(r.subject), int(r.night)): r.LightsOff for r in sc_meta.itertuples()}

rows = []
for path in sorted(glob.glob("datasets/sleep_edfx/*/*-Hypnogram.edf")):
    rec = os.path.basename(path).replace("-Hypnogram.edf", "")
    cohort = "cassette" if rec.startswith("SC") else "telemetry"
    try:
        start, anns = read_edf(path)
    except Exception as e:
        print(f"ERR {rec}: {e}", file=sys.stderr); continue
    ep = [(o, d, STAGE_MAP[l]) for o, d, l in anns if l in STAGE_MAP and STAGE_MAP[l] in SLEEP + ("W",)]
    if not ep: continue
    row = {"record": rec, "cohort": cohort, "start_datetime": start.isoformat()}
    row.update(summarise(ep, "full"))
    # nocturnal window
    if cohort == "cassette":
        subj, night = int(rec[3:5]), int(rec[5])
        lo_t = lights.get((subj, night))
        row["subj_num"], row["night_n"] = subj, night
        if lo_t is not None:
            lo_dt = dt.datetime.combine(start.date(), lo_t if isinstance(lo_t, dt.time)
                                        else dt.datetime.strptime(str(lo_t), "%H:%M:%S").time())
            if lo_dt < start: lo_dt += dt.timedelta(days=1)
            t0 = (lo_dt - start).total_seconds()
            win = [(o, d, s) for o, d, s in ep if o + d > t0 and o < t0 + 12 * 3600]
            row["lights_off"] = str(lo_t)
            row.update(summarise(win, "noct"))
    rows.append(row)

cols = sorted({k for r in rows for k in r})
lead = ["record","cohort","subj_num","night_n","start_datetime","lights_off"]
cols = lead + [c for c in cols if c not in lead]
out = "datasets/sleep_edfx/hypnogram_stages.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(f"wrote {out}: {len(rows)} recordings")
