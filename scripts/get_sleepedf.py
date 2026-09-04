#!/usr/bin/env python3
"""Download Sleep-EDF Expanded hypnograms (annotations only) + subject metadata."""
import httpx, re, os, sys
BASE = "https://physionet.org/files/sleep-edfx/1.0.0/"
OUT = "datasets/sleep_edfx"
os.makedirs(f"{OUT}/sleep-cassette", exist_ok=True)
os.makedirs(f"{OUT}/sleep-telemetry", exist_ok=True)
with httpx.Client(timeout=120, follow_redirects=True) as c:
    for f in ("SC-subjects.xls", "ST-subjects.xls", "RECORDS", "SHA256SUMS.txt"):
        r = c.get(BASE + f)
        if r.status_code == 200: open(f"{OUT}/{f}", "wb").write(r.content); print("meta", f, len(r.content))
    for sub in ("sleep-cassette", "sleep-telemetry"):
        idx = c.get(BASE + sub + "/").text
        files = re.findall(r'href="([^"]*Hypnogram\.edf)"', idx)
        print(f"{sub}: {len(files)} hypnograms")
        for i, fn in enumerate(files, 1):
            p = f"{OUT}/{sub}/{fn}"
            if os.path.exists(p) and os.path.getsize(p) > 100: continue
            r = c.get(f"{BASE}{sub}/{fn}")
            if r.status_code == 200:
                open(p, "wb").write(r.content)
            if i % 25 == 0: print(f"  {i}/{len(files)}")
print("done")
