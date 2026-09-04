#!/usr/bin/env python3
"""Pull sentences containing effect-size statistics near sleep/amyloid terms."""
import glob, os, re, sys
SLEEP = r"(slow[- ]wave|SWS|SWA|N3|NREM|non-?rapid|sleep efficiency|total sleep time|TST|sleep duration|delta power|slow oscillation)"
AMY   = r"(amyloid|A\W?\W?42|A\W?\W?40|PiB|SUVR|florbetapir|centiloid|CSF A)"
STAT  = r"([rR]\s?=\s?[-−–]?[0-9.]+|β\s?=\s?[-−–]?[0-9.]+|[Bb]\s?=\s?[-−–]?[0-9.]+|OR\s?=|95%\s?CI|[Pp]\s?[<=]\s?0?\.[0-9]+|partial r|rho)"
pat = re.compile(f"(?=.*{SLEEP})(?=.*{AMY}).*", re.I)
files = sys.argv[1:] or sorted(glob.glob("papers/fulltext/*.txt"))
for f in files:
    txt = open(f, errors="ignore").read()
    txt = re.sub(r"\s*\n\s*", " ", txt)
    sents = re.split(r"(?<=[.;])\s+", txt)
    hits = [s.strip() for s in sents if len(s) < 700 and pat.match(s) and re.search(STAT, s)]
    if hits:
        print("\n" + "="*100)
        print(os.path.basename(f))
        for h in dict.fromkeys(hits): print("  •", h)
