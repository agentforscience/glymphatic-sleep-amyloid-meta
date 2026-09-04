import re,sys,glob,os
STAT=re.compile(r"(β\s*=\s*[-−–]?[0-9.]+|\br\s*=\s*[-−–]?0?\.[0-9]+|\bB\s*=\s*[-−–]?[0-9.]+|\bOR\s*[=:]\s*[0-9.]+|95%\s*CI|\bd\s*=\s*[-−–]?[0-9.]+|\bt\s*=\s*[-−–]?[0-9.]+)")
SLEEP=re.compile(r"(slow[- ]?wave|SWS|\bSWA\b|\bN3\b|NREM|non-?REM|sleep efficien|total sleep time|\bTST\b|sleep duration|delta power|slow oscillat|sleep quality|PSQI|WASO)",re.I)
AMY=re.compile(r"(amyloid|A[βB]\s?4[02]|PiB|SUVR|centiloid|florbetap|flutemet|A\W?\W?42|\bAβ\b)",re.I)
for f in sys.argv[1:]:
    t=re.sub(r"\s+"," ",open(f,errors="ignore").read())
    out=[]
    for m in STAT.finditer(t):
        s=max(0,m.start()-260); e=min(len(t),m.end()+180)
        w=t[s:e]
        if SLEEP.search(w) and AMY.search(w): out.append(w)
    # dedupe overlapping
    ded=[]
    for w in out:
        if not any(w[:60] in d for d in ded): ded.append(w)
    if ded:
        print("\n"+"="*95); print(os.path.basename(f))
        for w in ded[:12]: print("  •",w.strip())
