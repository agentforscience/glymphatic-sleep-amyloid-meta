#!/usr/bin/env python3
"""Retry failed downloads by scraping the PMC article page for its real PDF link."""
import json, os, re, time, httpx
from urllib.parse import urljoin

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"}
fail = json.load(open("paper_search_results/failed.json"))
done = json.load(open("paper_search_results/downloaded.json"))

def slug(x):
    au = (x.get("authorString") or "unknown").split(",")[0].split()[0]
    au = re.sub(r"[^A-Za-z]", "", au) or "unknown"
    t = re.sub(r"[^a-z0-9]+", "_", (x.get("title") or "")[:55].lower()).strip("_")
    return f"{au}{x.get('pubYear','')}_{t}"

still, newly = [], []
with httpx.Client(timeout=90, follow_redirects=True, headers=UA) as c:
    for x in fail:
        fn = f"papers/{slug(x)}.pdf"
        if os.path.exists(fn) and os.path.getsize(fn) > 20000:
            newly.append(x); continue
        pmcid = x.get("pmcid")
        got = None
        pages = []
        if pmcid: pages.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
        if x.get("pmid"): pages.append(f"https://pmc.ncbi.nlm.nih.gov/articles/pmid/{x['pmid']}/")
        for page in pages:
            try:
                art = c.get(page)
                if art.status_code != 200: continue
                base = str(art.url)
                links = re.findall(r'href="([^"]*\.pdf[^"]*)"', art.text)
                links = [L for L in links if "supplement" not in L.lower() and "/bin/" not in L]
                for L in links:
                    u = urljoin(base, L)
                    r = c.get(u, headers={**UA, "Referer": base})
                    if r.status_code == 200 and r.content[:5] == b"%PDF-" and len(r.content) > 20000:
                        open(fn, "wb").write(r.content); got = u; break
                    time.sleep(0.4)
            except Exception: pass
            if got: break
            time.sleep(0.5)
        if got:
            x["_pdf"] = fn; x["_pdf_src"] = got; newly.append(x)
            print(f"OK   {os.path.basename(fn)[:70]}")
        else:
            still.append(x); print(f"MISS {(x.get('title') or '')[:70]}")
        time.sleep(0.5)

print(f"\nrecovered {len(newly)} ; still missing {len(still)}")
json.dump(done + newly, open("paper_search_results/downloaded.json", "w"), indent=1)
json.dump(still, open("paper_search_results/failed.json", "w"), indent=1)
