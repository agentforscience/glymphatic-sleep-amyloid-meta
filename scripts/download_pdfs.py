#!/usr/bin/env python3
"""Download OA PDFs for selected papers: Unpaywall (all locations) -> Europe PMC -> PMC."""
import json, os, re, sys, time, httpx

os.makedirs("papers", exist_ok=True)
papers = json.load(open("paper_search_results/selected.json"))
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36",
      "Accept": "application/pdf,text/html,*/*"}

def slug(x):
    au = (x.get("authorString") or "unknown").split(",")[0].split()[0]
    au = re.sub(r"[^A-Za-z]", "", au) or "unknown"
    t = re.sub(r"[^a-z0-9]+", "_", (x.get("title") or "")[:55].lower()).strip("_")
    return f"{au}{x.get('pubYear','')}_{t}"

def urls_for(c, x):
    out = []
    doi = x.get("doi")
    if doi:
        try:
            j = c.get(f"https://api.unpaywall.org/v2/{doi}?email=chicagohailab@gmail.com").json()
            for L in (j.get("oa_locations") or []):
                for u in (L.get("url_for_pdf"), L.get("url")):
                    if u: out.append(u)
        except Exception: pass
    pmcid = x.get("pmcid")
    if pmcid:
        out.append(f"https://europepmc.org/articles/{pmcid}?pdf=render")
        out.append(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextPdf")
    for ft in (x.get("fullTextUrlList") or {}).get("fullTextUrl", []):
        if ft.get("availabilityCode") == "OA": out.append(ft["url"])
    # normalise old PMC repository urls into pdf endpoints
    extra = []
    for u in out:
        m = re.search(r"pmc/articles/(?:PMC)?(\d+)", u)
        if m: extra.append(f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{m.group(1)}/pdf/")
    seen, final = set(), []
    for u in out + extra:
        if u not in seen: seen.add(u); final.append(u)
    return final

ok, fail = [], []
with httpx.Client(timeout=90, follow_redirects=True, headers=UA) as c:
    for i, x in enumerate(papers, 1):
        fn = f"papers/{slug(x)}.pdf"
        if os.path.exists(fn) and os.path.getsize(fn) > 20000:
            ok.append((fn, x)); continue
        got = None
        for u in urls_for(c, x):
            try:
                r = c.get(u)
                if r.status_code == 200 and r.content[:5] == b"%PDF-" and len(r.content) > 20000:
                    open(fn, "wb").write(r.content); got = u; break
            except Exception: pass
            time.sleep(0.25)
        if got:
            x["_pdf"] = fn; x["_pdf_src"] = got; ok.append((fn, x))
            print(f"[{i:>2}/{len(papers)}] OK   {os.path.basename(fn)[:70]}")
        else:
            fail.append(x)
            print(f"[{i:>2}/{len(papers)}] MISS {(x.get('title') or '')[:70]}")
        time.sleep(0.3)

print(f"\nDOWNLOADED {len(ok)} / {len(papers)}")
json.dump([x for _, x in ok], open("paper_search_results/downloaded.json", "w"), indent=1)
json.dump(fail, open("paper_search_results/failed.json", "w"), indent=1)
