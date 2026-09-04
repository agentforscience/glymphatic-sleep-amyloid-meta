#!/usr/bin/env python3
"""Recover remaining papers: S2 OA PDFs, then PMC full-text HTML -> text + rendered PDF."""
import json, os, re, time, httpx
from bs4 import BeautifulSoup
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from xml.sax.saxutils import escape

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"}
os.makedirs("papers/fulltext", exist_ok=True)
fail = json.load(open("paper_search_results/failed.json"))
done = json.load(open("paper_search_results/downloaded.json"))

def slug(x):
    au = (x.get("authorString") or "unknown").split(",")[0].split()[0]
    au = re.sub(r"[^A-Za-z]", "", au) or "unknown"
    t = re.sub(r"[^a-z0-9]+", "_", (x.get("title") or "")[:55].lower()).strip("_")
    return f"{au}{x.get('pubYear','')}_{t}"

S2_PDF = {
 "10.1002/ana.26604": "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/ana.26604",
 "10.1111/jon.12837": "http://minerva-access.unimelb.edu.au/bitstreams/334a0f27-2ad0-5efe-94e0-2be1d1c9a016/download",
}

def text_to_pdf(path, title, body):
    doc = SimpleDocTemplate(path, pagesize=letter, topMargin=48, bottomMargin=48,
                            leftMargin=54, rightMargin=54)
    ss = getSampleStyleSheet()
    flow = [Paragraph(escape(title), ss["Title"]), Spacer(1, 10)]
    for para in body.split("\n"):
        p = para.strip()
        if not p: continue
        style = ss["Heading2"] if (len(p) < 90 and p.isupper()) else ss["BodyText"]
        flow.append(Paragraph(escape(p), style))
        flow.append(Spacer(1, 4))
    doc.build(flow)

def pmc_fulltext(c, pmcid):
    r = c.get(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
    if r.status_code != 200: return None
    soup = BeautifulSoup(r.text, "lxml")
    for bad in soup(["script", "style", "nav", "header", "footer", "aside"]): bad.decompose()
    main = soup.find("section", {"aria-label": "Article content"}) or soup.find("main") or soup.find("article")
    if main is None: return None
    parts = []
    for el in main.find_all(["h1","h2","h3","p","li","td","caption"]):
        t = el.get_text(" ", strip=True)
        if t and len(t) > 2:
            parts.append(t.upper() if el.name in ("h1","h2","h3") and len(t) < 90 else t)
    txt = "\n".join(dict.fromkeys(parts))
    return txt if len(txt) > 6000 else None

still, newly = [], []
with httpx.Client(timeout=90, follow_redirects=True, headers=UA) as c:
    for x in fail:
        fn = f"papers/{slug(x)}.pdf"
        if os.path.exists(fn) and os.path.getsize(fn) > 20000:
            newly.append(x); continue
        # 1) direct S2 OA pdf
        u = S2_PDF.get(x.get("doi"))
        if u:
            try:
                r = c.get(u)
                if r.status_code == 200 and r.content[:5] == b"%PDF-" and len(r.content) > 20000:
                    open(fn, "wb").write(r.content)
                    x["_pdf"] = fn; x["_pdf_src"] = u; x["_mode"] = "publisher-pdf"
                    newly.append(x); print(f"PDF  {os.path.basename(fn)[:66]}"); continue
            except Exception: pass
        # 2) PMC full-text HTML
        if x.get("pmcid"):
            try:
                txt = pmc_fulltext(c, x["pmcid"])
            except Exception: txt = None
            if txt:
                tp = f"papers/fulltext/{slug(x)}.txt"
                open(tp, "w").write(txt)
                text_to_pdf(fn, x.get("title") or slug(x), txt)
                x["_pdf"] = fn; x["_text"] = tp; x["_mode"] = "pmc-html-rendered"
                newly.append(x); print(f"HTML {os.path.basename(fn)[:66]} ({len(txt)} chars)"); time.sleep(0.6); continue
        still.append(x); print(f"MISS {(x.get('title') or '')[:66]}")
        time.sleep(0.5)

print(f"\nrecovered {len(newly)} ; still abstract-only {len(still)}")
json.dump(done + newly, open("paper_search_results/downloaded.json", "w"), indent=1)
json.dump(still, open("paper_search_results/abstract_only.json", "w"), indent=1)
