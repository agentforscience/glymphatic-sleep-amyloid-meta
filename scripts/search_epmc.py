#!/usr/bin/env python3
"""Search Europe PMC for papers, dedupe, and save results."""
import httpx, json, sys, time, os

QUERIES = [
 ('swsl_abeta', 'ABSTRACT:"slow wave sleep" AND ABSTRACT:"amyloid"'),
 ('sleep_abeta_pet', '(ABSTRACT:"amyloid" AND ABSTRACT:"PET") AND ABSTRACT:"sleep"'),
 ('glymphatic_human', 'ABSTRACT:"glymphatic" AND (ABSTRACT:"human" OR ABSTRACT:"MRI" OR ABSTRACT:"sleep")'),
 ('sleep_ad_meta', 'ABSTRACT:"sleep" AND ABSTRACT:"Alzheimer" AND (TITLE:"meta-analysis" OR TITLE:"systematic review")'),
 ('psg_amyloid', '(ABSTRACT:"polysomnography" OR ABSTRACT:"actigraphy") AND ABSTRACT:"amyloid"'),
 ('sws_memory_abeta', 'ABSTRACT:"slow oscillation" AND ABSTRACT:"amyloid"'),
 ('dtialps', 'ABSTRACT:"DTI-ALPS" OR ABSTRACT:"along the perivascular space"'),
 ('sleep_tau', 'ABSTRACT:"sleep" AND ABSTRACT:"tau PET"'),
 ('osa_amyloid', 'ABSTRACT:"sleep apnea" AND ABSTRACT:"amyloid"'),
 ('suvr_harmon', 'ABSTRACT:"centiloid" AND ABSTRACT:"amyloid PET"'),
 ('csf_abeta_sleep', 'ABSTRACT:"sleep deprivation" AND ABSTRACT:"amyloid-beta"'),
 ('meta_method', '(TITLE:"meta-analysis" AND ABSTRACT:"random-effects") AND ABSTRACT:"neuroimaging"'),
]

FIELDS = "id,source,pmid,pmcid,doi,title,authorString,journalTitle,pubYear,abstractText,isOpenAccess,citedByCount,fullTextUrlList"

def search(q, pagesize=40):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": q, "format": "json", "pageSize": pagesize,
              "resultType": "core", "sort": "CITED desc"}
    r = httpx.get(url, params=params, timeout=90)
    r.raise_for_status()
    return r.json().get("resultList", {}).get("result", [])

out = {}
for tag, q in QUERIES:
    try:
        res = search(q)
        print(f"[{tag}] {len(res)} hits", file=sys.stderr)
        for d in res:
            key = d.get("doi") or d.get("pmid") or d.get("id")
            if key and key not in out:
                d["_query_tag"] = tag
                out[key] = d
    except Exception as e:
        print(f"[{tag}] FAILED {e}", file=sys.stderr)
    time.sleep(1)

os.makedirs("paper_search_results", exist_ok=True)
with open("paper_search_results/epmc_raw.json","w") as f:
    json.dump(list(out.values()), f, indent=1)
print(f"TOTAL UNIQUE: {len(out)}", file=sys.stderr)
