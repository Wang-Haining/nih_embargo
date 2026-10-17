"""Grant numbers listed in PubMed for the NIH-supported articles of the 2024, 2025 and 2026 cohorts.

Covers the NIH-supported articles of the primary sample (2025, 2026) and of the pre-policy sample's 2024
cohort.
Writes data/pubmed_grants.csv with one row per listed grant (pmid, grant_id, agency).

usage: python src/acquire/pubmed_grants.py
"""

import time

import pandas as pd
import requests
from lxml import etree

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

prim = pd.read_csv("data/matched_articles.csv", dtype={"pmid": str})
plac = pd.read_csv("data/placebo_articles.csv", dtype={"pmid": str})
pre_policy = plac.loc[plac["nih"] & (plac["year"] == 2024), "pmid"]
pmids = sorted(set(prim.loc[prim["nih"], "pmid"]) | set(pre_policy))
rows = []
for i in range(0, len(pmids), 200):
    data = {"db": "pubmed", "id": ",".join(pmids[i : i + 200]), "retmode": "xml", "tool": "nih_embargo"}
    r = requests.post(EFETCH, data=data, timeout=120)
    r.raise_for_status()
    for article in etree.fromstring(r.content).iterfind(".//PubmedArticle"):
        pmid = article.findtext(".//MedlineCitation/PMID")
        for grant in article.iterfind(".//GrantList/Grant"):
            row = {"pmid": pmid, "grant_id": grant.findtext("GrantID"), "agency": grant.findtext("Agency")}
            rows.append(row)
    time.sleep(0.4)
out = pd.DataFrame(rows)
out.to_csv("data/pubmed_grants.csv", index=False)
print(f"{len(pmids)} articles queried, {out['pmid'].nunique()} with grants, {len(out)} grant rows")
