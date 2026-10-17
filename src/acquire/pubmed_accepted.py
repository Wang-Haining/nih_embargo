"""PubMed accepted dates for the NIH-supported author manuscripts in the matched sample.

PMC rarely records an accepted date for author manuscripts, so the PubMed history date with status
"accepted" is used. Writes data/pubmed_accepted.csv. Requires network access (NCBI E-utilities).

usage: python src/acquire/pubmed_accepted.py
"""

import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

rec = pd.read_csv("data/pmc_records.csv")
art = pd.read_csv("data/matched_articles.csv", dtype={"pmid": str}).merge(rec, on="pmcid")
pmids = art.loc[art["nih"] & art["author_manuscript"].eq("yes"), "pmid"].unique().tolist()
rows = []
for i in range(0, len(pmids), 200):
    form = urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids[i : i + 200]), "retmode": "xml"}
    ).encode()
    req = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=form)
    root = ET.fromstring(urllib.request.urlopen(req, timeout=120).read())
    for article in root.iter("PubmedArticle"):
        accepted = None
        for d in article.iter("PubMedPubDate"):
            if d.get("PubStatus") == "accepted":
                accepted = f"{d.findtext('Year')}-{int(d.findtext('Month')):02d}-{int(d.findtext('Day')):02d}"
        rows.append({"pmid": article.findtext(".//MedlineCitation/PMID"), "pubmed_accepted_date": accepted})
    time.sleep(0.6)
pd.DataFrame(rows).to_csv("data/pubmed_accepted.csv", index=False)
