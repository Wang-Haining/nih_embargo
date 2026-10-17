"""Publisher of each matched article from its DOI prefix, and whether it was a Method D publisher.

The publisher is the Crossref member that registered the DOI prefix, which reflects the publisher of the
article itself rather than a journal-level record. Method D publishers volunteered to deposit author
manuscripts to the NIH Manuscript Submission system on the author's behalf; the author then completes the
submission. The list is no longer on the NIH website; this is the list as archived on 15 June 2025:
https://web.archive.org/web/20250615180042/https://sharing.nih.gov/public-access-policy/submitting-to-pubmed-central/method-d-publishers
Writes data/doi_publishers.csv (one row per DOI prefix).

usage: python src/acquire/doi_publishers.py
"""

import time

import pandas as pd
import requests

# Crossref member names of the Method D publishers that appear among the matched articles. The archived list
# prints "Karger Publishing Lippincott, Williams And Wilkins" on one line; both are treated as listed, and
# Lippincott Williams & Wilkins registers its DOIs as Ovid Technologies (Wolters Kluwer Health).
METHOD_D = {
    "Elsevier BV",
    "Wiley",  # John Wiley And Sons Publishing
    "SAGE Publications",
    "Informa UK Limited",  # Taylor And Francis Group
    "JMIR Publications Inc.",
    "BMJ",
    "American Association for Cancer Research (AACR)",
    "American Association for the Advancement of Science (AAAS)",
    "American Society of Neuroradiology (ASNR)",
    "Royal Society of Chemistry (RSC)",
    "IOP Publishing",  # Institute Of Physics
    "MDPI AG",  # Multidisciplinary Digital Publishing Institute
    "S. Karger AG",
    "Ovid Technologies (Wolters Kluwer Health)",
    "The Journal of Rheumatology",
}

dois = pd.read_csv("data/matched_articles.csv", usecols=["doi"], dtype=str)["doi"].dropna()
prefixes = sorted(dois.str.replace("https://doi.org/", "", regex=False).str.split("/").str[0].unique())
rows = []
for prefix in prefixes:
    r = requests.get(f"https://api.crossref.org/prefixes/{prefix}", timeout=30)
    if r.status_code == 404:  # prefix registered outside Crossref
        rows.append({"prefix": prefix, "publisher": None})
        continue
    r.raise_for_status()
    rows.append({"prefix": prefix, "publisher": r.json()["message"]["name"]})
    time.sleep(0.2)
out = pd.DataFrame(rows)
out["method_d"] = out["publisher"].isin(METHOD_D)
out.to_csv("data/doi_publishers.csv", index=False)
print(f"{len(out)} prefixes, {out['method_d'].sum()} from Method D publishers")
