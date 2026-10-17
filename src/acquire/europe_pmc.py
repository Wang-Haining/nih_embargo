"""Europe PMC full-text received date for a seeded random sample of 200 matched 2026 articles in PMC.

Used to check the PMC go-live date against an independent clock (SI Appendix, S4). Writes
data/europe_pmc_sample.csv. Requires network access.

usage: python src/acquire/europe_pmc.py
"""

import json
import time
import urllib.request

import pandas as pd

rec = pd.read_csv("data/pmc_records.csv")
art = pd.read_csv("data/matched_articles.csv", dtype={"pmid": str})
pool = art[(art["year"] == 2026) & art["pmcid"].notna()].merge(rec, on="pmcid")
pool = pool[pool["pmc_live_date"].notna()].sample(200, random_state=20261016)
rows = []
for pmcid in pool["pmcid"]:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=PMCID:{pmcid}&resultType=core&format=json"
    try:
        res = json.load(urllib.request.urlopen(url, timeout=60))["resultList"]["result"]
        rows.append(
            {"pmcid": pmcid, "fulltext_received_date": res[0].get("fullTextReceivedDate") if res else None}
        )
    except OSError:
        rows.append({"pmcid": pmcid, "fulltext_received_date": None})
    time.sleep(0.2)
pd.DataFrame(rows).to_csv("data/europe_pmc_sample.csv", index=False)
