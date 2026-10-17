"""Check the OAI-PMH go-live date against PMC's Live Date search field ([pmcrdat]) for a seeded sample.

PMC defines its Live Date as "the date the article was first available in PMC"; for embargoed articles it is
updated "to the date the full text became available" (PMC User Guide, https://pmc.ncbi.nlm.nih.gov/about/userguide/).
For each sampled record, two E-utilities searches ask whether PMC's Live Date equals the OAI-PMH go-live date
and whether it equals the release date. Sample: 50 final published versions whose release date precedes going
live by more than 30 days, 25 final versions whose two dates agree, and 25 author manuscripts. Writes
data/pmc_live_check.csv. About 200 requests; run outside PMC peak hours.

usage: python src/acquire/pmc_live_check.py
"""

import time

import pandas as pd
import requests

SEED = 20261017
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

rec = pd.read_csv("data/pmc_records.csv", parse_dates=["pmc_release_date", "pmc_live_date"])
rec = rec[rec["record_available"] & rec["pmc_live_date"].notna() & rec["pmc_release_date"].notna()]
final = rec["author_manuscript"].ne("yes")
gap = (rec["pmc_live_date"] - rec["pmc_release_date"]).dt.days
LATE = "final, release more than 30 days before going live"
strata = {
    LATE: rec[final & (gap > 30)].sample(50, random_state=SEED),
    "final, dates agree": rec[final & (gap == 0)].sample(25, random_state=SEED),
    "author manuscript": rec[~final].sample(25, random_state=SEED),
}


def live_date_is(uid: str, date: pd.Timestamp) -> bool:
    term = f"{uid}[uid] AND {date:%Y/%m/%d}[pmcrdat]"
    params = {"db": "pmc", "term": term, "retmode": "json", "tool": "nih_embargo"}
    r = requests.get(EUTILS, params=params, timeout=30)
    r.raise_for_status()
    time.sleep(0.4)
    return int(r.json()["esearchresult"]["count"]) == 1


rows = []
for stratum, sample in strata.items():
    for _, x in sample.iterrows():
        uid = x["pmcid"][3:]
        rows.append(
            {
                "pmcid": x["pmcid"],
                "stratum": stratum,
                "release_date": x["pmc_release_date"].date(),
                "live_date": x["pmc_live_date"].date(),
                "pmc_live_date_equals_live": live_date_is(uid, x["pmc_live_date"]),
                "pmc_live_date_equals_release": live_date_is(uid, x["pmc_release_date"]),
            }
        )
out = pd.DataFrame(rows)
out.to_csv("data/pmc_live_check.csv", index=False)
print(out.groupby("stratum")[["pmc_live_date_equals_live", "pmc_live_date_equals_release"]].sum())
