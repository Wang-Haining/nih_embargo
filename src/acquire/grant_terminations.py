"""NIH awards terminated since 2025, from Grant Witness (https://grantwitness.org/nih).

Grant Witness compiles NIH grant disruptions from NIH RePORTER ("Terminated - departmental authority"), the
HHS list of terminated awards, USAspending, court records and researcher reports; methods at
https://grantwitness.org/nih/updates/methodology. Its award key is the institute code and six-digit serial
number (e.g. CA123456). An award counts as terminated if its detailed status is Terminated or Reinstated;
awards whose funds were only frozen, and supplement-level records, are excluded. Only the key, the date of the
earliest disruption, the latest restoration and the status are kept. Writes data/grant_terminations.csv.

usage: python src/acquire/grant_terminations.py
"""

import io

import pandas as pd
import requests

URL = "https://data.grantwitness.org/nih/dl-table.csv"

r = requests.get(URL, timeout=300)
r.raise_for_status()
print("Grant Witness data last modified:", r.headers.get("last-modified"))
gw = pd.read_csv(io.StringIO(r.text), dtype=str)
terminated = gw["Detailed Status"].isin(["Terminated", "Reinstated (confirmed)", "Reinstated (unconfirmed)"])
parent = gw["Subaward ID"].isna()
out = gw.loc[
    terminated & parent, ["Award ID", "Earliest Disruption", "Latest Restoration", "Detailed Status"]
]
out.columns = ["award_key", "earliest_disruption", "latest_restoration", "detailed_status"]
out = out.sort_values(["award_key", "earliest_disruption"]).drop_duplicates("award_key")
out.to_csv("data/grant_terminations.csv", index=False)
print(f"{len(gw)} disrupted records, {len(out)} terminated parent awards")
