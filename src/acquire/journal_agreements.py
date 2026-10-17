"""Deposit agreement of each journal in the PMC journal list, linked to OpenAlex journals by ISSN-L.

Agreement to deposit is "All articles", "NIH funded articles", or none (no active PMC agreement). The list
is a current snapshot (16 October 2026). Writes data/journal_agreements.csv.

usage: python src/acquire/journal_agreements.py <journal_issn.csv with journal_id,issn_l>
"""

import sys

import pandas as pd

jlist = pd.read_csv("https://cdn.ncbi.nlm.nih.gov/pmc/home/jlist.csv", dtype=str)
jlist = jlist[jlist["Agreement Status"].eq("Active")]
by_issn = {}
for _, row in jlist.iterrows():
    for col in ("ISSN (print)", "ISSN (online)"):
        if isinstance(row[col], str):
            by_issn[row[col]] = row["Agreement to Deposit"]
journals = pd.read_csv(sys.argv[1]).dropna().drop_duplicates("journal_id")
journals["deposit_agreement"] = journals["issn_l"].map(by_issn).fillna("none")
journals.sort_values("journal_id").to_csv("data/journal_agreements.csv", index=False)
