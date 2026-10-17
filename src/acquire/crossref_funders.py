#!/usr/bin/env python3
"""Works carrying the NIH funder identifier in Crossref (10.13039/100000002), January to April 2024 to 2026.

Writes data/crossref_nih_dois.csv and data/crossref_date_audit.csv. The year check uses
data/crossref_openalex_links.csv: each DOI linked by DOI to its record in the OpenAlex snapshot of
26 June 2026 (publication date and year, PMID, journal, DOAJ status, NIH among the OpenAlex funders).
"""

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
URL = "https://api.crossref.org/funders/10.13039/100000002/works"


def fetch_year(year):
    rows, cursor = [], "*"
    while True:
        query = urllib.parse.urlencode(
            {
                "filter": f"from-pub-date:{year}-01-01,until-pub-date:{year}-04-30,type:journal-article",
                "rows": 1000,
                "cursor": cursor,
                "select": "DOI,published,type",
            }
        )
        request = urllib.request.Request(
            URL + "?" + query, headers={"User-Agent": "nih-embargo/1.0 (hw56@iu.edu)"}
        )
        with urllib.request.urlopen(request, timeout=90) as response:
            message = json.load(response)["message"]
        items = message["items"]
        rows.extend({"doi": x["DOI"].lower(), "crossref_year": year} for x in items)
        if len(items) < 1000:
            break
        cursor = message["next-cursor"]
        time.sleep(0.1)
    print(f"Crossref {year}: {len(rows)}", flush=True)
    return rows


def main():
    raw_path = ROOT / "work/crossref_nih_dois.parquet"
    if raw_path.exists():
        works = pd.read_parquet(raw_path)
    else:
        works = pd.DataFrame(sum((fetch_year(y) for y in (2024, 2025, 2026)), []))
        if works.empty or works.doi.duplicated().any():
            raise ValueError(
                f"Crossref NIH cohort invalid: rows={len(works)}, "
                f"duplicate_doi={works.doi.duplicated().sum()}"
            )
        works.to_parquet(raw_path, index=False)
    works[["doi"]].to_csv(ROOT.parents[1] / "data/crossref_nih_dois.csv", index=False)
    dates = pd.read_csv(ROOT.parents[1] / "data/crossref_openalex_links.csv", dtype={"pmid": "string"})
    dates["linked"], dates["same_year"] = dates.oa_year.notna(), dates.oa_year.eq(dates.crossref_year)
    dates["eligible"] = dates.same_year & dates.oa_nih.eq(True) & dates.pmid.notna() & dates.pmid.ne("")
    date_summary = (
        dates.groupby("crossref_year")
        .agg(
            n=("doi", "size"),
            linked=("linked", "sum"),
            same_year=("same_year", "sum"),
            eligible=("eligible", "sum"),
        )
        .reset_index()
    )
    date_path = ROOT.parents[1] / "data/crossref_date_audit.csv"
    date_summary.to_csv(date_path, index=False)
    print(f"Crossref NIH funder DOIs: n={len(works)}")
    print(date_summary.to_string(index=False))


if __name__ == "__main__":
    main()
