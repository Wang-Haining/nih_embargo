#!/usr/bin/env python3
"""NIH-supported articles identified from funder metadata (Crossref and OpenAlex), independent of PubMed.

Keeps DOIs whose OpenAlex publication year matches the Crossref year and whose OpenAlex record also lists
NIH funding, in the journal-month cells of the matched sample. Writes data/funder_replication_articles.csv;
outcomes are joined from data/pmc_records.csv by src/estimate.py.
"""

from pathlib import Path

import pandas as pd
from build_cohort import efetch, pmc_map

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parents[1] / "data"


def main():
    audit = pd.read_csv(DATA / "crossref_openalex_links.csv", dtype={"pmid": "string"})
    matched = pd.read_csv(DATA / "matched_articles.csv", dtype={"pmid": str})
    denominator = pd.read_parquet(ROOT / "work/nih_denominator.parquet")
    audit["year"], audit["month"] = audit.oa_year.astype("Int64"), pd.to_datetime(audit.oa_pub_date).dt.month
    audit = audit[
        audit.oa_year.eq(audit.crossref_year) & audit.oa_nih.eq(True) & audit.pmid.notna() & audit.pmid.ne("")
    ]
    audit = audit[audit.year.isin([2025, 2026])].merge(
        matched[["journal_id", "month"]].drop_duplicates(), on=["journal_id", "month"], validate="many_to_one"
    )
    audit["pubmed_cohort"] = audit.pmid.isin(set(denominator.pmid.astype(str)))
    metadata = efetch(audit.pmid.astype(str).tolist())[["pmid", "pub_date"]]
    audit = audit[audit.pmid.isin(set(metadata.pmid.astype(str)))]
    audit = audit.merge(metadata, on="pmid", how="left", validate="one_to_one")
    mapped, _ = pmc_map(audit.pmid.astype(str))
    audit = audit.merge(mapped, on="pmid", how="left", validate="one_to_one")
    cols = [
        "pmid",
        "doi",
        "journal_id",
        "is_in_doaj",
        "year",
        "month",
        "pub_date",
        "oa_pub_date",
        "pmcid",
        "pubmed_cohort",
    ]
    audit[cols].to_csv(DATA / "funder_replication_articles.csv", index=False)
    print(f"funder-defined NIH articles: {len(audit)}")


if __name__ == "__main__":
    main()
