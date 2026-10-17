#!/usr/bin/env python3
"""Build the NIH-supported cohort, the comparison candidates, and the matched samples.

PubMed (E-utilities) gives the NIH-supported articles and the candidates' publication types and grant
agencies; a local OpenAlex works snapshot (OPENALEX_WORKS) gives journals, DOAJ status, dates and the
comparison pool; the PMC identifier list gives PMCIDs. Writes data/matched_articles.csv,
data/placebo_articles.csv and data/nih_cohort.csv; with COHORT_SUFFIX=_oct it writes
data/matched_articles_oct.csv from a later retrieval. Outcomes come from src/acquire/pmc_oai.py.
"""

import gzip
import hashlib
import http.client
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parents[1] / "data"
SUFFIX = os.environ.get("COHORT_SUFFIX", "")  # "_oct" for the re-match on a later PubMed retrieval
WORK = ROOT / ("work" + SUFFIX)
WORKS = os.environ.get("OPENALEX_WORKS", "openalex/works/*.parquet")  # local OpenAlex works snapshot
PMC_IDS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz"
SUPPORT = (
    '("Research Support, N.I.H., Extramural"[pt] OR "Research Support, N.I.H., Intramural"[pt])'
    " AND journal article[pt]"
)
SEED = 20260912
HORIZONS = (0, 7, 30, 180)
LAST_REQUEST = 0.0


def stable(value):
    return hashlib.sha256(f"{SEED}:{value}".encode()).hexdigest()


def request(url, data=None):
    global LAST_REQUEST
    is_oai = "pmc.ncbi.nlm.nih.gov/api/oai" in url
    attempts = 6 if is_oai else 3
    interval = 0.8 if is_oai else 0.4
    for attempt in range(attempts):
        wait = interval - (time.monotonic() - LAST_REQUEST)
        if wait > 0:
            time.sleep(wait)
        req = urllib.request.Request(
            url, data=data, headers={"User-Agent": "nih_embargo/1.0", "Accept-Encoding": "gzip"}
        )
        try:
            with urllib.request.urlopen(req, timeout=90) as response:
                body = response.read()
                if response.headers.get("Content-Encoding", "").lower() == "gzip":
                    body = gzip.decompress(body)
            LAST_REQUEST = time.monotonic()
            return body
        except urllib.error.HTTPError as error:
            if error.code != 429 and error.code < 500:
                raise
            if attempt == attempts - 1:
                raise
            if error.code == 429:
                time.sleep(float(error.headers.get("Retry-After", 30)))
                continue
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead):
            if attempt == attempts - 1:
                raise
        time.sleep(2**attempt)
    raise RuntimeError(f"unreachable retry state for {url}")


MONTHS = {
    name.lower(): i
    for i, name in enumerate(
        ("", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    )
}


def publication_date_parts(article):
    candidates = article.findall("./MedlineCitation/Article/ArticleDate")
    candidates += article.findall("./MedlineCitation/Article/Journal/JournalIssue/PubDate")
    partial = (None, None)
    for node in candidates:
        values = {x.tag: (x.text or "").strip() for x in node}
        year = int(values["Year"]) if values.get("Year", "").isdigit() else None
        raw_month = values.get("Month", "").lower()
        month = int(raw_month) if raw_month.isdigit() else MONTHS.get(raw_month[:3])
        day = int(values["Day"]) if values.get("Day", "").isdigit() else None
        if year and month and day:
            return date(year, month, day), year, month
        if year and month and partial == (None, None):
            partial = (year, month)
        medline = values.get("MedlineDate", "")
        match = re.search(r"(19|20)\d{2}\s+([A-Za-z]{3})", medline)
        if match and partial == (None, None):
            partial = (int(match.group(0)[:4]), MONTHS.get(match.group(2).lower()))
    return None, partial[0], partial[1]


def parse_articles(xml, assigned=None):
    root = ET.fromstring(xml)
    rows = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext("./MedlineCitation/PMID")
        ids = {
            (x.attrib.get("IdType") or "").lower(): (x.text or "").strip()
            for x in article.findall("./PubmedData/ArticleIdList/ArticleId")
        }
        types = [
            (x.text or "").strip()
            for x in article.findall("./MedlineCitation/Article/PublicationTypeList/PublicationType")
        ]
        grants = [
            (x.findtext("Agency") or "").strip()
            for x in article.findall("./MedlineCitation/Article/GrantList/Grant")
        ]
        journal = article.findtext("./MedlineCitation/MedlineJournalInfo/ISSNLinking")
        if not pmid:
            continue
        pub_date, pub_year, pub_month = publication_date_parts(article)
        rows.append(
            {
                "pmid": pmid,
                "doi": ids.get("doi"),
                "issn_l": journal,
                "pub_date": pub_date,
                "pub_year": pub_year,
                "pub_month": pub_month,
                "publication_types": types,
                "grant_agencies": grants,
                "query_month": assigned.get(pmid) if assigned else None,
            }
        )
    return rows


def efetch(pmids, assigned=None):
    rows = []
    for start in range(0, len(pmids), 500):
        form = urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(pmids[start : start + 500]), "retmode": "xml"}
        ).encode()
        for attempt in range(3):
            xml = request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", form)
            try:
                rows.extend(parse_articles(xml, assigned))
                break
            except ET.ParseError:
                if attempt == 2:
                    raise
                time.sleep(2**attempt)
    return pd.DataFrame(rows)


def search_ids(term):
    query = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmode": "json", "retmax": 9_999})
    result = json.loads(request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + query))[
        "esearchresult"
    ]
    return int(result["count"]), result["idlist"]


def query_pmids(term, low=0, high=99_999_999):
    bounded = f"{term} AND {low}:{high}[uid]"
    n, pmids = search_ids(bounded)
    if n > 9_999:
        midpoint = (low + high) // 2
        pmids = query_pmids(term, low, midpoint) + query_pmids(term, midpoint + 1, high)
    elif n != len(pmids):
        raise ValueError(f"PubMed ID capture incomplete for UID {low}:{high}: count={n}, ids={len(pmids)}")
    if len(pmids) != len(set(pmids)):
        raise ValueError(f"duplicate PMID inside UID partition {low}:{high}")
    return pmids


def nih_universe():
    counts = {}
    frames = []
    for year in (2024, 2025, 2026):
        term = f'{SUPPORT} AND ("{year:04d}/01/01"[dp] : "{year:04d}/04/30"[dp])'
        pmids = query_pmids(term)
        rows = efetch(pmids)
        if len(pmids) != len(rows):
            raise ValueError(
                f"PubMed explicit-ID fetch incomplete for {year}: ids={len(pmids)}, rows={len(rows)}"
            )
        counts[f"{year}-01-01_to_04-30"] = len(pmids)
        rows["query_year"] = year
        frames.append(rows)
    frame = pd.concat(frames, ignore_index=True)
    duplicate_hits = int(frame.pmid.duplicated().sum())
    doi_conflicts = frame.dropna(subset=["doi"]).groupby("pmid").doi.nunique()
    if (doi_conflicts > 1).any():
        raise ValueError(f"duplicate PMIDs have conflicting DOIs: {(doi_conflicts > 1).sum()}")
    year_hits = frame.groupby("pmid").query_year.agg(lambda x: ",".join(map(str, sorted(set(x)))))
    frame = frame.sort_values(["pmid", "query_year"]).drop_duplicates("pmid").copy()
    frame["query_year_hits"] = frame.pmid.map(year_hits)
    counts["cross_year_duplicate_hits"] = duplicate_hits
    counts["unique_pmids"] = len(frame)
    frame["nih"] = True
    return frame, counts


def normalize_doi(series):
    return series.astype("string").str.lower().str.replace(r"^https?://doi.org/", "", regex=True).str.strip()


def attach_openalex(nih):
    sample = nih[["pmid", "doi"]].copy()
    sample["doi_norm"] = normalize_doi(sample.doi)
    con = duckdb.connect()
    con.execute("SET threads=32")
    con.execute("SET memory_limit='150GB'")
    con.register("sample", sample)
    hits = con.execute(f"""
      WITH x AS (
        SELECT s.pmid sample_pmid, w.id openalex_id, w.publication_date oa_pub_date,
               w.primary_location.source.id journal_id,
               w.primary_location.source.is_in_doaj is_in_doaj,
               regexp_extract(w.ids.pmid, '([0-9]+)$', 1) oa_pmid,
               count(*) OVER (PARTITION BY s.pmid) doi_n,
               sum((regexp_extract(w.ids.pmid, '([0-9]+)$', 1)=s.pmid)::INTEGER)
                 OVER (PARTITION BY s.pmid) pmid_match_n
        FROM read_parquet('{WORKS}', union_by_name=true) w
        JOIN sample s ON lower(regexp_replace(w.doi, '^https://doi.org/', ''))=s.doi_norm
        WHERE s.doi_norm IS NOT NULL
      )
      SELECT sample_pmid, openalex_id, oa_pub_date, journal_id, is_in_doaj FROM x
      WHERE doi_n=1 OR (pmid_match_n=1 AND oa_pmid=sample_pmid)
    """).df()
    con.close()
    if hits.sample_pmid.duplicated().any():
        raise ValueError("OpenAlex DOI+PMID resolution is not unique")
    out = nih.merge(hits, how="left", left_on="pmid", right_on="sample_pmid", validate="one_to_one").drop(
        columns="sample_pmid"
    )
    out["pub_date"] = pd.to_datetime(out.pub_date).fillna(pd.to_datetime(out.oa_pub_date))
    out["year"] = out.pub_year.fillna(out.pub_date.dt.year).astype("Int64")
    out["month"] = out.pub_month.fillna(out.pub_date.dt.month).astype("Int64")
    return out


def control_candidates(nih):
    valid = nih.dropna(subset=["journal_id"])
    needs = valid.groupby(["journal_id", "year", "month"]).size().rename("nih_n").reset_index()
    treated = valid[["pmid"]].drop_duplicates()
    candidate_path = WORK / "control_candidates_openalex.parquet"
    if candidate_path.exists():
        candidates = pd.read_parquet(candidate_path)
    else:
        con = duckdb.connect()
        con.execute("SET threads=32")
        con.execute("SET memory_limit='150GB'")
        con.register("needs", needs)
        con.register("treated", treated)
        candidates = con.execute(f"""
      WITH pool AS (
        SELECT regexp_extract(w.ids.pmid, '([0-9]+)$', 1) pmid,
               w.doi, w.id openalex_id, w.publication_date oa_pub_date,
               w.primary_location.source.id journal_id,
               w.primary_location.source.is_in_doaj is_in_doaj,
               w.publication_year match_year, month(w.publication_date) match_month,
               row_number() OVER (
                 PARTITION BY w.primary_location.source.id, w.publication_year, month(w.publication_date)
                 ORDER BY hash(w.id, {SEED})
               ) candidate_rank
        FROM read_parquet('{WORKS}', union_by_name=true) w
        JOIN needs n ON w.primary_location.source.id=n.journal_id
                    AND w.publication_year=n."year" AND month(w.publication_date)=n."month"
        LEFT JOIN treated t ON regexp_extract(w.ids.pmid, '([0-9]+)$', 1)=t.pmid
        WHERE w.type='article' AND w.is_paratext=false AND coalesce(w.is_retracted,false)=false
          AND regexp_extract(w.ids.pmid, '([0-9]+)$', 1)<>'' AND t.pmid IS NULL
      )
      SELECT p.* FROM pool p JOIN needs n ON p.journal_id=n.journal_id
        AND p.match_year=n."year" AND p.match_month=n."month"
      WHERE candidate_rank<=greatest(12, 6*n.nih_n)
        """).df()
        con.close()
        candidates = candidates.rename(columns={"match_year": "year", "match_month": "month"})
        if candidates.pmid.duplicated().any():
            candidates = candidates.sort_values(["candidate_rank", "openalex_id"]).drop_duplicates("pmid")
        candidates.to_parquet(candidate_path, index=False)
    meta = efetch(candidates.pmid.tolist())
    meta["us_gov"] = meta.publication_types.map(
        lambda x: any("Research Support, U.S. Gov't" in v or "N.I.H." in v or "P.H.S." in v for v in x)
    )
    meta["us_gov"] |= meta.grant_agencies.map(
        lambda x: any(re.search(r"\b(NIH|NCI|NHLBI|NIAID|NIGMS|NSF)\b", v, re.I) for v in x)
    )
    controls = candidates.merge(meta.drop(columns=["query_month"]), on="pmid", validate="one_to_one")
    controls["pub_date"] = pd.to_datetime(controls.pub_date).fillna(pd.to_datetime(controls.oa_pub_date))
    controls = controls[~controls.us_gov & controls.pub_date.notna()].copy()
    controls["nih"] = False
    return controls


def select_comparison(nih, controls, years, name):
    a, b = years
    selected = []
    eligible = nih.year.isin(years) & nih.month.isin([1, 2, 3, 4]) & nih.pub_date.notna()
    for (journal, month), group in nih[eligible].groupby(["journal_id", "month"]):
        na, nb = (group.year == a).sum(), (group.year == b).sum()
        ca = controls[(controls.journal_id == journal) & (controls.year == a) & (controls.month == month)]
        cb = controls[(controls.journal_id == journal) & (controls.year == b) & (controls.month == month)]
        n = min(na, nb, len(ca) // 2, len(cb) // 2)
        if n == 0:
            continue
        for source in (group[group.year == a], group[group.year == b]):
            take = (
                source.assign(_key=source.pmid.map(lambda x: stable(f"{name}:nih:{x}")))
                .sort_values("_key")
                .head(n)
                .drop(columns="_key")
            )
            take["comparison"], take["control_rank"] = name, 0
            selected.append(take)
        for source in (ca, cb):
            take = (
                source.assign(_key=source.pmid.map(lambda x: stable(f"{name}:control:{x}")))
                .sort_values("_key")
                .head(2 * n)
                .drop(columns="_key")
            )
            take["comparison"] = name
            take["control_rank"] = np.tile((1, 2), n)
            selected.append(take)
    if not selected:
        raise ValueError(f"no common-support cells for {name}")
    result = pd.concat(selected, ignore_index=True)
    counts = result.groupby(["journal_id", "month", "year", "nih"]).size().unstack(fill_value=0)
    if (counts.min(axis=1) == 0).any():
        raise ValueError(f"unbalanced treatment support in {name}")
    return result


def pmc_map(pmids):
    path = WORK / "pmc_ids.csv.gz"
    if not path.exists():
        with urllib.request.urlopen(PMC_IDS_URL, timeout=300) as response, path.open("wb") as target:
            while chunk := response.read(1024 * 1024):
                target.write(chunk)
    con = duckdb.connect()
    con.register("wanted", pd.DataFrame({"pmid": list(set(pmids))}))
    mapped = con.execute(f"""
      SELECT cast(ids."PMID" AS VARCHAR) pmid, ids."PMCID" pmcid
      FROM read_csv_auto('{path}', compression='gzip', all_varchar=true, delim=',', quote='"', escape='"') ids
      JOIN wanted w ON cast(ids."PMID" AS VARCHAR)=w.pmid
    """).df()
    con.close()
    if mapped.pmid.duplicated().any():
        raise ValueError("PMC bulk map has duplicate PMID")
    return mapped, path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if not list(Path(WORKS.split("updated_date=")[0]).glob("updated_date=*/*.parquet")):
        raise FileNotFoundError("OpenAlex 2026-06-26 snapshot is missing")
    WORK.mkdir(exist_ok=True)
    metadata_path, metadata_info_path = WORK / "design_metadata.parquet", WORK / "design_metadata.json"
    if metadata_path.exists() != metadata_info_path.exists():
        raise ValueError("design metadata checkpoint is incomplete")
    if metadata_path.exists():
        metadata = pd.read_parquet(metadata_path)
        metadata_info = json.loads(metadata_info_path.read_text())
        query_counts = metadata_info["query_counts"]
        placebo, primary = (
            metadata[metadata.comparison == "placebo"],
            metadata[metadata.comparison == "primary"],
        )
    else:
        raw_path, raw_info_path = WORK / "nih_pubmed_raw.parquet", WORK / "nih_pubmed_raw.json"
        if raw_path.exists() != raw_info_path.exists():
            raise ValueError("raw PubMed checkpoint is incomplete")
        if raw_path.exists():
            nih = pd.read_parquet(raw_path)
            query_counts = json.loads(raw_info_path.read_text())["query_counts"]
        else:
            nih, query_counts = nih_universe()
            nih.to_parquet(raw_path, index=False)
            raw_info_path.write_text(
                json.dumps({"query_counts": query_counts, "sha256": sha(raw_path)}, indent=2) + "\n"
            )
        nih = attach_openalex(nih)
        if nih.openalex_id.notna().mean() < 0.90:
            raise ValueError(f"OpenAlex linkage below 90%: {nih.openalex_id.notna().mean():.3f}")
        denominator_path = WORK / "nih_denominator.parquet"
        nih.to_parquet(denominator_path, index=False)
        denominator_audit = [
            {
                "canonical_year": "all",
                "captured": len(nih),
                "doi_present": int(nih.doi.notna().sum()),
                "openalex_linked": int(nih.openalex_id.notna().sum()),
                "exact_or_openalex_date": int(nih.pub_date.notna().sum()),
                "eligible_jan_apr": int(
                    (
                        nih.year.isin([2024, 2025, 2026])
                        & nih.month.isin([1, 2, 3, 4])
                        & nih.pub_date.notna()
                    ).sum()
                ),
            }
        ]
        for year in (2024, 2025, 2026):
            group = nih[nih.year.eq(year)]
            denominator_audit.append(
                {
                    "canonical_year": year,
                    "captured": len(group),
                    "doi_present": int(group.doi.notna().sum()),
                    "openalex_linked": int(group.openalex_id.notna().sum()),
                    "exact_or_openalex_date": int(group.pub_date.notna().sum()),
                    "eligible_jan_apr": int(
                        (group.year.eq(year) & group.month.isin([1, 2, 3, 4]) & group.pub_date.notna()).sum()
                    ),
                }
            )
        controls = control_candidates(nih)
        placebo = select_comparison(nih, controls, (2024, 2025), "placebo")
        primary = select_comparison(nih, controls, (2025, 2026), "primary")
        metadata = pd.concat([placebo, primary], ignore_index=True, sort=False)
        metadata.to_parquet(metadata_path, index=False)
        metadata_info = {
            "query_counts": query_counts,
            "nih_n": len(nih),
            "control_candidates_n": len(controls),
            "denominator_audit": denominator_audit,
            "denominator_sha256": sha(denominator_path),
            "metadata_sha256": sha(metadata_path),
        }
        metadata_info_path.write_text(json.dumps(metadata_info, indent=2) + "\n")
    ids, _ = pmc_map(metadata.pmid.astype(str))
    tables = metadata.merge(ids, on="pmid", how="left", validate="many_to_one")
    tables["pub_date"] = pd.to_datetime(tables.pub_date).dt.date
    join = lambda x: "; ".join(x) if hasattr(x, "__len__") and not isinstance(x, str) else x  # noqa: E731
    for col in ("publication_types", "grant_agencies"):
        tables[col] = tables[col].map(join)
    primary = tables[(tables.comparison == "primary") & (tables.control_rank <= 1)]
    placebo = tables[(tables.comparison == "placebo") & (tables.control_rank <= 1)]
    keep = ["pmid", "doi", "journal_id", "is_in_doaj", "year", "month", "nih", "pub_date", "oa_pub_date"]
    primary[keep + ["publication_types", "grant_agencies", "pmcid"]].to_csv(
        DATA / f"matched_articles{SUFFIX}.csv", index=False
    )
    if not SUFFIX:
        placebo[[c for c in keep if c != "doi"] + ["pmcid"]].to_csv(
            DATA / "placebo_articles.csv", index=False
        )
        nih = pd.read_parquet(WORK / "nih_denominator.parquet")
        nih = nih[nih.year.isin([2025, 2026]) & nih.month.isin([1, 2, 3, 4]) & nih.pub_date.notna()]
        nih = nih.merge(pmc_map(nih.pmid.astype(str))[0], on="pmid", how="left")
        nih["pub_date"] = pd.to_datetime(nih.pub_date).dt.date
        nih[[c for c in keep if c != "nih"] + ["pmcid"]].to_csv(DATA / "nih_cohort.csv", index=False)
    print(f"BUILD COMPLETE: nih={metadata_info['nih_n']}, controls={metadata_info['control_candidates_n']}")


if __name__ == "__main__":
    main()
