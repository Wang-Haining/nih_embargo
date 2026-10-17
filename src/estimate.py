"""Estimate the effect of the 2025 NIH public-access policy on availability in PubMed Central.

Availability is dated by the day each article went live in PMC (the pmc-live event in the PMC OAI-PMH
front matter); the release date PMC reports is used as a comparison outcome. Publication date is the
earlier of the PubMed and OpenAlex dates. Reads the article-level tables in data/ and writes every reported
estimate to results/.

Run from the repository root:  python src/estimate.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260912
BOOTSTRAP_DRAWS = 2000
HORIZONS = (0, 7, 30, 60, 90, 120)
CURVE_DAYS = 120
AS_OF = pd.Timestamp("2026-10-16")
POLICY = pd.Timestamp("2025-07-01")
IN = Path("data")
OUT = Path("results")
rng = np.random.default_rng(SEED)


def reseed() -> None:
    """Each analysis stage starts from the same seed, so stages can be rerun independently."""
    global rng
    rng = np.random.default_rng(SEED)


# ---------------------------------------------------------------- estimator


def cell_did(frame: pd.DataFrame, outcome: str, years: tuple[int, int]) -> tuple[pd.Series, pd.Series]:
    """Per-cell difference-in-differences and cell weights (matched NIH articles in the first year)."""
    before, after = years
    means = frame.groupby(["journal_id", "month", "year", "nih"])[outcome].mean().unstack(["year", "nih"])
    means = means.dropna()
    weights = (
        frame[frame["nih"] & (frame["year"] == before)]
        .groupby(["journal_id", "month"])
        .size()
        .reindex(means.index)
        .fillna(0)
    )
    did = (means[(after, True)] - means[(after, False)]) - (means[(before, True)] - means[(before, False)])
    return did, weights


def journal_bootstrap(did: pd.Series, weights: pd.Series) -> tuple[float, float, float, np.ndarray]:
    journals = did.index.get_level_values(0).unique()
    numer = did.groupby(level=0).apply(lambda s: (s * weights[s.index]).sum())
    denom = weights.groupby(level=0).sum()
    draws = np.empty(BOOTSTRAP_DRAWS)
    for b in range(BOOTSTRAP_DRAWS):
        sample = rng.choice(journals, len(journals))
        draws[b] = numer[sample].sum() / denom[sample].sum()
    point = float(np.average(did, weights=weights))
    lower, upper = np.percentile(draws, [2.5, 97.5])
    return point, float(lower), float(upper), draws


def did(frame: pd.DataFrame, outcome: str, years=(2025, 2026), weights_cap: float | None = None) -> dict:
    """Weighted journal-month DID in percentage points with a journal-cluster bootstrap interval."""
    d, w = cell_did(frame, outcome, years)
    if weights_cap is not None:
        journal_weight = w.groupby(level=0).transform("sum")
        w = w * np.minimum(1.0, weights_cap * w.sum() / journal_weight)
    point, lower, upper, draws = journal_bootstrap(d, w)
    return {
        "est": 100 * point,
        "lo": 100 * lower,
        "hi": 100 * upper,
        "se": 100 * draws.std(ddof=1),
        "journals": d.index.get_level_values(0).nunique(),
        "cells": len(d),
        "w": float(w.sum()),
    }


def point_did(frame: pd.DataFrame, outcome: str, years=(2025, 2026)) -> float:
    d, w = cell_did(frame, outcome, years)
    return float(100 * np.average(d, weights=w))


def curves_with_bands(frame: pd.DataFrame, days: int) -> pd.DataFrame:
    """Share live in PMC by each day after publication, by DOAJ stratum, year and group, with bands."""
    grid = list(range(days + 1))
    indicator = np.stack(
        [(frame["pmcid"].notna() & frame["lag_live"].le(d)).to_numpy(float) for d in grid], axis=1
    )
    rows = []
    for (in_doaj, year, nih), group in frame.groupby(["is_in_doaj", "year", "nih"]):
        values = indicator[frame.index.get_indexer(group.index)]
        journal = group["journal_id"].to_numpy()
        unique = np.unique(journal)
        sums = np.stack([values[journal == j].sum(0) for j in unique])
        counts = np.array([(journal == j).sum() for j in unique])
        draws = np.empty((BOOTSTRAP_DRAWS, len(grid)))
        for b in range(BOOTSTRAP_DRAWS):
            sample = rng.integers(0, len(unique), len(unique))
            draws[b] = sums[sample].sum(0) / counts[sample].sum()
        lower, upper = np.percentile(draws, [2.5, 97.5], axis=0)
        point = values.mean(0)
        for k, day in enumerate(grid):
            rows.append(
                {
                    "in_doaj": bool(in_doaj),
                    "year": int(year),
                    "nih": bool(nih),
                    "day": day,
                    "share": point[k],
                    "ci_lower": lower[k],
                    "ci_upper": upper[k],
                    "n": len(group),
                    "journals": len(unique),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- data


def load_records() -> tuple[pd.DataFrame, set]:
    rec = pd.read_csv(IN / "pmc_records.csv")
    unavailable = set(rec.loc[~rec["record_available"], "pmcid"])
    rec = rec[rec["record_available"]].drop(columns="record_available")
    for c in ("pmc_release_date", "pmc_live_date", "nihms_submitted_date", "accepted_date"):
        rec[c] = pd.to_datetime(rec[c])
    return rec, unavailable


def attach(frame: pd.DataFrame, rec: pd.DataFrame, unavailable: set, agreements: pd.Series) -> pd.DataFrame:
    """Join PMC records and derive publication date, lags, outcomes by horizon, route and agreement."""
    f = frame[~frame["pmcid"].isin(unavailable)]  # records the OAI-PMH service did not return
    f = f.merge(rec, on="pmcid", how="left")
    f["early"] = pd.concat(
        [pd.to_datetime(f["pub_date"], errors="coerce"), pd.to_datetime(f["oa_pub_date"], errors="coerce")],
        axis=1,
    ).min(axis=1)
    f["release"] = f["pmc_release_date"]
    fallback = f["status_live"].eq("yes") & f["pmc_live_date"].isna()
    f["public"] = f["pmc_live_date"].where(~fallback, f["release"])
    f["lag_release"] = (f["release"] - f["early"]).dt.days
    f["lag_live"] = (f["public"] - f["early"]).dt.days
    f["follow_up"] = (AS_OF - f["early"]).dt.days
    has = f["pmcid"].notna()
    for h in HORIZONS:
        f[f"rel_{h}d"] = (has & f["lag_release"].le(h)).astype(float)
        f[f"live_{h}d"] = (has & f["lag_live"].le(h)).astype(float)
    f["route"] = np.where(~has, "none", np.where(f["author_manuscript"].eq("yes"), "manuscript", "final"))
    f["agreement"] = f["journal_id"].map(agreements).fillna("none")
    f["covered"] = f["accepted_date"].ge(POLICY)
    return f


# ---------------------------------------------------------------- go-live date checks


def restamp_flags(rec: pd.DataFrame) -> pd.DataFrame:
    """A record is re-stamped if it went live >14 days after its release date and its 10 PMCID neighbours."""
    o = rec[rec["pmc_live_date"].notna()].copy()
    o["n"] = o["pmcid"].str[3:].astype(int)
    o = o.sort_values("n").reset_index(drop=True)
    live = o["pmc_live_date"].astype("datetime64[ns]").astype("int64")
    neigh = [
        np.median(live.iloc[[k for k in range(max(0, i - 5), min(len(o), i + 6)) if k != i]])
        for i in range(len(o))
    ]
    o["neigh"] = pd.to_datetime(np.array(neigh, dtype="int64"), unit="ns")
    o["restamped_without_release"] = (o["pmc_live_date"] - o["neigh"]).dt.days > 14
    later = pd.concat([o["neigh"], o["pmc_release_date"]], axis=1).max(axis=1)
    o["restamped"] = (o["pmc_live_date"] - later).dt.days > 14
    return o[["pmcid", "restamped", "restamped_without_release"]]


# ---------------------------------------------------------------- articles on terminated NIH awards


def terminated_award_articles(prim: pd.DataFrame, plac: pd.DataFrame) -> dict:
    """Availability of NIH-supported articles that list an NIH award terminated since 2025.

    Compared within journal-month cells, in 2026 (after most terminations) and in 2024 (the same kind
    of award, before any termination); the difference between the two is the change. Journal-cluster
    bootstrap.
    """
    grants = pd.read_csv(IN / "pubmed_grants.csv", dtype=str)
    agency = grants["agency"].fillna("")
    nih_agency = agency.str.contains(
        r"\bNIH\b|National Institutes? (?:of|on)|National (?:Cancer|Eye|Heart|Human|Library|Center)|Fogarty|"
        r"Office of the Director|\bOBSSR\b|\bNCATS\b"
    ) & ~agency.str.contains("Foundation for the")
    keys = grants["grant_id"].fillna("").str.upper().str.findall(r"(?<![A-Z])([A-Z]{2})\s*-?\s*(\d{6})(?!\d)")
    grants["key"] = keys.apply(lambda ks: [a + b for a, b in ks])  # every award key in the grant field
    grants = grants[nih_agency].explode("key").dropna(subset=["key"])
    term = pd.read_csv(IN / "grant_terminations.csv", dtype=str).set_index("award_key")
    grants = grants.join(term, on="key")
    listed = grants.groupby("pmid")["earliest_disruption"].apply(lambda d: d.notna().any())
    arts = pd.concat([prim[prim["nih"] & (prim["year"] == 2026)], plac[plac["nih"] & (plac["year"] == 2024)]])
    arts = arts[arts["pmid"].isin(listed.index)].copy()  # articles listing at least one NIH award number
    arts["terminated"] = arts["pmid"].map(listed).astype(float)
    arts["cell"] = arts["journal_id"] + "|" + arts["year"].astype(str) + "|" + arts["month"].astype(str)
    out: dict = {}
    # 2026 articles listing an award terminated before publication, and one still terminated on that date
    linked = grants[grants["earliest_disruption"].notna()].merge(arts[["pmid", "year", "early"]], on="pmid")
    linked = linked[linked["year"] == 2026]
    start = pd.to_datetime(linked["earliest_disruption"])
    restored = pd.to_datetime(linked["latest_restoration"])
    linked["before"] = start < linked["early"]
    linked["at_pub"] = linked["before"] & (restored.isna() | (restored > linked["early"]))
    out["2026_terminated_before_pub_n"] = int(linked.groupby("pmid")["before"].any().sum())
    out["2026_terminated_at_pub_n"] = int(linked.groupby("pmid")["at_pub"].any().sum())
    for y, g in arts.groupby("year"):
        out[f"{y}_n"] = len(g)
        out[f"{y}_terminated_n"] = int(g["terminated"].sum())
        for t, gg in g.groupby("terminated"):
            out[f"{y}_{'terminated' if t else 'other'}_live_120d"] = round(float(gg["live_120d"].mean()), 6)
    for h in (7, 120):
        y_col = f"live_{h}d"
        d = arts[["journal_id", "year", "cell", "terminated", y_col]].copy()
        for c in ("terminated", y_col):  # deviations from cell means: within-cell comparison
            d[c] = d[c] - d.groupby("cell")[c].transform("mean")
        d["sxy"], d["sxx"] = d["terminated"] * d[y_col], d["terminated"] ** 2
        sums = d.groupby(["journal_id", "year"])[["sxy", "sxx"]].sum().unstack("year", fill_value=0.0)
        journals = sums.index.to_numpy()

        def slopes(sample, sums=sums):
            s = sums.loc[sample].sum()
            return 100 * s["sxy"][2026] / s["sxx"][2026], 100 * s["sxy"][2024] / s["sxx"][2024]

        b26, b24 = slopes(journals)
        draws = np.array([slopes(rng.choice(journals, len(journals))) for _ in range(BOOTSTRAP_DRAWS)])
        for name, point, col in (
            ("2026", b26, draws[:, 0]),
            ("2024", b24, draws[:, 1]),
            ("change", b26 - b24, draws[:, 0] - draws[:, 1]),
        ):
            lo, hi = np.percentile(col, [2.5, 97.5])
            out[f"{name}_{h}d"] = {
                "est": round(float(point), 6),
                "lo": round(float(lo), 6),
                "hi": round(float(hi), 6),
            }
    ms = arts[(arts["year"] == 2026) & arts["route"].eq("manuscript")]
    lag = (ms["nihms_submitted_date"] - ms["early"]).dt.days
    for t, g in lag.groupby(ms["terminated"]):
        tag = "terminated" if t else "other"
        out[f"ms26_{tag}_n"] = int(g.size)
        out[f"ms26_{tag}_pub_to_sub_q50"] = float(g.median())
    return out


# ---------------------------------------------------------------- analysis of one matched sample


def analyse(
    prim: pd.DataFrame, plac: pd.DataFrame, rec: pd.DataFrame, num: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def put(key, value):
        if isinstance(value, dict):
            value = {k: round(float(value[k]), 6) for k in ("est", "lo", "hi")}
        elif value is not None:
            value = round(float(value), 6)
        num[key] = value

    cov = prim["pmcid"].notna()

    # Check 1: re-stamping, and estimates excluding flagged records
    flags = restamp_flags(rec)
    fl = prim.merge(flags, on="pmcid", how="left")
    by_arm = (
        prim.merge(flags, on="pmcid", how="inner")
        .groupby(["year", "nih"])[["restamped_without_release", "restamped"]]
        .mean()
    )
    for rule, key in (("restamped_without_release", "check1_without_release"), ("restamped", "check1")):
        keep = ~fl[rule].fillna(False).astype(bool).to_numpy()
        for h in (7, 120):
            put(f"{key}_{h}d", did(prim[keep], f"live_{h}d"))
        put(f"{key}_excluded", int((~keep).sum()))
    pl = plac[plac["year"] == 2024]
    put(
        "check2_tail_share",
        ((pl["pmc_live_date"] - pl["pmc_release_date"]).dt.days > 180)[pl["pmc_live_date"].notna()].mean(),
    )
    for (y, n), v in by_arm["restamped"].items():
        put(f"check1_{y}_{'nih' if n else 'cmp'}", v)
    for (y, n), v in by_arm["restamped_without_release"].items():
        put(f"check1_without_release_{y}_{'nih' if n else 'cmp'}", v)

    # Release date earlier than going live by more than 7 days
    x = prim[cov & prim["public"].notna()]
    early_release = (x["public"] - x["release"]).dt.days > 7
    for (y, n), v in early_release.groupby([x["year"], x["nih"]]).mean().items():
        put(f"release_early_{y}_{'nih' if n else 'cmp'}", v)
    for (y, r), v in early_release.groupby([x["year"], x["route"]]).mean().items():
        put(f"release_early_{y}_{r}", v)
    final = x["route"] == "final"
    at_pub = (x["release"] - x["early"]).dt.days.abs() <= 3  # release date within 3 days of publication
    put("release_early_final_n", (final & early_release).sum())
    put("release_early_final_at_pub", at_pub[final & early_release].mean())
    put("release_other_final_at_pub", at_pub[final & ~early_release].mean())
    put("release_early_final_gap_q50", (x["public"] - x["release"]).dt.days[final & early_release].median())

    # DID by horizon and stratum, go-live and release dates
    rows = []
    strata = {"all": prim, "non_doaj": prim[~prim["is_in_doaj"]], "doaj": prim[prim["is_in_doaj"]]}
    for a in ("All articles", "NIH funded articles", "none"):
        strata[f"agreement={a}"] = prim[prim["agreement"].eq(a)]
    for s, f in strata.items():
        for h in HORIZONS:
            for kind in ("rel", "live"):
                if s.startswith("agreement") and h not in (7, 120):
                    continue
                r = did(f, f"{kind}_{h}d")
                m = f.groupby(["year", "nih"])[f"{kind}_{h}d"].mean() * 100
                rows.append(
                    {
                        "stratum": s,
                        "horizon": h,
                        "outcome": kind,
                        **r,
                        "nih25": m[(2025, True)],
                        "cmp25": m[(2025, False)],
                        "nih26": m[(2026, True)],
                        "cmp26": m[(2026, False)],
                    }
                )
    by_horizon = pd.DataFrame(rows)

    # Route split (author manuscript, final published version); the two parts sum to the total
    rr = []
    for h in (0, 7, 30, 120):
        for route in ("manuscript", "final"):
            col = f"{route}_{h}"
            prim[col] = prim[f"live_{h}d"] * prim["route"].eq(route)
            for s, f in (("all", prim), ("non_doaj", prim[~prim["is_in_doaj"]])):
                rr.append({"h": h, "route": route, "stratum": s, **did(f, col)})
    routes = pd.DataFrame(rr)

    # Timing by route; author manuscripts: publication -> NIH Manuscript Submission -> live
    nih = prim[prim["nih"] & cov]
    for (y, r), g in nih.groupby(["year", "route"]):
        put(f"n_{y}_{r}", len(g))
        put(f"live_n_{y}_{r}", g["public"].notna().sum())
        put(f"lag_live_median_{y}_{r}", g["lag_live"].median())
        put(f"live_day0_{y}_{r}", (g["lag_live"] <= 0).mean())
        put(f"live_wk_{y}_{r}", (g["lag_live"] <= 7).mean())
    ms = nih[(nih["year"] == 2026) & nih["route"].eq("manuscript")].copy()
    ms["pub_to_sub"] = (ms["nihms_submitted_date"] - ms["early"]).dt.days
    ms["sub_to_live"] = (ms["public"] - ms["nihms_submitted_date"]).dt.days
    put("ms26_n", len(ms))
    put("ms26_live_n", ms["public"].notna().sum())
    put("ms26_embargoed_n", ms["status_embargo"].eq("yes").sum())
    put("ms26_after_pub_n", (ms["pub_to_sub"] > 0).sum())
    put("ms26_on_or_before_pub_n", (ms["pub_to_sub"] <= 0).sum())
    for c in ("pub_to_sub", "sub_to_live", "lag_live"):
        for q in (0.25, 0.5, 0.75):
            put(f"ms26_{c}_q{int(q * 100)}", ms[c].quantile(q))
    put("ms26_after_pub", (ms["pub_to_sub"] > 0).mean())
    put("ms26_before_pub_n", (ms["pub_to_sub"] <= 0).sum())
    put("ms26_before_pub_wk", (ms.loc[ms["pub_to_sub"] <= 0, "lag_live"] <= 7).mean())
    monthly = ms[ms["early"].dt.year == 2026].groupby(ms["early"].dt.month)["sub_to_live"].median()
    put("ms26_sub_to_live_monthly_min", monthly.min())
    put("ms26_sub_to_live_monthly_max", monthly.max())
    acc = pd.read_csv(IN / "pubmed_accepted.csv", dtype=str)
    k = ms.merge(acc, on="pmid", how="left")
    k["acc_pm"] = pd.to_datetime(k["pubmed_accepted_date"])
    k = k[k["acc_pm"] >= POLICY]
    k["acc_to_sub"] = (k["nihms_submitted_date"] - k["acc_pm"]).dt.days
    k["acc_to_pub"] = (k["early"] - k["acc_pm"]).dt.days
    put("acc_n", len(k))
    for c in ("acc_to_pub", "acc_to_sub", "pub_to_sub", "sub_to_live"):
        for q in (0.25, 0.5, 0.75):
            put(f"acc_{c}_q{int(q * 100)}", k[c].quantile(q))
    put("acc_within7", (k["acc_to_sub"] <= 7).mean())
    put("acc_before_pub", (k["pub_to_sub"] <= 0).mean())
    # Manuscripts from publishers on NIH's Method D list (publisher-initiated deposit), by DOI registrant
    pubs = pd.read_csv(IN / "doi_publishers.csv", dtype={"prefix": str}).set_index("prefix")
    prefix = ms["doi"].str.replace("https://doi.org/", "", regex=False).str.split("/").str[0]
    known = prefix.map(pubs["publisher"]).notna()
    put("ms26_publisher_unknown_n", (~known).sum())
    method_d = prefix[known].map(pubs["method_d"]).astype(bool)
    for flag, g in ms[known].groupby(method_d):
        tag = "method_d" if flag else "other"
        put(f"ms26_{tag}_n", len(g))
        put(f"ms26_{tag}_pub_to_sub_q50", g["pub_to_sub"].median())
        put(f"ms26_{tag}_after_pub", (g["pub_to_sub"] > 0).mean())
    ms5 = nih[(nih["year"] == 2025) & nih["route"].eq("manuscript")].copy()
    ms5["pub_to_sub"] = (ms5["nihms_submitted_date"] - ms5["early"]).dt.days
    put("lag_live_median_2026_all", nih.loc[nih["year"] == 2026, "lag_live"].median())
    put("ms25_pub_to_sub_q50", ms5["pub_to_sub"].median())
    put("ms25_lag_live_q50", ms5["lag_live"].median())
    put("comparison_manuscripts_in_pmc", int((~prim["nih"] & prim["route"].eq("manuscript")).sum()))
    put("comparison_in_pmc", int((~prim["nih"] & cov).sum()))

    # Pre-policy contrast, 2024 vs 2025: all journals, the primary journals, the primary cells
    pj = set(prim["journal_id"])
    pc = set(zip(prim["journal_id"], prim["month"], strict=True))
    incell = pd.Series(
        [c in pc for c in zip(plac["journal_id"], plac["month"], strict=True)], index=plac.index
    )
    for name, f in (
        ("all", plac),
        ("same_journals", plac[plac["journal_id"].isin(pj)]),
        ("same_cells", plac[incell]),
    ):
        for h in (7, 30, 120):
            for kind in ("rel", "live"):
                r = did(f, f"{kind}_{h}d", (2024, 2025))
                put(f"placebo_{name}_{h}d_{kind}", r)
                num[f"placebo_{name}_journals"] = r["journals"]
                num[f"placebo_{name}_cells"] = r["cells"]

    # Concentration: leave one journal out, and per-journal weight capped at 2%
    for kind in ("rel_7", "live_7", "rel_120", "live_120"):
        d, w = cell_did(prim, f"{kind}d", (2025, 2026))
        loo = [
            100 * np.average(d[d.index.get_level_values(0) != j], weights=w[d.index.get_level_values(0) != j])
            for j in d.index.get_level_values(0).unique()
        ]
        put(f"lojo_min_{kind}", min(loo))
        put(f"lojo_max_{kind}", max(loo))
        put(f"capped_{kind}", did(prim, f"{kind}d", weights_cap=0.02))
    return by_horizon, routes


def main() -> None:
    OUT.mkdir(exist_ok=True)
    rec, unavailable = load_records()
    agreements = pd.read_csv(IN / "journal_agreements.csv").set_index("journal_id")["deposit_agreement"]
    num: dict = {"pmc_records_unavailable": len(unavailable)}
    check = pd.read_csv(IN / "pmc_live_check.csv")  # PMC Live Date search field against the OAI-PMH dates
    num["live_check"] = {
        stratum: {
            "n": len(g),
            "equals_live": int(g["pmc_live_date_equals_live"].sum()),
            "equals_release": int(g["pmc_live_date_equals_release"].sum()),
        }
        for stratum, g in check.groupby("stratum")
    }

    def load(name: str) -> pd.DataFrame:
        return attach(pd.read_csv(IN / name, dtype={"pmid": str}), rec, unavailable, agreements)

    prim, plac = load("matched_articles.csv"), load("placebo_articles.csv")
    num["matched_rows"] = len(prim)
    num["analysed_n"] = {
        f"{y}_{'nih' if n else 'cmp'}": int(v) for (y, n), v in prim.groupby(["year", "nih"]).size().items()
    }
    num["analysed_n_non_doaj"] = {
        f"{y}_{'nih' if n else 'cmp'}": int(v)
        for (y, n), v in prim[~prim["is_in_doaj"]].groupby(["year", "nih"]).size().items()
    }
    num["min_follow_up_days"] = int(prim["follow_up"].min())

    # Primary sample
    reseed()
    by_horizon, routes = analyse(prim, plac, rec, num)
    by_horizon.to_csv(OUT / "did_by_horizon.csv", index=False)
    routes.to_csv(OUT / "route_decomposition.csv", index=False)

    # Funder-metadata replication: NIH support from Crossref and OpenAlex; controls lack Crossref NIH funding
    audit = load("funder_replication_articles.csv").assign(nih=True)
    crossref = set(pd.read_csv(IN / "crossref_nih_dois.csv")["doi"].str.lower())
    doi = prim["doi"].astype("string").str.lower().str.replace(r"^https?://doi.org/", "", regex=True)
    contaminated = (~prim["nih"]) & doi.isin(crossref)
    rep = pd.concat([audit, prim[~prim["nih"] & ~contaminated]], ignore_index=True, sort=False)
    for h in (7, 30, 120):
        for kind in ("rel", "live"):
            r = did(rep, f"{kind}_{h}d")
            num[f"crossref_{h}d_{kind}"] = {k: round(float(r[k]), 6) for k in ("est", "lo", "hi")}
            num["crossref_journals"], num["crossref_cells"] = r["journals"], r["cells"]
    d_rep, _ = cell_did(rep, "live_7d", (2025, 2026))
    in_cells = rep.set_index(["journal_id", "month"]).index.isin(d_rep.index)
    num["crossref_counts"] = {
        f"{y}_{'nih' if n else 'cmp'}": int(v)
        for (y, n), v in rep[in_cells].groupby(["year", "nih"]).size().items()
    }
    contamination = (
        prim[~prim["nih"]]
        .assign(crossref_nih=contaminated[~prim["nih"]])
        .groupby("year")
        .agg(n=("pmid", "size"), crossref_nih=("crossref_nih", "sum"))
        .assign(share=lambda t: t["crossref_nih"] / t["n"])
    )
    contamination.reset_index().to_csv(OUT / "control_contamination.csv", index=False)

    # Unmatched NIH cohort, cohort construction, persisting embargoes
    full = load("nih_cohort.csv")
    full["doaj"] = full["is_in_doaj"].fillna(False).astype(bool)
    cols = [f"{k}_{h}d" for k in ("rel", "live") for h in (0, 7, 30, 120)]
    for y, row in (full.groupby("year")[cols].mean() * 100).iterrows():
        for c, v in row.items():
            num[f"full_{y}_{c}"] = round(float(v), 6)
    for (y, dj), row in (full.groupby(["year", "doaj"])[cols].mean() * 100).iterrows():
        for c, v in row.items():
            num[f"full_{y}_{'doaj' if dj else 'nondoaj'}_{c}"] = round(float(v), 6)
    num["full_n"] = {str(k): int(v) for k, v in full.groupby("year").size().items()}
    emb = full[full["pmcid"].notna()]
    for (y, dj), v in (
        emb.groupby(["year", "doaj"])["status_embargo"].apply(lambda s: s.eq("yes").mean()).items()
    ):
        num[f"embargoed_{y}_{'doaj' if dj else 'nondoaj'}"] = round(float(v), 6)
    covered = emb[emb["covered"]]
    long = (covered["release"] - covered["early"]).dt.days >= 300
    num.update(
        {
            "covered_in_pmc_n": len(covered),
            "covered_long_n": int(long.sum()),
            "covered_long_share": round(float(long.mean()), 6),
            "covered_long_all_final": bool(covered.loc[long, "route"].eq("final").all()),
        }
    )
    raw = pd.read_csv(
        IN / "nih_cohort.csv", dtype={"pmid": str}
    )  # cohort as retrieved, before any outcome join
    raw["doaj"] = raw["is_in_doaj"].fillna(False).astype(bool)
    cells = raw.groupby(["journal_id", "month", "year"]).size().unstack(fill_value=0)
    both = cells[(cells[2025] > 0) & (cells[2026] > 0)]
    matched_nih = prim[prim["nih"]]
    num["funnel"] = {
        "pool_2025": int((raw["year"] == 2025).sum()),
        "pool_2026": int((raw["year"] == 2026).sum()),
        "cells_both_years": len(both),
        "journals_both_years": both.index.get_level_values(0).nunique(),
        "in_cells_2025": int(both[2025].sum()),
        "in_cells_2026": int(both[2026].sum()),
        "balanced_max": int(both.min(axis=1).sum()),
        "matched_per_year": int((matched_nih["year"] == 2025).sum()),
    }
    num["funnel"]["lost_for_candidates"] = num["funnel"]["balanced_max"] - num["funnel"]["matched_per_year"]
    raw_doi = raw["doi"].astype("string").str.lower().str.replace(r"^https?://doi.org/", "", regex=True)
    coverage = raw.assign(in_crossref=raw_doi.isin(crossref)).groupby("year")["in_crossref"].mean()
    num["crossref_coverage"] = {str(y): round(float(v), 6) for y, v in coverage.items()}
    built = pd.read_csv(IN / "matched_articles.csv", dtype={"pmid": str})
    m26 = set(built.loc[built["nih"] & (built["year"] == 2026), "pmid"])
    f26 = raw[raw["year"] == 2026].assign(
        matched=lambda t: t["pmid"].isin(m26), in_pmc=lambda t: t["pmcid"].notna()
    )
    comp = f26.groupby("matched")[["doaj", "in_pmc"]].mean()
    num["matched_vs_unmatched_2026"] = {
        "doaj_matched": round(float(comp.loc[True, "doaj"]), 6),
        "doaj_unmatched": round(float(comp.loc[False, "doaj"]), 6),
        "in_pmc_matched": round(float(comp.loc[True, "in_pmc"]), 6),
        "in_pmc_unmatched": round(float(comp.loc[False, "in_pmc"]), 6),
    }
    nih_acc = prim[prim["nih"] & prim["accepted_date"].notna()]
    for y in (2025, 2026):
        g = nih_acc[nih_acc["year"] == y]
        num[f"accepted_{y}"] = {
            "n": len(g),
            "post_policy": int(g["covered"].sum()),
            "share": round(float(g["covered"].mean()), 6),
        }
    num["share_pubmed_date_outside_window"] = round(
        float(
            (
                (pd.to_datetime(prim["pub_date"]).dt.month > 4)
                | (pd.to_datetime(prim["pub_date"]).dt.year != prim["year"])
            ).mean()
        ),
        6,
    )

    # Europe PMC full-text received date against the go-live and release dates (check 3)
    e = pd.read_csv(IN / "europe_pmc_sample.csv", parse_dates=["fulltext_received_date"]).merge(
        rec, on="pmcid"
    )
    e = e[e["fulltext_received_date"].notna()]
    ft_live = (e["fulltext_received_date"] - e["pmc_live_date"]).dt.days
    ft_rel = (e["fulltext_received_date"] - e["pmc_release_date"]).dt.days
    discordant = (e["pmc_live_date"] - e["pmc_release_date"]).dt.days > 14
    num["epmc"] = {
        "n": len(e),
        "discordant": int(discordant.sum()),
        "closer_to_live": int((ft_live.abs() < ft_rel.abs())[discordant].sum()),
        "median_ft_minus_live": float(ft_live.median()),
    }

    # Sensitivity analyses: date rules and exposure definitions, monthly estimates
    reseed()
    pub, oa = (
        pd.to_datetime(prim["pub_date"], errors="coerce"),
        pd.to_datetime(prim["oa_pub_date"], errors="coerce"),
    )
    has = prim["pmcid"].notna()
    for h in (7, 120):
        prim[f"pubmed_{h}"] = (has & ((prim["public"] - pub).dt.days <= h)).astype(float)
        prim[f"openalex_{h}"] = (has & ((prim["public"] - oa).dt.days <= h)).astype(float)
    explicit = prim["grant_agencies"].fillna("").str.contains("NIH")
    unresolved = prim["accepted_date"].gt(prim["early"])
    sens: dict = {"unresolved_n": int(unresolved.sum())}
    checks = {
        "pubmed_date": (prim, "pubmed_{h}"),
        "openalex_date": (prim, "openalex_{h}"),
        "dates_agree_within_7d": (prim[(pub - oa).dt.days.abs().le(7)], "live_{h}d"),
        "excluding_unresolved_dates": (prim[~unresolved], "live_{h}d"),
        "drop_crossref_nih_controls": (prim[~contaminated], "live_{h}d"),
        "explicit_nih_grant_only": (prim[~prim["nih"] | explicit], "live_{h}d"),
    }
    for name, (frame, col) in checks.items():
        for h in (7, 120):
            r = did(frame, col.format(h=h))
            sens[f"{name}_{h}d"] = {k: round(float(r[k]), 6) for k in ("est", "lo", "hi")}
    for m in (1, 2, 3, 4):
        r = did(prim[prim["month"] == m], "live_120d")
        sens[f"month_{m}_120d"] = {k: round(float(r[k]), 6) for k in ("est", "lo", "hi")}
        g = prim[prim["month"] == m].groupby(["year", "nih"])["live_120d"].mean() * 100
        sens[f"month_{m}_nih_2025"], sens[f"month_{m}_nih_2026"] = (
            round(float(g[(2025, True)]), 6),
            round(float(g[(2026, True)]), 6),
        )
    us_federal = (  # US federal funders in PubMed grant fields, curated by hand
        r"(?:\bCDC\b|Centers for Disease Control|\bVA$|\(VA\)$|Department of Veterans Affairs|"
        r"Department of Defense|AHRQ HHS|OASH HHS|^National Science Foundation$|\bUSDA\b|"
        r"United States Department of Agriculture|NASA|Food and Drug Administration|FDA HHS|"
        r"U\.S\. Fish and Wildlife Service|Biomedical Advanced Research and Development|"
        r"Administration for Strategic Preparedness|\bNOAA\b|National Institutes of Health)"
    )
    agencies = prim["grant_agencies"].fillna("").str.split("; ")
    other_federal = (~prim["nih"]) & agencies.apply(lambda xs: any(pd.Series(xs).str.contains(us_federal)))
    sens["comparison_other_federal_n"] = int(other_federal.sum())
    for h in (7, 120):
        r = did(prim[~other_federal], f"live_{h}d")
        sens[f"drop_other_federal_controls_{h}d"] = {k: round(float(r[k]), 6) for k in ("est", "lo", "hi")}
    for m in (1, 2, 3, 4):  # monthly 120-day estimate by route, DOAJ share, and pre-policy contrast
        f = prim[prim["month"] == m]
        for r in ("manuscript", "final"):
            by_route = f.assign(o=((f["live_120d"] == 1) & (f["route"] == r)).astype(float))
            sens[f"month_{m}_120d_{r}"] = round(point_did(by_route, "o"), 6)
        nih25 = f[f["nih"] & (f["year"] == 2025)]
        sens[f"month_{m}_doaj_share"] = round(float(nih25["is_in_doaj"].mean()), 6)
        pre = point_did(plac[plac["month"] == m], "live_120d", (2024, 2025))
        sens[f"month_{m}_placebo_120d"] = round(pre, 6)
    (OUT / "sensitivity.json").write_text(json.dumps(sens, indent=1) + "\n")

    # October re-match: the same analysis on a later PubMed retrieval
    reseed()
    octs = load("matched_articles_oct.csv")
    oct_num: dict = {}
    oct_horizon, _ = analyse(octs, plac, rec, oct_num)
    oct_horizon.to_csv(OUT / "oct_did_by_horizon.csv", index=False)
    num["oct_month_120d"] = {
        str(m): round(point_did(octs[octs["month"] == m], "live_120d"), 6) for m in (1, 2, 3, 4)
    }
    oct_nih26 = octs[octs["nih"] & (octs["year"] == 2026)]
    sep_nih26 = set(prim.loc[prim["nih"] & (prim["year"] == 2026), "pmid"])
    num["oct_overlap_2026_nih"] = round(float(oct_nih26["pmid"].isin(sep_nih26).mean()), 6)

    # Fig. 1 curves
    reseed()
    curves_with_bands(prim.reset_index(drop=True), CURVE_DAYS).to_csv(OUT / "curves.csv", index=False)

    # Articles on terminated NIH awards (exploratory)
    reseed()
    num["terminated_awards"] = terminated_award_articles(prim, plac)

    (OUT / "numbers.json").write_text(json.dumps(num, indent=1) + "\n")
    print(
        by_horizon[by_horizon["stratum"] == "all"][["horizon", "outcome", "est", "lo", "hi"]]
        .round(1)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
