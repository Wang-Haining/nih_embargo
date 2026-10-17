# Provenance

## Sources

| Source | What was taken | Retrieved / snapshot |
|---|---|---|
| PubMed (E-utilities) | NIH-supported journal articles: publication types `Research Support, N.I.H., Extramural` or `Research Support, N.I.H., Intramural` with `Journal Article`, published 1 January to 30 April of 2024, 2025 and 2026; publication dates, publication types and grant agencies for these and for comparison candidates | 12 September 2026; repeated 16 October 2026 for the October re-match |
| PubMed (E-utilities) | Accepted dates (history status "accepted") for NIH-supported author manuscripts in the matched sample | 16 October 2026 |
| OpenAlex works | Journal identity, DOAJ membership, publication dates, the comparison-candidate pool (same journal, year and month, with a PMID), funder records | Snapshot of 26 June 2026 |
| PMC identifier list | PMCID for each PMID | 12 September 2026; 16 October 2026 for the October re-match |
| PMC OAI-PMH (`pmc_fm`) | Go-live date (pmc-live event), release date (pmc-release event), NIH Manuscript Submission date, accepted date, author-manuscript, live and embargo flags, licence | 16 and 17 October 2026 |
| Crossref REST API | Works with the NIH funder identifier 10.13039/100000002 issued January to April 2024 to 2026 | 12 September 2026 |
| Europe PMC REST API | Full-text received date for a seeded random sample of 200 matched 2026 articles | 16 October 2026 |
| PMC journal list | Deposit agreement of each journal ("All articles", "NIH funded articles", none) | 16 October 2026 |
| Crossref REST API (`/prefixes`) | Publisher that registered each DOI prefix of the matched articles | 17 October 2026 |
| NIH list of Method D publishers | Publishers that volunteered to start NIH Manuscript Submission deposits for authors, as archived on 15 June 2025 (https://web.archive.org/web/20250615180042/https://sharing.nih.gov/public-access-policy/submitting-to-pubmed-central/method-d-publishers) | 17 October 2026 |
| PubMed (E-utilities, efetch) | Grant numbers and agencies listed for the NIH-supported articles of the 2024, 2025 and 2026 cohorts | 17 October 2026 |
| Grant Witness (https://grantwitness.org/nih) | NIH awards terminated since 2025, including those later reinstated (award key, earliest disruption, latest restoration and status only); data last updated 13 October 2026 | 17 October 2026 |
| PMC E-utilities (esearch, `[pmcrdat]`) | Whether PMC's Live Date equals the OAI-PMH go-live date or the release date, for a seeded sample of 100 records | 17 October 2026 |

Publication date throughout: the earlier of the PubMed date and the OpenAlex date. Sampling seed 20260912. The 7-day horizon and the 1:1 comparison ratio were fixed before any 2026 outcome was retrieved.

## `data/` (inputs to `src/estimate.py`)

| File | Rows | Content | Built by |
|---|---|---|---|
| `matched_articles.csv` | 7,520 | 1,880 NIH-supported and 1,880 comparison articles in each of 2025 and 2026, matched within journal-month cells | `src/acquire/build_cohort.py` |
| `placebo_articles.csv` | 19,028 | Pre-policy sample (2024 vs 2025), same construction | `src/acquire/build_cohort.py` |
| `matched_articles_oct.csv` | 8,952 | The same match on the PubMed retrieval of 16 October 2026 | `src/acquire/build_cohort.py` (`COHORT_SUFFIX=_oct`) |
| `nih_cohort.csv` | 13,888 | All NIH-supported articles published January to April 2025 and 2026, matched or not | `src/acquire/build_cohort.py` |
| `funder_replication_articles.csv` | 2,455 | Articles with NIH support in Crossref and OpenAlex funder metadata, in the primary journal-month cells | `src/acquire/funder_cohort.py` |
| `pmc_records.csv` | 24,119 | One row per PMCID: go-live date, release date, NIH Manuscript Submission date, accepted date, flags; `record_available` is false for the 2 records the service did not return | `src/acquire/pmc_oai.py` |
| `pubmed_accepted.csv` | 1,090 | PubMed accepted dates for NIH-supported author manuscripts | `src/acquire/pubmed_accepted.py` |
| `europe_pmc_sample.csv` | 200 | Europe PMC full-text received date | `src/acquire/europe_pmc.py` |
| `journal_agreements.csv` | 842 | PMC deposit agreement by OpenAlex journal | `src/acquire/journal_agreements.py` |
| `doi_publishers.csv` | 67 | Publisher of each DOI prefix in the matched sample (Crossref member) and whether it was on NIH's Method D list | `src/acquire/doi_publishers.py` |
| `pubmed_grants.csv` | 34,241 | Grant numbers listed in PubMed for 8,186 NIH-supported articles | `src/acquire/pubmed_grants.py` |
| `grant_terminations.csv` | 2,850 | NIH parent awards terminated since 2025 (Grant Witness): key, earliest disruption, latest restoration, status | `src/acquire/grant_terminations.py` |
| `pmc_live_check.csv` | 100 | PMC Live Date compared with the go-live and release dates, seeded sample | `src/acquire/pmc_live_check.py` |
| `crossref_nih_dois.csv` | 74,741 | DOIs Crossref returns for the NIH funder identifier in the study windows | `src/acquire/crossref_funders.py` |
| `crossref_openalex_links.csv` | 74,741 | Each Crossref NIH DOI linked by DOI to the OpenAlex snapshot: publication date and year, PMID, journal, DOAJ status, and whether OpenAlex lists NIH funding | DOI join against the OpenAlex works snapshot |
| `crossref_date_audit.csv` | 3 | Crossref records retrieved, linked to OpenAlex and passing the year check, by year | `src/acquire/crossref_funders.py` |
| `manifest.json` | | Query, counts, windows, snapshot dates, denominator audit, October retrieval counts | acquisition scripts |

## `results/` (outputs of `src/estimate.py`)

| File | Content |
|---|---|
| `did_by_horizon.csv` | Difference-in-differences at 0, 7, 30, 60, 90 and 120 days, dated by going live and by the release date, for all journals, DOAJ strata and deposit-agreement strata, with 95% intervals and group means |
| `route_decomposition.csv` | The same estimate split by route (author manuscript, final published version) |
| `oct_did_by_horizon.csv` | Difference-in-differences by horizon for the October re-match |
| `curves.csv` | Share live in PMC by day 0 to 120, by stratum, year and group, with journal-cluster bootstrap bands (Fig. 1A, 1B) |
| `numbers.json` | Checks of the go-live date, release-date gaps, timing by route, pre-policy contrasts, leave-one-journal-out and capped weights, funder-metadata replication, cohort construction, unmatched-cohort availability, persisting embargoes, articles on terminated NIH awards |
| `sensitivity.json` | Date-rule and exposure sensitivity analyses and monthly estimates |
| `control_contamination.csv` | Matched comparison articles carrying NIH funding in Crossref, by year |
