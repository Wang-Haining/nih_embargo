# nih_embargo

Code and data for **Ending the NIH embargo accelerated public access to funded research, but substantial delays remain** (Haining Wang, 2026).

On 1 July 2025 the NIH Public Access Policy ended the 12-month embargo on funded articles. Matching each NIH-supported article to an article in the same journal and month that reported no US federal support, and dating availability by the day each article went live in PubMed Central (PMC), the policy did not detectably change availability within a week of publication (difference-in-differences 0.3 percentage points, 95% CI −1.7 to 1.8). By 120 days, availability had risen by 19.4 points more for NIH-supported articles (14.3 to 25.2). Most of the gain came through author manuscripts, which went live a median 13 days after submission; 93% were submitted after publication.

## Reproduce the results

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python src/estimate.py       # writes every reported estimate to results/
python src/make_figure.py    # writes figures/fig1.pdf and fig1.png from results/
```

`src/estimate.py` reads only `data/` and runs in about a minute on a laptop (fixed seed, 2,000 journal-cluster bootstrap draws).

## Layout

| Path | Content |
|---|---|
| `data/` | Article-level inputs: the matched sample, the pre-policy sample, the funder-defined sample, the October re-match, the full NIH-supported cohort, one PMC record per PMCID, and the supporting tables listed in [PROVENANCE.md](PROVENANCE.md) |
| `results/` | Every estimate with its confidence interval, the cumulative-availability curves with bootstrap bands, and the sensitivity analyses |
| `src/estimate.py` | Outcome dating, difference-in-differences estimator, journal-cluster bootstrap, checks of the go-live date, sensitivity analyses |
| `src/make_figure.py` | Figure 1 |
| `src/acquire/` | Scripts that built `data/` from PubMed, PMC, OpenAlex, Crossref, Europe PMC and the PMC journal list. They need network access and a local OpenAlex works snapshot (`OPENALEX_WORKS`), and are included for documentation |

## Acquisition order

The acquisition scripts are included for documentation; `data/` already holds their outputs. Run from the repository root with network access and `OPENALEX_WORKS` set: `src/acquire/build_cohort.py` (cohort, matching, PMC identifiers; intermediate files go to `src/acquire/work/`), `src/acquire/crossref_funders.py`, `src/acquire/funder_cohort.py`, then `src/acquire/pmc_oai.py` on the PMCIDs of every article table, followed by `src/acquire/pubmed_accepted.py`, `src/acquire/europe_pmc.py`, `src/acquire/journal_agreements.py`, `src/acquire/doi_publishers.py`, `src/acquire/pmc_live_check.py`, `src/acquire/pubmed_grants.py` and `src/acquire/grant_terminations.py`. With `COHORT_SUFFIX=_oct`, `build_cohort.py` builds the October re-match from a later PubMed retrieval.

## Outcome

PMC's OAI-PMH service reports, for each article, a release date and the date the record went live. For author manuscripts the two coincide; for 42 to 53% of final published versions the release date precedes the go-live date by more than a week. PMC defines its Live Date as the date an article was first available in PMC; for a seeded sample of 100 records, PMC's Live Date equals the go-live date in every case (`data/pmc_live_check.csv`). The outcome is whether an article was live in PMC within a given number of days of its publication date (the earlier of the PubMed and OpenAlex dates). Estimates using the release date are reported alongside.

## Citation

```bibtex
@misc{wang2026embargo,
  author        = {Wang, Haining},
  title         = {Ending the {NIH} embargo accelerated public access to funded research, but substantial delays remain},
  year          = {2026},
  eprint        = {2609.14842},
  archivePrefix = {arXiv},
  primaryClass  = {cs.DL},
  doi           = {10.48550/arXiv.2609.14842},
  url           = {https://arxiv.org/abs/2609.14842}
}
```

## Contact

Haining Wang, Indiana University School of Medicine, hw56@iu.edu

## License

MIT; see [LICENSE](LICENSE). Lint settings are in `pyproject.toml`.
