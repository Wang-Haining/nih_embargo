# nih_embargo

Code and data for a study of the 2025 NIH Public Access Policy and the availability of NIH-funded articles
in PubMed Central. The repository holds the article-level inputs, the scripts that built them, the
estimation script and its outputs, and the figure script. The manuscript text is not here.

## Requirements

Python 3.10 or later. With [uv](https://docs.astral.sh/uv/), `uv sync` installs the locked versions.
With pip, `pip install -r requirements.txt` in a virtual environment.

## Reproduce

```
uv run python src/estimate.py       # data/ -> results/      about one minute
uv run python src/make_figure.py    # results/ -> figures/   seconds
```

`estimate.py` reads only `data/`, uses a fixed seed and 2,000 journal-cluster bootstrap draws, and
rewrites every file in `results/`. `make_figure.py` reads only `results/` and writes `figures/fig1.pdf`
and `figures/fig1.png`. With pip, activate the environment and drop the `uv run` prefix.

## Layout

| Path | Content |
|---|---|
| `data/` | Inputs: matched samples, PMC records, accepted dates, journal agreements, publisher and funder tables |
| `results/` | Outputs of `estimate.py`: estimates by horizon and stratum, route decomposition, curves, checks |
| `figures/` | Fig. 1 (A, estimates; B, availability curves; C, deposit-route timeline) |
| `src/` | `estimate.py`, `make_figure.py`, and `acquire/` (the scripts that built `data/`) |
| `docs/PROVENANCE.md` | Every file's source, retrieval date, row count and builder; the outcome definition; how to rerun acquisition |

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

Haining Wang, Indiana University School of Medicine, hw56@iu.edu

## License

MIT; see [LICENSE](LICENSE).
