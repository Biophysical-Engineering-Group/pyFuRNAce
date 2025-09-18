# Contributing to pyFuRNAce

Thanks for your interest in improving **pyFuRNAce**!

This guide covers setup, coding style, tests, docs, and how to submit changes.

> Verified project facts used here: Python >= **3.10+**, Streamlit GUI (`pyfurnace` / `python -m pyfurnace`), optional **OAT** (oxDNA analysis tools) and **Perl** (for sequence generation), docs on **Read the Docs**, `.pre-commit-config.yaml`, `tests/`, and CI workflows. License: **GPL-3.0**.
> Sources: GitHub repo & README; docs site; Actions overview.

---

## Code of Conduct
Be kind, curious, and constructive. Harassment or disrespectful behavior isn’t tolerated. If something makes you uncomfortable, contact the maintainers privately via GitHub.

---

## What contributions are welcome?
- Bug reports & minimal reproductions
- Small fixes (typos, docstrings)
- New features or UI polish (Streamlit app)
- Performance / reliability improvements
- Tests (unit/functional) and coverage bumps
- Docs & examples (tutorials/how-tos)
- New motifs additions

Before starting a **large** feature, open an Issue to discuss scope and approach.

---

## Local development setup

### Prerequisites
- **Python**: 3.10+ (Linux/macOS/Windows)
- **Git**
- (Optional) **OAT** (oxDNA analysis) for PDB↔oxDNA conversion & force files
- (Optional) **Perl** (used by Revolvr for sequence generation)

### Clone & create a virtual environment (suggested via conda)
```bash
git clone https://github.com/Biophysical-Engineering-Group/pyFuRNAce.git
conda create -n pyfurnace_dev python=3.10
conda activate pyfurnace_dev
```

### Install Perl for sequence generation (if needed)
```bash
conda install -c bioconda perl
```

### Install pyFuRNAce and dependencies
```bash
cd pyFuRNAce
python -m pip install -e ".[dev]"
```

### OAT (oxDNA analysis tools)
```bash
python -m pip install "git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis"
```


### Running the app locally
```bash
pyfurnace
# or
python -m pyfurnace
```
This launches the Streamlit GUI in your browser.

### Tests & coverage
The repository includes tests/ and a .coveragerc.

```bash
pytest
# with coverage (if pytest-cov installed)
pytest --cov=pyfurnace --cov-report=term-missing
```
Please keep (or improve) coverage if you change behavior.

### Linting, formatting & pre-commit
This repo uses pre-commit. Install hooks once:

```bash
pip install pre-commit
pre-commit install
```

### Run checks:

```bash
pre-commit run --all-files
```
Hooks (formatters/linters) are defined in .pre-commit-config.yaml. Ensure your changes pass locally before opening a PR.

---

## Style & typing

Let the formatter decide (no manual bikeshedding).

Prefer type hints (typing), docstrings, and clear APIs

Keep functions focused; write tests for edge cases.

---

## Documentation

Docs live in docs/ and are published on Read the Docs.

Typical Sphinx workflow:

```bash
pip install -r docs/requirements.txt || pip install sphinx myst-parser sphinx-autobuild
sphinx-build -b html docs docs/_build/html
# or auto-reload while writing:
sphinx-autobuild docs docs/_build/html
```
Add examples under examples/ where helpful.

---

## Branching, commits & PRs

Commits: clear messages; Conventional Commits encouraged (fix: …, feat: …, docs: …, test: …, chore: …)

Update tests and docs with the code

PR checklist

 - [ ] pre-commit run --all-files passes

 - [ ] pytest passes; coverage not reduced

 - [ ] Docs updated & build clean locally

 - [ ] Screenshots/GIFs for UI changes

---

## Notes on any new dependencies

CI (GitHub Actions) will run on your PR; please address failures before requesting review.

---

## Licensing

By contributing, you agree that your contributions are licensed under GPL-3.0 (see LICENSE).

---

## Getting help

Search existing Issues/PRs/Discussions or get in touch with maintainers.

Happy hacking! 🧬
