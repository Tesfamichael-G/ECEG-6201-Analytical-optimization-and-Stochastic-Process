# ECEG-6201 — Analytical Optimization and Stochastic Processes

> **📖 This repository is the *source* of the course website.**
> **Read the rendered site here → <https://tesfamichael-g.github.io/ECEG-6201-Analytical-optimization-and-Stochastic-Process/>**

The rendered site is the canonical course material — it includes properly
typeset math, worked examples, code, quizzes, and solutions. The `.md` files
in `docs/` are the source those pages are built from; they don't render
cleanly on github.com because the MkDocs extensions (MathJax, admonitions,
Material icons, grid cards) are not understood by GitHub's markdown renderer.

**Looking for a correction to file?**
- Use the [feedback form](https://forms.gle/…) for anonymous or non-technical reports.
- Open a [content-correction issue](issues/new?template=content-correction.yml) if you have a GitHub account.

## Repository layout

| Path | Purpose |
|---|---|
| `docs/` | Markdown source for each page |
| `mkdocs.yml` | MkDocs configuration (theme, extensions, nav) |
| `.github/workflows/deploy.yml` | GitHub Actions workflow (builds + deploys to Pages) |
| `requirements.txt` | Python dependencies for local build |

## For maintainers

```bash
pip install -r requirements.txt
mkdocs serve          # preview at http://127.0.0.1:8000
