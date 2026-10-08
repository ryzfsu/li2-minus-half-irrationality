# Irrationality of the dilogarithm at -1/2

Ruiyi Zhang

DOI: [10.5281/zenodo.23247749](https://doi.org/10.5281/zenodo.23247749)

This repository contains the manuscript proving that Li_2(-1/2) is irrational, its LaTeX source, and the ancillary verification files.

- `main.pdf`: the manuscript (23 pages)
- `main.tex`: LaTeX source (needs `anc/analytic/depth_sequences.tex`)
- `anc/`: scripts and data for the checks listed in Appendix F
- `SHA256SUMS`: SHA-256 checksums of all files above

Build: `pdflatex main.tex` (twice).

Checks: see `anc/README.md`.

## Citation

Ruiyi Zhang, *Irrationality of the dilogarithm at -1/2*, preprint, Zenodo, 2026. https://doi.org/10.5281/zenodo.23247749

```bibtex
@misc{zhang2026li2minushalf,
  author    = {Zhang, Ruiyi},
  title     = {Irrationality of the dilogarithm at -1/2},
  year      = {2026},
  publisher = {Zenodo},
  note      = {Preprint},
  doi       = {10.5281/zenodo.23247749},
  url       = {https://doi.org/10.5281/zenodo.23247749}
}
```
