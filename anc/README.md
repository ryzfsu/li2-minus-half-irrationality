# Ancillary files

This directory contains the data and scripts accompanying the manuscript.

Only Lemma 4.2 uses interval arithmetic: the five-support potential bounds,
the quadratic-energy enclosure, and the endpoint-disk bounds. The other
finite computations are exact rational or modular checks of data printed in
the paper. The scripts do not replace the written uniform arguments.

## Requirements

Python 3.12. The interval checker requires `python-flint` 0.9.0; the other
checkers use the standard library.

```sh
python3 -m venv /tmp/li2-paper-venv
/tmp/li2-paper-venv/bin/python -m pip install python-flint==0.9.0
```

The included PDF was built with pdfTeX 1.40.25 (TeX Live 2023/Debian), in
three passes; the third pass confirmed stable cross-references.

## Checks

Run from the directory containing `main.tex`:

```sh
mkdir -p replay/analytic
python3 anc/check_determinant.py --tex main.tex --output replay/determinant.json
python3 anc/check_arith.py --tex main.tex --output replay/arithmetic.json
python3 anc/check_mod11.py --tex main.tex --json replay/mod11.json
python3 anc/check_analytic.py --tex main.tex --out replay/analytic
```

The checkers read the printed tables. `check_analytic.py` also compares the
printed rational measure with `anc/analytic/reference.json`, rechecks the
reference potentials, and verifies that the printed depth sequences
reconstruct all five subdivisions. Its full leafwise Arb upper bounds are
stored in `anc/verification/analytic/leaves.json.gz`.

The included `main.pdf` matches `main.tex`. The interval checker regenerates
the run-length depth table and leafwise bounds from the included source data;
the other checkers verify the finite exact data printed in the manuscript.
Machine-readable results for all four checkers are included in
`anc/verification/`; each report records the source hash it checked.

## Scope

- `check_determinant.py` checks the printed polynomial for $G_4$, the Laurent
  moments, and the small determinant identities.
- `check_arith.py` checks Table A.1, its cut lines and integrals, the rank-one
  correction, the prime-price constants, and the formal dyadic exponent.
- `check_mod11.py` checks the displayed $9\times9$ matrix, its determinant,
  the polynomial $g(T)$, and the printed residues.
- `check_analytic.py` checks Lemma 4.2, Tables B.1--B.2, the energy and
  endpoint enclosures, the depth-sequence table, and the constants in
  Lemma 4.4.

The passage from these finite data to all required matrix sizes is proved in
the manuscript.
