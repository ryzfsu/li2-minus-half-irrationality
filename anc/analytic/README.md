# Analytic certificate data

Run the checker from the manuscript directory:

```sh
python3 anc/check_analytic.py --tex main.tex --out replay/analytic
```

The checker uses `python-flint` 0.9.0. It reads the rational measure from
`reference.json`, verifies that it agrees with the printed tables, recomputes
the energy and five-support potential bounds, and checks the endpoint disk
and smoothing constants. It also compares the printed run-length depth
sequences with the accepted Arb partition and verifies exact support
coverage.

`depth_sequences.tex` supplies the depth table included in the manuscript.
`leaves.json.gz` records every accepted rational interval and its outward
Arb upper bound. Circle intervals use the coordinate $t=\theta/\pi$ on
$[0,1]$; exact conjugation symmetry covers the lower semicircle.

These calculations establish the finite enclosures on the full supports.
The concavity, potential modulus, smoothing, and passage to all matrix sizes
are proved in Section 4 of the manuscript.
