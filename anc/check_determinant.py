#!/usr/bin/env python3
"""Exact finite checks for Section 2; the TeX is the numerical input.

The printed G_4 identity is proved by exact interpolation using the degree
bound proved in the paper. The other finite checks are consistency tests for
the displayed, uniformly proved identities, not replacements for their proofs.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial, gcd, prod
from pathlib import Path
import platform
import re
import sys
import time

from tex_tables import compact, keyed_rows, number, require_literal


def det(matrix):
    matrix = [[F(x) for x in row] for row in matrix]
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    value = F(1)
    for k in range(size):
        pivot = next((i for i in range(k, size) if matrix[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            value = -value
        divisor = matrix[k][k]
        value *= divisor
        for i in range(k + 1, size):
            factor = matrix[i][k] / divisor
            for j in range(k + 1, size):
                matrix[i][j] -= factor * matrix[k][j]
            matrix[i][k] = F(0)
    return value


def source_matrix(n, target):
    a, m, s, k, h = -2, 4 * n, n, n, n
    rows = range(s + k, s + k + m)
    columns = range(s, s + 2 * m + h)
    # These sums implement the definitions, independently of the recurrence.
    t = [sum((F(a ** (j - u), u) for u in range(1, j + 1)), F(0))
         for j in range(max(columns) + 1)]
    b = [sum((F(a ** (j - u), u * u) for u in range(1, j + 1)), F(0))
         for j in range(max(columns) + 1)]
    top = [[F(a ** i - a ** j, j - i) if i != j else F(0)
            for j in columns] for i in rows]
    middle = [[(t[i] - t[j]) / (j - i) if i != j
               else b[i] - a ** i * target
               for j in columns] for i in rows]
    bottom = [[F(1, j + ell) for j in columns] for ell in range(1, h + 1)]
    return top + middle + bottom


def polynomial_coefficients(n):
    degree = 4 * n
    values = [det(source_matrix(n, F(t))) for t in range(degree + 1)]
    result = [F(0)] * (degree + 1)
    falling = [F(1)]
    for k in range(degree + 1):
        coefficient = values[0] / factorial(k)
        for j, term in enumerate(falling):
            result[j] += coefficient * term
        values = [values[j + 1] - values[j] for j in range(len(values) - 1)]
        following = [F(0)] * (len(falling) + 1)
        for j, term in enumerate(falling):
            following[j] -= k * term
            following[j + 1] += term
        falling = following
    return result


def evaluate(coefficients, target):
    return sum((coefficient * target ** j
                for j, coefficient in enumerate(coefficients)), F(0))


def vandermonde(points):
    return prod((points[j] - points[i] for i in range(len(points))
                 for j in range(i + 1, len(points))), start=F(1))


def determinant_factor_test(m, s, k, h):
    a, length, e = F(-2), m + h, m * (m - 1) // 2
    indices_i = range(s + k, s + k + m)
    indices_j = range(s, s + 2 * m + h)
    x = [F(i + 1, m + 2) for i in range(m)]
    v = [F(b + 1, length + 2) for b in range(length)]
    z = [F(-2) + F(7 * (d + 1), 5 * (m + 1)) for d in range(m)]
    dz = det([[a ** i * q ** (-i - 1) for q in z] for i in indices_i])
    dzv = det([[q ** j for j in indices_j] for q in z + v])
    dx = det([[q ** i for q in x] for i in indices_i])
    dxv = det([[1 / (2 + q * t) for t in v] for q in x]
              + [[t ** ell for t in v] for ell in range(h)])
    denominator = prod((2 + q * t for q in x for t in v), start=F(1))
    cross = prod((q - t for q in z for t in v), start=F(1))
    expected_dz = ((-1) ** e * a ** (m * (s + k) + e) * vandermonde(z)
                   * prod((q ** (-s - k - m) for q in z), start=F(1)))
    expected_dzv = ((-1) ** (m * length) * vandermonde(z) * vandermonde(v)
                    * prod((q ** s for q in z + v), start=F(1)) * cross)
    expected_dx = vandermonde(x) * prod((q ** (s + k) for q in x), start=F(1))
    expected_dxv = (a ** e * vandermonde(x) * vandermonde(v)
                    * prod((q ** h for q in x), start=F(1)) / denominator)
    assert (dz, dzv, dx, dxv) == (expected_dz, expected_dzv, expected_dx, expected_dxv)
    integrand = (vandermonde(x) ** 2 * vandermonde(v) ** 2 * vandermonde(z) ** 2
                 * cross * prod((q ** (s + k + h) for q in x), start=F(1))
                 * prod((q ** s for q in v), start=F(1))
                 * prod((q ** (-m - k) for q in z), start=F(1)) / denominator)
    scalar = (-1) ** (e + m * length) * a ** (m * (s + k) + 2 * e)
    assert dz * dzv * dx * dxv == scalar * integrand
    if m == 4 * s and s == k == h:
        assert scalar == F(2) ** (3 * m * m // 2 - m)


def load_tex(path, seen=None):
    if seen is None:
        seen = {}
    path = path.resolve()
    if path in seen:
        return "", seen
    data = path.read_text(encoding="utf-8")
    seen[path] = hashlib.sha256(data.encode()).hexdigest()
    pieces = [data]
    for filename in re.findall(r"\\input\{([^}]+)\}", data):
        child = path.parent / filename
        if child.suffix != ".tex":
            child = child.with_suffix(".tex")
        if child.is_file():
            piece, seen = load_tex(child, seen)
            pieces.append(piece)
    return "\n".join(pieces), seen


def run(tex_path):
    started = time.monotonic()
    tex, tex_hashes = load_tex(tex_path)
    tex = re.sub(r"(?<!\\)%[^\n]*", "", tex)
    require_literal(tex)
    denominators = re.findall(r"(?<![0-9])([0-9]+)\s*\\,\s*G_4\(T\)\s*=\s*\\sum_\{j=0\}\^4p_jT\^j", tex)
    cells = keyed_rows(tex, "tab:g4-coefficients", 2)
    parsed = [(number(j), number(cell[0])) for j, cell in cells.items()]
    assert all(j.denominator == value.denominator == 1 for j, value in parsed), "Nonintegral G_4 table entry"
    rows = [(int(j), int(value)) for j, value in parsed]
    assert len(denominators) == 1, "Expected one printed G_4 multiplier"
    assert len(rows) == 5, "Expected the five printed G_4 coefficients"
    denominator = int(denominators[0])
    printed = dict(rows)
    assert sorted(printed) == list(range(5)), "Repeated or missing table row"
    assert denominator > 0 and gcd(*printed.values()) == 1
    actual = polynomial_coefficients(1)
    assert actual == [F(printed[j], denominator) for j in range(5)], "Printed G_4 differs"
    extra_targets = [F(-3), F(1, 2), F(7)]
    for target in extra_targets:
        assert evaluate(actual, target) == det(source_matrix(1, target))
    residues = [(value.numerator * pow(value.denominator, -1, 11)) % 11 for value in actual]
    assert residues == [4, 2, 7, 3, 9]
    root_values = [int(evaluate(residues, target) % 11) for target in range(11)]
    assert 0 not in root_values

    radius_match = re.search(r"Put \$r=([0-9]+)/([0-9]+)\$", tex)
    assert radius_match, "The printed contour radius was not found"
    radius = F(int(radius_match[1]), int(radius_match[2]))
    assert radius == F(3, 5)
    # c times the sum of the two oriented arc primitives is exactly zero.
    for d in range(-16, 17):
        if d == 0:
            continue
        real_part = ((-radius) ** d - F(-2) ** d + 1 - radius ** d) / d
        upper = ((-radius) ** d - radius ** d) / d
        lower = (radius ** d - (-radius) ** d) / d
        assert upper + lower == 0
        assert real_part - upper == (1 - F(-2) ** d) / d
    for i in range(13):
        for j in range(13):
            moment = F(0) if i == j else (1 - F(-2) ** (j - i)) / (j - i)
            entry = F(0) if i == j else (F(-2) ** i - F(-2) ** j) / (j - i)
            assert F(-2) ** i * moment == entry
    t, b = F(0), F(0)
    for j in range(1, 13):
        t, b = -2 * t + F(1, j), -2 * b + F(1, j * j)
        assert t == sum((F((-2) ** (j - u), u) for u in range(1, j + 1)), F(0))
        assert b == sum((F((-2) ** (j - u), u * u) for u in range(1, j + 1)), F(0))

    test_parameters = [(m, s, k, h) for m in (1, 2, 3)
                       for s, k in ((0, 2), (1, 0)) for h in (0, 1, 2)]
    test_parameters.append((4, 1, 1, 1))
    for parameters in test_parameters:
        determinant_factor_test(*parameters)
    return {
        "ok": True,
        "scope": "Exact G4 identity; finite rational tests of uniformly proved identities",
        "tex_sha256": {str(path): digest for path, digest in tex_hashes.items()},
        "G4_multiplier": denominator,
        "G4_integer_coefficients": [printed[j] for j in range(5)],
        "G4_rational_coefficients": [str(value) for value in actual],
        "G4_mod11": residues,
        "G4_values_mod11": root_values,
        "extra_exact_targets": [str(target) for target in extra_targets],
        "contour_radius_from_tex": str(radius),
        "nonzero_laurent_moment_range": [-16, 16],
        "zero_laurent_moment": "Retained as zero; the uniform path calculation is in the paper",
        "moment_recurrence_test_range": [1, 12],
        "four_determinant_factor_tests": len(test_parameters),
        "seconds": round(time.monotonic() - started, 6),
        "python": sys.version,
        "platform": platform.platform(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, help="Section 2 or a main TeX file including it")
    parser.add_argument("--output", type=Path, help="Optional JSON result path")
    args = parser.parse_args()
    base = Path(__file__).resolve().parent.parent
    tex = args.tex
    if tex is None:
        candidate = base / "parts" / "sec2.tex"
        tex = candidate if candidate.is_file() else base / "main.tex"
    result = run(tex)
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
