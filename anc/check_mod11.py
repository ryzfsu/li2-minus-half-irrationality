#!/usr/bin/env python3
"""Recompute G_4 modulo 11 and check the printed matrix, polynomial and table.

This is an exact finite computation.  The uniform congruence and valuation
argument are proved in Section 5; this program does not infer them from tests.
Only the Python standard library is needed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import sys
import time

from tex_tables import number, require_literal, table

P = 11


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_tex(path: Path, active: tuple[Path, ...] = ()) -> str:
    path = path.resolve()
    require(path not in active, f"cyclic TeX input: {path}")
    text = path.read_text()
    text = re.sub(r"(?<!\\)%[^\n]*", "", text)

    def include(match: re.Match[str]) -> str:
        child = path.parent / match.group(1)
        if not child.suffix:
            child = child.with_suffix(".tex")
        return read_tex(child, active + (path,))

    return re.sub(r"\\(?:input|include)\{([^}]+)\}", include, text)


def printed_polynomial(tex: str) -> list[int]:
    matches = re.findall(r"g\(T\):=(.*?)\\pmod\s*\{11\}", tex, re.S)
    require(len(matches) == 1, "expected one printed definition of g(T) modulo 11")
    expression = re.sub(r"\s+", "", matches[0])
    terms = re.findall(r"[+-]?[^+-]+", expression)
    require("".join(terms) == expression, "cannot parse printed polynomial g")
    coefficients: dict[int, int] = {}
    for term in terms:
        match = re.fullmatch(r"([+-]?)(\d*)(T(?:\^(\d+|\{\d+\}))?)?", term)
        require(match is not None, f"unsupported printed term of g: {term}")
        sign, digits, variable, exponent = match.groups()
        require(bool(digits or variable), f"empty term of g: {term}")
        coefficient = int(digits) if digits else 1
        if sign == "-":
            coefficient = -coefficient
        degree = int(exponent.strip("{}")) if exponent else int(bool(variable))
        coefficients[degree] = coefficients.get(degree, 0) + coefficient
    require(max(coefficients) == 4, "printed g must have degree four")
    return [coefficients.get(i, 0) % P for i in range(5)]


def printed_table(tex: str) -> tuple[list[int], list[int]]:
    printed = table(tex, "tab:g-mod-eleven")
    require(len(printed) == 2, "expected two printed modulo-11 rows")
    rows = []
    for row, label in zip(printed, ("$w$", "$g(w)$")):
        require(row[0] == label, f"missing table row {label}")
        parsed = [number(cell) for cell in row[1:]]
        require(all(value.denominator == 1 for value in parsed), "nonintegral residue table")
        rows.append([int(value) for value in parsed])
    require(rows[0] == list(range(P)), "printed table of g does not cover F_11")
    require(len(rows[1]) == P, "printed table of g has the wrong number of values")
    return rows[0], rows[1]


def determinant(matrix: list[list[int]]) -> int:
    matrix = [[x % P for x in row] for row in matrix]
    result = 1
    for i in range(len(matrix)):
        pivot = next((j for j in range(i, len(matrix)) if matrix[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
            result = -result
        diagonal = matrix[i][i]
        result = result * diagonal % P
        inverse = pow(diagonal, -1, P)
        for j in range(i + 1, len(matrix)):
            multiplier = matrix[j][i] * inverse % P
            for k in range(i + 1, len(matrix)):
                matrix[j][k] = (matrix[j][k] - multiplier * matrix[i][k]) % P
            matrix[j][i] = 0
    return result % P


def literal_matrix(value: int) -> list[list[int]]:
    a = -2
    t = [0] + [sum(pow(a, j-u, P) * pow(u, -1, P)
                   for u in range(1, j+1)) % P for j in range(1, 10)]
    b = [0] + [sum(pow(a, j-u, P) * pow(u, -2, P)
                   for u in range(1, j+1)) % P for j in range(1, 10)]
    rows = []
    for i in range(2, 6):
        rows.append([0 if i == j else
                     (pow(a, i, P)-pow(a, j, P))*pow(j-i, -1, P) % P
                     for j in range(1, 10)])
    for i in range(2, 6):
        rows.append([(b[i]-value*pow(a, i, P)) % P if i == j else
                     (t[i]-t[j])*pow(j-i, -1, P) % P
                     for j in range(1, 10)])
    rows.append([pow(j+1, -1, P) for j in range(1, 10)])
    return rows


def literal_determinant(value: int) -> int:
    return determinant(literal_matrix(value))


def printed_matrix(tex: str) -> list[list[tuple[int, int]]]:
    """Entries of the printed matrix as (constant, coefficient of T) mod 11."""
    match = re.search(r"\\label\{eq:g4-mod-matrix\}.*?\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", tex, re.S)
    require(match is not None, "missing printed matrix eq:g4-mod-matrix")
    rows = [row for row in re.split(r"\\\\", match[1]) if row.strip()]
    require(len(rows) == 9, "printed matrix must have nine rows")
    result = []
    for row in rows:
        cells = [re.sub(r"\s+", "", cell) for cell in row.split("&")]
        require(len(cells) == 9, "printed matrix must have nine columns")
        parsed = []
        for cell in cells:
            entry = re.fullmatch(r"(\d+)(?:\+(\d*)T)?", cell)
            require(entry is not None, f"unsupported matrix entry: {cell}")
            constant = int(entry[1])
            coefficient = 0 if entry[0] == entry[1] else int(entry[2] or "1")
            require(0 <= constant < P and 0 <= coefficient < P, f"entry not reduced: {cell}")
            parsed.append((constant, coefficient))
        result.append(parsed)
    return result


def main() -> None:
    start = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tex", nargs="?", type=Path)
    parser.add_argument("--tex", dest="tex_option", type=Path)
    parser.add_argument("--json", type=Path, help="write a compact exact result")
    args = parser.parse_args()
    require(not (args.tex and args.tex_option), "give the TeX path only once")
    path = args.tex_option or args.tex or Path(__file__).resolve().parents[1]/"main.tex"
    tex = read_tex(path)
    require_literal(tex)
    coefficients = printed_polynomial(tex)
    arguments, values = printed_table(tex)
    determinants = [literal_determinant(w) for w in arguments]
    evaluations = [sum(c*pow(w, j, P) for j, c in enumerate(coefficients)) % P
                   for w in arguments]
    require(determinants == evaluations,
            "printed polynomial differs from the determinant modulo 11")
    require(values == evaluations, "printed table of g is incorrect")
    require(all(values), "g has a root in F_11")
    matrix = printed_matrix(tex)
    for w in range(P):
        literal = literal_matrix(w)
        printed = [[(c + d * w) % P for c, d in row] for row in matrix]
        require(printed == literal, f"printed matrix differs from the reduction at T={w}")
        require(determinant(printed) == evaluations[w], "printed matrix determinant")
    multiplier = re.findall(r"(\d+)\\equiv(\d+)\\pmod\{11\}", re.sub(r"\s+", "", tex))
    require(len(multiplier) == 1, "missing printed residue of the G_4 multiplier")
    require(int(multiplier[0][0]) % P == int(multiplier[0][1]) and int(multiplier[0][1]) != 0,
            "printed residue of the G_4 multiplier")
    # The determinant has degree at most four (only four rows involve T).
    # Agreement at eleven arguments proves equality of these polynomials.
    result = {
        "ok": True,
        "scope": "exact computation modulo 11; the uniform argument is Section 5",
        "tex": str(path.resolve()),
        "expanded_tex_sha256": hashlib.sha256(tex.encode()).hexdigest(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "modulus": P,
        "degree_bound": 4,
        "printed_coefficients_ascending": coefficients,
        "printed_arguments": arguments,
        "printed_values": values,
        "literal_determinants": determinants,
        "printed_matrix_checked": True,
        "root_free": True,
        "seconds": round(time.perf_counter()-start, 6),
    }
    if args.json:
        args.json.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
