"""Read the manuscript's literal TeX tables without executing TeX macros.

This module only locates printed cells and parses their explicit numerals.
It contains no replacement numerical dataset. Unsupported cell syntax,
duplicate rows, and missing tables are errors rather than silent fallbacks.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import re


MACHINE_MACROS = (
    "AnalyticValue", "AnalyticCapRow", "GFourInteger", "GFourCoefficient",
    "CheckTime", "ArithRow", "ArithAffine", "ArithTotals",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_tex(path: Path, active=()) -> str:
    path = path.resolve()
    require(path not in active, f"cyclic TeX input: {path}")
    source = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))

    def include(match):
        child = path.parent / match.group(1)
        if not child.suffix:
            child = child.with_suffix(".tex")
        return read_tex(child, active + (path,))

    return re.sub(r"\\(?:input|include)\s*\{([^}]+)\}", include, source)


def require_literal(source):
    for name in MACHINE_MACROS:
        require(not re.search(r"\\" + name + r"(?![A-Za-z])", source),
                f"machine numerical macro remains: {name}")


def compact(cell):
    return re.sub(r"\s+", "", re.sub(r"\\[,;! ]", "", cell.replace("$", "")))


def number(cell):
    cell = compact(cell)
    fraction = re.fullmatch(r"\\(?:t?frac)\{(-?\d+)\}\{(\d+)\}", cell)
    if fraction:
        return Fraction(int(fraction[1]), int(fraction[2]))
    require(re.fullmatch(r"[+-]?\d+(?:\.\d+|/\d+)?", cell) is not None,
            f"not a literal rational cell: {cell}")
    return Fraction(cell)


def numbers(cell, count):
    value = compact(cell).replace(r"\,", "").replace(r"\ ", "")
    if value.startswith("(") and value.endswith(")"):
        value = value[1:-1]
    parts = value.split(",")
    require(len(parts) == count, f"expected {count} literal values: {cell}")
    return [number(part) for part in parts]


def upper_bound(cell):
    value = compact(cell)
    require(value.startswith("<"), f"expected a strict printed upper bound: {cell}")
    return number(value[1:])


def table(source, label):
    pattern = r"\\begin\{(table|longtable|minipage)\}(.*?)\\end\{\1\}"
    matches = [match for match in re.finditer(pattern, source, re.S)
               if r"\label{" + label + "}" in match[2]]
    require(len(matches) == 1, f"expected one printed table {label}")
    environment, body = matches[0][1], matches[0][2]
    if environment == "minipage":
        require(r"\captionof{table}" in body, f"uncaptioned table minipage: {label}")
    if environment != "longtable":
        begin = re.search(r"\\begin\{tabular\}(?:\[[^]]*\])?", body)
        require(begin is not None, f"missing tabular in {label}")
        body = body[begin.end():]
        body = body.split(r"\end{tabular}", 1)[0]
    # Discard the column-format argument, including its nested braces.
    start = body.find("{")
    require(start >= 0, f"missing column specification in {label}")
    depth, end = 0, None
    for i in range(start, len(body)):
        if body[i] == "{":
            depth += 1
        elif body[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    require(end is not None, f"unterminated columns in {label}")
    body = body[end:]
    rows = []
    for row in re.split(r"\\\\(?:\[[^]]*\])?", body):
        row = re.sub(r"\\(?:hline|endfirsthead|endhead|endfoot|endlastfoot)\b", "", row).strip()
        if "&" in row:
            rows.append([cell.strip() for cell in row.split("&")])
    return rows


def keyed_rows(source, label, columns):
    rows = table(source, label)
    require(rows, f"empty table {label}")
    result = {}
    for row in rows[1:]:
        require(len(row) == columns, f"wrong column count in {label}: {row}")
        key = compact(row[0])
        require(key not in result, f"duplicate row {key} in {label}")
        result[key] = row[1:]
    return result


def affine(cell):
    """Parse a literal affine expression c+b*rho, terms in either order."""
    text = compact(cell)
    require(text != "", "empty affine expression")
    terms = re.findall(r"[+-]?[^+-]+", text)
    require("".join(terms) == text, f"cannot parse affine expression: {cell}")
    intercept, slope = Fraction(0), Fraction(0)
    for term in terms:
        match = re.fullmatch(r"([+-]?)(\d+(?:/\d+)?)?(\\rho)?", term)
        require(match is not None and (match[2] or match[3]), f"unsupported affine term: {term}")
        sign = -1 if match[1] == "-" else 1
        coefficient = number(match[2]) if match[2] else Fraction(1)
        if match[3]:
            slope += sign * coefficient
        else:
            intercept += sign * coefficient
    return intercept, slope


def arithmetic_rows(source):
    """Rows (lo, hi, lines, values, intercept, slope, integral) of Table A.1."""
    result = []
    for row in table(source, "tab:arithmetic"):
        interval = compact(row[0])
        if not interval.startswith("("):
            continue
        require(len(row) == 5, "arithmetic table must have five columns")
        lo, hi = numbers(row[0], 2)
        lines = [affine(part) for part in compact(row[1]).split(",")]
        values = [number(part) for part in compact(row[2]).split(",")]
        require(all(v.denominator == 1 for v in values), "nonintegral rank value")
        intercept, slope = affine(row[3])
        result.append((lo, hi, lines, [int(v) for v in values], intercept, slope, number(row[4])))
    require(len(result) == 44, "expected exactly 44 literal arithmetic rows")
    return result


def analytic_data(source):
    require_literal(source)
    rows = {}
    labels = {
        r"x\in[0,1]": "x", r"v\in[0,1]": "v",
        r"z\in[-2,-3/5]": "zn", r"z\in[3/5,1]": "zp",
        r"|z|=3/5,\operatorname{Im}z\ge0": "zc",
    }
    for row in table(source, "tab:analytic-bounds")[1:]:
        require(len(row) == 5, "analytic support table must have five columns")
        label = compact(row[0])
        require(label in labels, f"unknown printed support: {row[0]}")
        key = labels[label]
        require(key not in rows, f"duplicate support row: {key}")
        count, depth, upper, bound = map(number, row[1:])
        require(count.denominator == depth.denominator == 1, "nonintegral partition count/depth")
        rows[key] = dict(label=row[0], count=int(count), depth=int(depth), upper=upper, cap=bound)
    require(set(rows) == {"x", "v", "zn", "zp", "zc"}, "five printed support rows required")
    cells = keyed_rows(source, "tab:analytic-constants", 2)
    values = {}

    def take(label, keys, *, upper=False):
        key = compact(label)
        require(key in cells, f"missing analytic constants row: {label}")
        cell = cells.pop(key)[0]
        parsed = [upper_bound(cell)] if upper else numbers(cell, len(keys))
        for name, value in zip(keys, parsed):
            values[name] = str(value)

    take(r"intervals on $x$, $v$, $z$; Poisson kernels", ("n_x", "n_v", "n_zreal", "n_poisson"))
    take(r"(K_x,K_v,K_z)", ("Kx", "Kv", "Kz"))
    take(r"\delta", ("delta",))
    take(r"\sup_{|q|\le\delta}F_x(q)", ("endpoint_x",), upper=True)
    take(r"\sup_{|q|\le\delta}F_v(q)", ("endpoint_v",), upper=True)
    take(r"(H_x,H_v,H_{z,\mathrm{real}},H_{z,\mathrm{circle}})", ("Mx", "Mv", "Mzreal", "Mz_circle_ceiling"))
    take(r"subintervals in Table~\ref{tab:analytic-bounds}", ("leaves",))
    # Accept either combined or separate endpoint precision rows.
    precision_label = ("precision (decimal digits)"
                       if compact("precision (decimal digits)") in cells
                       else ("endpoint-disk precision(bits)"
                             if compact("endpoint-disk precision (bits)") in cells
                             else "precision"))
    take(precision_label, ("precision",))
    # Consume a separate endpoint-disk precision row when present.
    if compact("endpoint-disk precision (bits)") in cells:
        take("endpoint-disk precision (bits)", ("endpoint_precision_bits",))
    take("series terms", ("series_N",))
    take(r"H_z", ("Mz",))
    take(r"(C_x,C_v,C_z)", ("Cx", "Cv", "Cz"))
    take(r"C_x+\tfrac54C_v+C_z", ("weighted_density_coeff",))
    take(r"(L_D,L_z)", ("LD", "Lz"))
    take(r"A_0", ("A0",))
    take(r"\tfrac{17}{16}\cdot8\epsilon/\delta^2+A_0$ at $m=64", ("allowance",), upper=True)
    take(r"3(C_x+\tfrac54C_v+C_z)+200", ("delta_coeff",))
    take(r"C_*", ("Cstar_upper",), upper=True)
    remainder = compact(cells.pop(compact("series remainder per pair"))[0])
    match = re.fullmatch(r"<([0-9.]+)\\times10\^\{(-?\d+)\}", remainder)
    require(match is not None, "unparseable printed series remainder")
    values["series_tail_mantissa"], values["series_tail_power"] = match.groups()
    require(not cells, f"unchecked analytic constant rows: {sorted(cells)}")

    # Constants printed in the text of Lemmas 4.3 and 4.4, read literally.
    flat = compact(source)
    def single(pattern, what):
        found = re.findall(pattern, flat)
        require(len(found) >= 1 and len(set(found)) == 1, f"printed {what} not found uniquely")
        return found[0]
    values["Kmod"] = single(r"\|\\mathcalU_i\(q\)-\\mathcalU_i\(q'\)\|\\le(\d+)H_id\(1\+\|\\logd\|\)", "modulus constant")
    values["m0"] = single(r"Ifm=4n\\ge(\d+),then", "matrix size threshold")
    rate, coeff = single(r"\\log\|G_m\(\\beta\)\|\\le-(\d+\.\d+)m\^2\+(\d+)m\\logm\.", "final bound")
    values["rate"], values["error_coeff"] = rate, coeff
    values["B"] = single(r"=24566\\epsilon\\omega<(\d+)\\epsilon\\omega", "smoothing coefficient")
    values["self_coeff"] = single(r"\\log\\epsilon=\\frac\{(\d+)\}\{2\}m\\logm", "diagonal coefficient")
    values["self_coeff"] = str(Fraction(int(values["self_coeff"]), 2))

    numeric = r"(-?\d+\.\d+)"
    expressions = [
        ("qlo", "qhi", r"(?:q_-\s*=\s*)?" + numeric + r"\s*&?\s*<\s*\\mathcal Q\(\\nu\)\s*<\s*" + numeric),
        ("alpha_lower", "alpha_upper", numeric + r"\s*&?\s*<\s*\\alpha\s*<\s*" + numeric),
        ("Tlo", "Thi", numeric + r"\s*&?\s*<\s*\\varkappa\s*<\s*" + numeric),
    ]
    for lower, upper, pattern in expressions:
        matches = re.findall(pattern, source)
        require(len(matches) == 1, f"expected one printed enclosure: {lower}/{upper}")
        values[lower], values[upper] = matches[0]
    return values, rows


def depth_sequence_runs(source):
    """Read the printed run-length encoded dyadic depths by support."""
    labels = {"x": "x", "v": "v", r"z_-": "zn", r"z_+": "zp", r"z_\circ": "zc"}
    runs, chunks = {}, {}
    for row in table(source, "tab:depth-sequences"):
        if compact(row[0]) not in labels:
            continue
        require(len(row) == 3, f"wrong depth-sequence row: {row}")
        key = labels[compact(row[0])]
        chunk = int(number(row[1]))
        require(chunk == chunks.get(key, 0) + 1, (key, "nonconsecutive depth chunk"))
        chunks[key] = chunk
        cell = compact(row[2])
        match = re.fullmatch(r"\\texttt\{([0-9:,]+)\}", cell)
        require(match is not None, (key, "invalid run-length sequence"))
        parsed = []
        for token in match[1].split(","):
            pair = re.fullmatch(r"(\d+):(\d+)", token)
            require(pair is not None, (key, "invalid depth:count token", token))
            depth, count = map(int, pair.groups())
            require(depth >= 0 and count > 0, (key, "invalid run"))
            parsed.append((depth, count))
        runs.setdefault(key, []).extend(parsed)
    require(set(runs) == {"x", "v", "zn", "zp", "zc"}, "five depth-sequence supports required")
    return runs


def check_times(source):
    result = {}
    for row in table(source, "tab:checks")[1:]:
        require(len(row) == 3, "computer-check table must have three columns")
        script = re.fullmatch(r"\\texttt\{check\\_([a-z0-9]+)\.py\}", compact(row[0]))
        require(script is not None, f"unknown script cell: {row[0]}")
        key = "arithmetic" if script[1] == "arith" else script[1]
        require(key not in result, f"duplicate timing row: {key}")
        result[key] = float(upper_bound(row[2]))
    require(set(result) == {"determinant", "arithmetic", "analytic", "mod11"}, "four printed times required")
    return result


def display_rows(body):
    """Split an align/gather at top-level TeX line breaks only."""
    rows, start, depth, i = [], 0, 0, 0
    while i < len(body):
        if body.startswith(r"\\", i) and depth == 0:
            rows.append(body[start:i])
            i += 2
            match = re.match(r"\[[^]]*\]", body[i:])
            if match:
                i += match.end()
            start = i
            continue
        if body[i] == "{" and (i == 0 or body[i-1] != "\\"):
            depth += 1
        elif body[i] == "}" and (i == 0 or body[i-1] != "\\"):
            depth -= 1
        i += 1
    if body[start:].strip():
        rows.append(body[start:])
    return rows


def numbered_equations(source):
    references = set(re.findall(r"\\(?:eqref|ref)\{([^{}]+)\}", source))
    numbered = []
    for match in re.finditer(r"\\begin\{(equation|align|gather)\}(.*?)\\end\{\1\}", source, re.S):
        rows = [match[2]] if match[1] == "equation" else display_rows(match[2])
        for row in rows:
            if re.search(r"\\(?:notag|nonumber)\b", row):
                continue
            labels = re.findall(r"\\label\{([^{}]+)\}", row)
            require(len(labels) == 1, "each numbered equation needs exactly one label")
            require(labels[0] in references, f"unreferenced numbered equation: {labels[0]}")
            numbered.append(labels[0])
    return numbered
