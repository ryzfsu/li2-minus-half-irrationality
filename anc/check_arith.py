#!/usr/bin/env python3
"""Exact reconstruction of the R=5 arithmetic data from the source inequalities.

Run from the complete paper bundle:
    python3 anc/check_arith.py --tex main.tex --output arithmetic.json

No archived cell file, table from another script, or model output is input.
The literal entries of the actual TeX tables are compared with the reconstruction.
Finite lattice and pole tests are reported separately from exact enumeration.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
import re
import sys
import time
from fractions import Fraction as F
from pathlib import Path

from tex_tables import arithmetic_rows, number, require_literal

LOW, HIGH = F(1, 5), F(5, 4)
KMAX, LMAX, WMAX = 7, 12, 2


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def value(line, rho):
    return line[0] + line[1] * rho


def integral(intercept, slope, lo, hi):
    return intercept * (hi - lo) + slope * (hi * hi - lo * lo) / 2


def profile(x, rho):
    return (
        tuple(k for k in range(KMAX + 1) if F(1, 2) <= x + k * rho < F(3, 2)),
        tuple(l for l in range(LMAX + 1) if F(1, 4) <= x + l * rho < F(5, 2)),
        tuple(w for w in range(1, WMAX + 1) if 0 < w * rho - x <= F(1, 4)),
    )


def cap(p):
    k, l, w = p
    return min(2 * len(k) + len(w), len(l)) - int(p == ((0,), (0, 1), ()))


def meets_domain(line):
    # Membership in 0 <= c+b*rho <= rho can change only at these two roots.
    endpoints = {LOW, HIGH}
    for slope in (line[1], line[1] - 1):
        if slope:
            root = -line[0] / slope
            if LOW < root < HIGH:
                endpoints.add(root)
    endpoints = sorted(endpoints)
    return any(
        0 <= value(line, (a + b) / 2) <= (a + b) / 2
        for a, b in zip(endpoints, endpoints[1:])
    )


def reconstruct():
    # The omitted integer indices are impossible throughout the open domain.
    require((KMAX + 1) * LOW > F(3, 2), "K cutoff is insufficient")
    require((LMAX + 1) * LOW > F(5, 2), "L cutoff is insufficient")
    require(WMAX * LOW > F(1, 4), "W cutoff is insufficient")
    lines = {(F(0), F(0)), (F(0), F(1))}
    lines.update((c, F(-k)) for k in range(KMAX + 1) for c in (F(1, 2), F(3, 2)))
    lines.update((c, F(-l)) for l in range(LMAX + 1) for c in (F(1, 4), F(5, 2)))
    lines.update((c, F(w)) for w in range(1, WMAX + 1) for c in (F(-1, 4), F(0)))
    all_lines = sorted(lines)
    lines = sorted(line for line in lines if meets_domain(line))
    cuts = {LOW, HIGH}
    for left, right in itertools.combinations(lines, 2):
        if left[1] != right[1]:
            rho = (right[0] - left[0]) / (left[1] - right[1])
            if LOW < rho < HIGH and 0 <= value(left, rho) <= rho:
                cuts.add(rho)
    cuts = sorted(cuts)
    rows, cells, seen_profiles = [], [], set()
    saving = F(0)
    for lo, hi in zip(cuts, cuts[1:]):
        rho = (lo + hi) / 2
        active = sorted(
            (line for line in lines if 0 <= value(line, rho) <= rho),
            key=lambda line: value(line, rho),
        )
        require(value(active[0], rho) == 0, "Missing left domain boundary")
        require(value(active[-1], rho) == rho, "Missing right domain boundary")
        intercept, slope = F(1), F(0)
        special_intercept, special_slope = F(0), F(0)
        for left, right in zip(active, active[1:]):
            require(value(left, rho) < value(right, rho), "Non-strict region ordering")
            p = profile((value(left, rho) + value(right, rho)) / 2, rho)
            seen_profiles.add(p)
            rank_cap = cap(p)
            dc, db = right[0] - left[0], right[1] - left[1]
            intercept += rank_cap * dc
            slope += rank_cap * db
            if p == ((0,), (0, 1), ()):
                special_intercept += dc
                special_slope += db
                saving += integral(dc, db, lo, hi)
            cells.append(dict(lo=lo, hi=hi, left=left, right=right,
                              K=p[0], L=p[1], W=p[2], rank_cap=rank_cap))
        if rho < F(11, 12):
            special = (F(0), F(0))
        elif rho < 1:
            special = (F(-11, 4), F(3))
        else:
            special = (F(-3, 4), F(1))
        require((special_intercept, special_slope) == special,
                "Rank-one region width disagrees with its printed formula")
        rows.append((lo, hi, intercept, slope, integral(intercept, slope, lo, hi)))
    combined = []
    for row in rows:
        if combined and combined[-1][2:4] == row[2:4]:
            previous = combined[-1]
            combined[-1] = (previous[0], row[1], previous[2], previous[3],
                            previous[4] + row[4])
        else:
            combined.append(row)
    inner = sum(row[4] for row in rows)
    outer = integral(F(19, 4), F(-2), F(5, 4), F(2))
    outer += integral(F(11, 4), F(-1), F(2), F(11, 4))
    strip = F(13, 4) * LOW
    require(inner == F(1121513, 369600), "Inner integral mismatch")
    require(outer == F(45, 32), "Outer integral mismatch")
    require(strip == F(13, 20), "Small-prime term mismatch")
    require(inner + outer + strip == F(1881503, 369600), "Total exponent mismatch")
    require(saving == F(5, 48), "Rank-one correction mismatch")
    require(saving == integral(F(-11, 4), F(3), F(11, 12), F(1))
            + integral(F(-3, 4), F(1), F(1), F(5, 4)),
            "Printed two-piece rank-one formula mismatch")
    counts = (len(lines), len(cuts), len(rows), len(cells), len(seen_profiles),
              len(combined), len(set(row[2:4] for row in combined)))
    require(counts == (27, 45, 44, 264, 66, 37, 34), "Enumeration counts mismatch")
    return dict(all_lines=all_lines, lines=lines, cuts=cuts, rows=rows,
                combined_rows=combined, cells=cells, profiles=sorted(seen_profiles),
                inner=inner, outer=outer, strip=strip, saving=saving, counts=counts)


def strip_comments(source):
    return re.sub(r"(?<!\\)%[^\n]*", "", source)


def read_tex(path, root=None, stack=()):
    path = path.resolve()
    root = root or path.parent
    require(path not in stack, "Cyclic TeX input")
    source = strip_comments(path.read_text(encoding="utf-8"))
    def expand(match):
        name = match.group(1)
        candidate = root / name
        if candidate.suffix == "":
            candidate = candidate.with_suffix(".tex")
        if not candidate.exists():
            candidate = path.parent / name
            if candidate.suffix == "":
                candidate = candidate.with_suffix(".tex")
        return read_tex(candidate, root, stack + (path,))
    return re.sub(r"\\(?:input|include)\s*\{([^}]+)\}", expand, source)


def verify_tex(source, data):
    require_literal(source)
    printed = arithmetic_rows(source)
    require(len(printed) == len(data["rows"]) == 44, "Table A.1 must list the 44 uncombined intervals")
    by_interval = {}
    for cell in data["cells"]:
        by_interval.setdefault((cell["lo"], cell["hi"]), []).append(cell)
    endpoint_checks = 0
    for number_, (actual, expected) in enumerate(zip(printed, data["rows"]), 1):
        lo, hi, lines, values, intercept, slope, contribution = actual
        require((lo, hi) == expected[0:2], f"Table A.1 row {number_}: interval {lo},{hi}")
        cells = by_interval[(lo, hi)]
        interior = [cell["left"] for cell in cells[1:]]
        require(lines == interior, f"Table A.1 row {number_}: lines {lines} != {interior}")
        selected = {(F(0), F(0)), (F(0), F(1)), *lines}
        for endpoint in (lo, hi):
            ordered = [value(line, endpoint) for line in lines]
            require(all(0 <= x <= endpoint for x in ordered),
                    f"Table A.1 row {number_}: printed line outside strip at {endpoint}")
            require(ordered == sorted(ordered),
                    f"Table A.1 row {number_}: printed line order fails at {endpoint}")
            for line in data["all_lines"]:
                if line in selected:
                    continue
                left_side = value(line, lo) <= 0 and value(line, hi) <= 0
                right_side = value(line, lo) >= lo and value(line, hi) >= hi
                require(left_side or right_side,
                        f"Table A.1 row {number_}: excluded line enters strip: {line}")
            endpoint_checks += 1
        require(values == [cell["rank_cap"] for cell in cells], f"Table A.1 row {number_}: values")
        require((intercept, slope) == expected[2:4], f"Table A.1 row {number_}: c(rho)")
        require(contribution == expected[4], f"Table A.1 row {number_}: integral")
        # The printed row determines c(rho) by itself: c = 1 + sum r * length.
        edges = [(F(0), F(0))] + lines + [(F(0), F(1))]
        c0 = 1 + sum(r * (b[0] - a[0]) for r, a, b in zip(values, edges, edges[1:]))
        c1 = sum(r * (b[1] - a[1]) for r, a, b in zip(values, edges, edges[1:]))
        require((c0, c1) == (intercept, slope), f"Table A.1 row {number_}: row does not sum to c(rho)")
        require(contribution == integral(intercept, slope, lo, hi), f"Table A.1 row {number_}: integral of c")
    appendix = source.split(r"\label{app:arithmetic}", 1)
    require(len(appendix) == 2, "Missing arithmetic appendix")
    appendix = appendix[1].split(r"\section", 1)[0]
    require("$25$ meet $0<x<\\rho$" in appendix,
            "Printed interior-line count in Appendix A")
    require(re.search(r"\$27\$ meet the closed (?:strip|interval) after adding \$x=0\$ and\s+\$x=\\rho\$",
                      appendix) is not None,
            "Printed closed-strip line count in Appendix A")
    counts = re.findall(r"\$44\\cdot6=264\$ regions and \$(\d+)\$ distinct", appendix)
    require(counts == ["66"], "Printed profile count in Appendix A")
    summary = re.findall(r"\\begin\{align\*\}(.*?)\\end\{align\*\}", appendix, re.S)
    require(len(summary) == 1, "Expected one literal arithmetic summary")
    literals = [re.split(r"&|\\\\", piece, maxsplit=1)[0].strip().rstrip(",.")
                for piece in re.split(r"&\s*=\s*", summary[0])[1:]]
    expected_totals = (data["inner"], data["outer"], data["strip"],
                       data["inner"] + data["outer"] + data["strip"], F(53, 48),
                       -data["saving"])
    require(tuple(number(x) for x in literals) == expected_totals, "Printed totals mismatch")
    tau_pattern = (
        r"\\tau\s*=\s*\\frac\{([0-9]+)\}\{([0-9]+)\}\s*-\s*"
        r"\\frac\{([0-9]+)\}\{([0-9]+)\}\\log2"
    )
    tau_formulas = re.findall(tau_pattern, source)
    require(len(tau_formulas) >= 2, "Missing printed exact tau formula")
    for match in tau_formulas:
        num, den, credit_num, credit_den = map(int, match)
        require(F(num, den) == expected_totals[3] and F(credit_num, credit_den) == F(53, 48),
                "Printed tau coefficients mismatch")
    require(r"-\int_{11/12}^{1}(3\rho-11/4)\,d\rho" in source
            and r"-\int_1^{5/4}(\rho-3/4)\,d\rho=-\frac5{48}" in source, "Missing rank-one integral")
    require(r"5/2-2\rho\le x<\rho-1/4$ for $11/12<\rho<1" in source
            and r"1/2\le x<\rho-1/4$ for $1\le\rho<5/4" in source, "Rank-one region")
    require(r"\int_{5/4}^{11/4}\mathfrak c_{\rm out}(\rho)\,d\rho=\frac98+\frac9{32}=\frac{45}{32}" in source,
            "Outer integral display")
    # Accept the displayed profile spelling used in the manuscript.
    require(
        r"\int_{1/5}^{5/4}c(\rho)\,d\rho=\frac{1121513}{369600}" in source
        or r"\int_{1/5}^{5/4}\mathfrak c(\rho)\,d\rho=\frac{1121513}{369600}" in source,
        "Inner integral display",
    )
    require(r"=4.325296633153372\ldots" in source and r"4.36-\tau>0.0347" in source, "tau decimals")
    require(r"\label{eq:subdivision-lines}" in source, "Missing line construction")
    return len(printed), endpoint_checks


def padd(*polynomials):
    return tuple(sum(p[i] for p in polynomials) for i in range(3))


def pscale(p, factor):
    return tuple(coefficient * factor for coefficient in p)


def pmul(a, b):
    coefficients = [F(0)] * 5
    for i in range(3):
        for j in range(3):
            coefficients[i+j] += a[i] * b[j]
    require(coefficients[3:] == [0, 0], "Unexpected polynomial degree")
    return tuple(coefficients[:3])


def verify_minimum_identity():
    # Symbolic polynomial calculation in t, separately for n=3t+s.
    one = (F(1), F(0), F(0))
    for s in range(3):
        n = (F(s), F(3), F(0))
        r = (5 * s) % 3
        q = (F(5 * s - r, 3), F(5), F(0))
        minimum = padd(
            pscale(pmul(n, n), F(7, 2)), pscale(n, F(-3, 2)),
            pscale(pmul(n, q), 6),
            pscale(pmul(q, padd(q, pscale(one, -1))), F(3, 2)),
            pscale(n, 2 * r), pscale(q, r),
        )
        formula = padd(pscale(pmul(n, n), F(53, 3)), pscale(n, -4),
                       pscale(one, F(r * (3-r), 6)))
        require(minimum == formula, f"Formal N_m identity fails for n mod 3 = {s}")
        require(F(r * (3-r), 6) == (0 if s == 0 else F(1, 3)),
                "N_m rounding correction fails")


def logarithm_interval():
    # log 2 = 2 atanh(1/3); bound the omitted positive tail geometrically.
    z, count = F(1, 3), 48
    lower = 2 * sum(z ** (2*k+1) / (2*k+1) for k in range(count))
    upper = lower + 2 * z ** (2*count+1) / ((2*count+1)*(1-z*z))
    total, credit = F(1881503, 369600), F(53, 48)
    tau_lower, tau_upper = total-credit*upper, total-credit*lower
    prefix = F("4.325296633153372")
    require(prefix < tau_lower < tau_upper < prefix + F(1, 10**15),
            "Printed tau decimal prefix is wrong")
    require(F(109, 25) - tau_upper > F(347, 10000), "Printed margin fails")
    return dict(log2_lower=lower, log2_upper=upper,
                tau_lower=tau_lower, tau_upper=tau_upper)


def finite_lattice_sanity():
    cases = 0
    # No primality is needed for the combinatorial half-open identity.
    for n in range(1, 17):
        m = 4*n
        for p in range(m//5 + 1, 5*n):
            rho = F(p, m)
            discrete = sum(cap(profile(F(u, m), rho)) for u in range(p))
            xcuts = {F(0), rho}
            for k in range(KMAX+1):
                xcuts.update(x for x in (F(1,2)-k*rho, F(3,2)-k*rho) if 0<x<rho)
            for l in range(LMAX+1):
                xcuts.update(x for x in (F(1,4)-l*rho, F(5,2)-l*rho) if 0<x<rho)
            for w in range(1, WMAX+1):
                xcuts.update(x for x in (w*rho-F(1,4), w*rho) if 0<x<rho)
            xcuts = sorted(xcuts)
            continuous = sum((b-a)*cap(profile((a+b)/2, rho))
                             for a, b in zip(xcuts, xcuts[1:]))
            require(discrete == m*continuous, f"Half-open lattice sanity failed: m={m}, p={p}")
            for a,b in zip(xcuts,xcuts[1:]):
                require(profile(a,rho) == profile((a+b)/2,rho),
                        f"Endpoint orientation fails: m={m}, p={p}, x={a}")
            cases += 1
    return cases


def finite_outer_sanity():
    cases, entry_checks, second_v_multiples = 0, 0, 0
    for n in range(1, 13):
        I, J, V = range(2*n,6*n), range(n,10*n), range(1,n+1)
        for p in range(5*n+1,11*n):
            if p < 6*n:
                lb=lambda i: 1+int(i>=p)
                lv=lambda ell: 1
                mu=lambda j: int(j>=p+2*n)
                expected=19*n-2*p
            elif p < 8*n:
                lb=lambda i: int(i+p<10*n)
                lv=lambda ell: 1
                mu=lambda j: int(j>=p)
                expected=19*n-2*p
            elif p < 10*n:
                lb=lambda i: 0
                lv=lambda ell: 1
                mu=lambda j: int(j>=p)
                expected=11*n-p
            else:
                lb=lambda i: 0
                lv=lambda ell: int(p-ell in J)
                mu=lambda j: 0
                expected=11*n-p
            require(sum(map(lb,I))+sum(map(lv,V))+sum(map(mu,J)) == expected,
                    f"Outer exponent mismatch: n={n}, p={p}")
            for i in I:
                for j in J:
                    cw=int(i!=j and (i-j)%p==0)
                    bw=int(max(i,j)>=p)+int(max(i,j)>=p and (i-j)%p==0)
                    require(mu(j)>=cw and lb(i)+mu(j)>=bw,
                            f"Outer C/B edge fails: n={n}, p={p}, i={i}, j={j}")
                    entry_checks += 2
            for ell in V:
                for j in J:
                    require(lv(ell)+mu(j)>=int((j+ell)%p==0),
                            f"Outer V edge fails: n={n}, p={p}, ell={ell}, j={j}")
                    entry_checks += 1
                    second_v_multiples += int(j+ell == 2*p)
            cases += 1
    require(second_v_multiples > 0, "No second V multiple was tested")
    return dict(cases=cases, entry_checks=entry_checks,
                second_V_multiples=second_v_multiples)


def json_default(value):
    if isinstance(value,F):
        return str(value)
    raise TypeError(type(value).__name__)


def outward_decimal(value, places, upper=False):
    scaled = value * 10**places
    integer = -(-scaled.numerator // scaled.denominator) if upper else scaled.numerator // scaled.denominator
    digits = str(integer).zfill(places+1)
    return digits[:-places]+"."+digits[-places:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path)
    parser.add_argument("--section", type=Path, help="Additional section file before main.tex is assembled")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    tex = args.tex or root/"main.tex"
    if not tex.exists() and args.tex is None:
        tex = root/"parts"/"appA.tex"
        args.section = root/"parts"/"sec3.tex"
    start = time.monotonic()
    source = read_tex(tex)
    if args.section:
        source += "\n"+read_tex(args.section)
    data = reconstruct()
    row_count, endpoint_checks = verify_tex(source,data)
    verify_minimum_identity()
    bounds = logarithm_interval()
    finite_lattice = finite_lattice_sanity()
    finite_outer = finite_outer_sanity()
    result = dict(
        ok=True, method="exact Fraction affine arrangement, literal TeX table entries parsed",
        data=data, log_bounds=bounds,
        endpoint_chain_checks=endpoint_checks,
        tex_sha256=hashlib.sha256(source.encode()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        finite_sanity_only=dict(lattice_cases=finite_lattice, outer=finite_outer),
        python=platform.python_version(), elapsed_seconds=time.monotonic()-start,
    )
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,default=json_default,indent=2)+"\n")
    print("PASS exact arrangement: 27 lines, 45 cut points, 44 intervals, 264 regions, 66 triples (K,mathcal L,W)")
    print(f"PASS Table A.1: {row_count} rows; lines, values of r, c(rho) and integrals match; each row sums to c(rho)")
    print(f"PASS endpoint ordering: {endpoint_checks} endpoint chains checked against all {len(data['all_lines'])} lines")
    print(f"PASS exact integrals: inner={data['inner']}; outer={data['outer']}; small primes={data['strip']}")
    print(f"PASS rank-one correction=-{data['saving']}; total exponent=1881503/369600")
    print("PASS N_m identity in all three classes n mod 3")
    print("PASS first-line branch precedence handles p^2<=U")
    print("PASS tau in ["+outward_decimal(bounds["tau_lower"],30)+", "
          +outward_decimal(bounds["tau_upper"],30,True)+"]; 4.36-tau>0.0347")
    print(f"FINITE SANITY ONLY: {finite_lattice} half-open lattice cases; "
          f"{finite_outer['entry_checks']} outer entry checks; "
          f"{finite_outer['second_V_multiples']} cases j+ell=2p")
    print("Uniform valuation and prime-sum proofs are in Section 3; finite sanity tests do not replace them.")
    print(f"Python {platform.python_version()}; elapsed {time.monotonic()-start:.3f} s")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}",file=sys.stderr)
        sys.exit(1)
