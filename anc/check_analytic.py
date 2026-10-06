#!/usr/bin/env python3
"""Interval-arithmetic proof of Lemma 4.2 and check of the printed constants.

The subdivision and endpoint routines (analytic/original_*.py) are unchanged
copies of the original implementation.  This wrapper also checks the exact
geometry and masses of the measure, the literal entries of Tables B.1 and B.2,
the enclosures (B.1)-(B.2), and the rational constants of Lemma 4.4.
The passage to all matrix sizes is the written proof of Lemma 4.4.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import sys
import time

from tex_tables import analytic_data, compact, depth_sequence_runs, number, read_tex, table

sys.dont_write_bytecode = True

import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
DATA = HERE / "analytic"
EXPECTED_SHA = {
    "reference.json": "4e72d95479eb0a9632852aabe644fb17741ecc56ee9f960c8edc2fe193b7201b",
    "original_p3.py": "a557130cf07c85ec7df8e158e5537a3d147193a9b741ccaf45d4a8b495e07303",
    "original_endpoint.py": "13d7900affa247a432d48e921137279ace1d99b146f95fb9b34c22d76dfd2973",
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def encode_depth_runs(depths):
    runs = []
    for depth in depths:
        if runs and runs[-1][0] == depth:
            runs[-1] = (depth, runs[-1][1] + 1)
        else:
            runs.append((depth, 1))
    return runs


def depth_table_tex(leaves):
    labels = {"x": r"$x$", "v": r"$v$", "zn": r"$z_-$",
              "zp": r"$z_+$", "zc": r"$z_\circ$"}
    out = [r"\begin{longtable}{@{}l r p{0.79\textwidth}@{}}",
           r"\caption{Run-length encoded dyadic depths of the certificate partitions.}\label{tab:depth-sequences}\\",
           r"\hline support & block & depth:repeat \\", r"\hline", r"\endfirsthead",
           r"\hline support & block & depth:repeat \\", r"\hline", r"\endhead"]
    for key in ("x", "v", "zn", "zp", "zc"):
        runs = encode_depth_runs([item["depth"] for item in leaves[key]])
        for block, start in enumerate(range(0, len(runs), 7), 1):
            tokens = ",".join(f"{d}:{count}" for d, count in runs[start:start+7])
            out.append(f"{labels[key]} & {block} & \\texttt{{{tokens}}} \\\\")
    out.extend([r"\hline", r"\end{longtable}"])
    return "\n".join(out) + "\n"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def A(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def tex_data(root, explicit=None):
    if explicit is not None:
        files = [explicit]
    elif (root / "main.tex").is_file():
        files = [root / "main.tex"]
    else:
        files = [root / "parts" / "sec4.tex", root / "parts" / "appB.tex"]
    text = "\n".join(read_tex(p) for p in files)
    values, rows = analytic_data(text)
    return files, values, rows


def check_reference(data):
    masses = {}
    for key, mass, bounds in [("x", F(1), (0, 1)), ("y", F(5, 4), (0, 1)),
                               ("z_real", None, (-2, F(-3, 5)))]:
        rows = sorted((F(t["a"]), F(t["b"]), F(t["w"])) for t in data[key])
        require(all(bounds[0] <= a < b <= bounds[1] and w > 0 for a, b, w in rows), (key, "support or positivity"))
        require(all(b <= c for (_, b, _), (c, _, _) in zip(rows, rows[1:])), (key, "overlapping density intervals"))
        masses[key] = sum((w for a, b, w in rows), F())
        if mass is not None:
            require(masses[key] == mass, (key, "mass"))
    r = F(data["radius"])
    require(r == F(3, 5), "reference radius")
    poisson = Counter((F(t["re"]), F(t["im"]), F(t["w"])) for t in data["z_poisson"])
    require(all(w > 0 and a*a + b*b < r*r for a, b, w in poisson), "Poisson positivity or pole location")
    require(poisson == Counter((a, -b, w) for a, b, w in poisson.elements()), "exact conjugation symmetry")
    masses["z_poisson"] = sum((w for a, b, w in poisson.elements()), F())
    require(masses["z_real"] + masses["z_poisson"] == 1, "z mass")
    require([len(data[k]) for k in ("x", "y", "z_real", "z_poisson")] == [50, 52, 62, 157], "original reference component counts")
    return {k: str(v) for k, v in masses.items()}


def literal_fraction(cell):
    """Parse a rational cell from the printed reference tables."""
    text = compact(cell)
    match = re.fullmatch(r"-\\frac\{(\d+)\}\{(\d+)\}", text)
    if match:
        return -F(int(match[1]), int(match[2]))
    return number(cell)


def check_reference_tables(source, data):
    """Require the printed Appendix D rows to equal reference.json exactly."""
    interval_rows = {"x": [], "v": [], "z": []}
    for row in table(source, "tab:reference-intervals"):
        if len(row) != 5 or not compact(row[0]).isdigit():
            continue
        support = compact(row[1])
        require(support in interval_rows, ("printed interval support", support))
        interval_rows[support].append(tuple(literal_fraction(x) for x in row[2:]))
    for support, key in (("x", "x"), ("v", "y"), ("z", "z_real")):
        expected = [(F(x["a"]), F(x["b"]), F(x["w"])) for x in data[key]]
        require(interval_rows[support] == expected,
                ("printed interval table differs from reference.json", support))

    poisson_rows = []
    for row in table(source, "tab:reference-poisson"):
        if len(row) != 5 or not compact(row[0]).isdigit():
            continue
        poisson_rows.append(tuple(literal_fraction(x) for x in row[2:]))
    expected_poisson = [(F(x["re"]), F(x["im"]), F(x["w"]))
                       for x in data["z_poisson"]]
    require(poisson_rows == expected_poisson,
            "printed Poisson table differs from reference.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, help="Read the delivered standalone TeX file.")
    parser.add_argument("--tex-root", type=Path, default=HERE.parent)
    parser.add_argument("--out", type=Path, default=Path("/tmp/li2-analytic-replay"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    files, printed, rows = tex_data(args.tex_root, args.tex)
    for name, sha in EXPECTED_SHA.items():
        require(hashlib.sha256((DATA/name).read_bytes()).hexdigest() == sha, (name, "SHA-256"))
    data = json.loads((DATA / "reference.json").read_text())
    source = "\n".join(read_tex(path) for path in files)
    check_reference_tables(source, data)
    print("Printed reference tables == reference.json: PASS", flush=True)
    masses = check_reference(data)
    os.environ["P3_REFERENCE"] = str(DATA / "reference.json")
    p3 = load_module("analytic_original_p3", DATA / "original_p3.py")
    endpoint = load_module("analytic_original_endpoint", DATA / "original_endpoint.py")
    ctx.dps = 90
    require(ctx.dps == int(printed["precision"]), "printed precision")
    print("Original reference geometry and exact masses: PASS", flush=True)
    q, parts = p3.compute_Q()
    qlo, qhi = F(printed["qlo"]), F(printed["qhi"])
    require(A(qlo) < q.lower() and q.upper() < A(qhi), "printed rational Q enclosure")
    print("Q enclosure:", q, flush=True)
    caps = {"x": F(printed["Kx"]), "v": F(printed["Kv"]), "z": F(printed["Kz"])}
    require((caps["x"], caps["v"], caps["z"]) == (p3.CAPX, p3.CAPV, p3.CAPZ), "original caps")
    leaves = {}
    reports = {}
    ranges = {"x": (F(0), F(1)), "v": (F(0), F(1)), "zn": (F(-2), F(-3,5)),
              "zp": (F(3,5), F(1)), "zc": (F(0), F(1))}
    labels = {"x": r"$x\in[0,1]$", "v": r"$v\in[0,1]$",
              "zn": r"$z\in[-2,-3/5]$", "zp": r"$z\in[3/5,1]$",
              "zc": r"$|z|=3/5,\ \operatorname{Im}z\ge0$"}
    for kind in ("x", "v", "zn", "zp", "zc"):
        target = caps[kind] if kind in ("x", "v") else caps["z"]
        accepted, maxdepth, rejected = p3.bisect_circle(target) if kind == "zc" else p3.bisect_real(kind, target)
        accepted.sort(key=lambda t: t[0])
        require(accepted[0][0] == ranges[kind][0] and accepted[-1][1] == ranges[kind][1], (kind, "coverage endpoints"))
        require(all(a < b for a, b, d, v in accepted), (kind, "nonempty leaves"))
        require(all(b == c for (a, b, d, v), (c, e, f, w) in zip(accepted, accepted[1:])), (kind, "complete partition"))
        require(all(v.upper() < A(target) for a, b, d, v in accepted), (kind, "strict leaf comparison"))
        worst = max((v for a, b, d, v in accepted), key=lambda v: v.upper())
        row = rows[kind]
        require(re.sub(r"\s+", "", row["label"]) == re.sub(r"\s+", "", labels[kind]), (kind, "printed support"))
        require((len(accepted), maxdepth) == (row["count"], row["depth"]), (kind, "printed leaf counts/depths"))
        require(worst.upper() < A(row["upper"]) < A(row["cap"]), (kind, "printed outward bound"))
        require(row["cap"] == target, (kind, "printed cap"))
        reports[kind] = dict(leaves=len(accepted), depth=maxdepth, upper=str(worst), rejected=rejected)
        leaves[kind] = [dict(a=str(a), b=str(b), depth=d, bound=str(v)) for a, b, d, v in accepted]
        print(kind, len(accepted), "leaves, strict cap", target, "PASS", flush=True)
    printed_runs = depth_sequence_runs(source)
    for kind in ("x", "v", "zn", "zp", "zc"):
        items = leaves[kind]
        expected_runs = encode_depth_runs([item["depth"] for item in items])
        require(printed_runs[kind] == expected_runs, (kind, "printed depth sequence"))
        require(sum(F(count, 2**depth) for depth, count in printed_runs[kind]) == 1,
                (kind, "normalized dyadic lengths sum to one"))
        lo, hi = ranges[kind]
        cursor = lo
        for item in items:
            a, b, depth = F(item["a"]), F(item["b"]), item["depth"]
            require(a == cursor and b-a == (hi-lo)/2**depth,
                    (kind, "depth reconstructs contiguous interval"))
            cursor = b
        require(cursor == hi, (kind, "depth sequence covers full support"))
    (args.out / "depth_sequences.tex").write_text(depth_table_tex(leaves), encoding="utf-8")
    require(sum(r["leaves"] for r in reports.values()) == int(printed["leaves"]) == 3791, "total printed subinterval count")
    require((F(printed["rate"]), int(printed["error_coeff"]), int(printed["B"])) == (F(109, 25), 870, 26000),
            "printed final constants")
    alpha = q - A(F(3,2))*arb(2).log() - A(caps["x"] + F(5,4)*caps["v"] + caps["z"])
    alpha_lo = A(qlo) - A(F(3,2))*arb(2).log() - A(caps["x"] + F(5,4)*caps["v"] + caps["z"])
    require(alpha_lo > A(printed["alpha_lower"]) > A(printed["rate"]), "dual bound uses Q lower endpoint")
    require(alpha < A(printed["alpha_upper"]), "alpha upper endpoint")
    endpoint.IN = DATA
    require(int(printed["endpoint_precision_bits"]) == 100,
            "printed endpoint-disk precision")
    ctx.prec = 100
    ep = endpoint.run("original", "reference.json", tuple(str(caps[k]) for k in ("x", "v", "z")))
    # The original endpoint routine contains an unused diagnostic C_* field
    # for an optimized reference.  It is discarded; all constants used here
    # are recomputed below with directed Arb arithmetic.
    ep.pop("C_star_check")
    require(arb(ep["endpoint_values"]["Fx_disk_upper"]) < A(printed["endpoint_x"]), "printed x endpoint upper bound")
    require(arb(ep["endpoint_values"]["Fv_disk_upper"]) < A(printed["endpoint_v"]), "printed v endpoint upper bound")
    db = ep["density_budget"]
    for key in ("Mx", "Mv", "Mzreal", "Mz_circle_ceiling", "Mz", "Kmod", "Cx", "Cv", "Cz", "weighted_density_coeff", "delta_coeff"):
        require(db[key] == int(printed[key]), (key, "printed density constant"))
    require(db["delta_coeff"] == 3*db["weighted_density_coeff"] + 200 < int(printed["B"]), "B budget")
    ctx.dps = 90
    m0, delta, eps = F(printed["m0"]), F(printed["delta"]), F(1, int(printed["m0"])**2)
    LD, Lz = F(printed["LD"]), F(printed["Lz"])
    require((m0, delta, LD, Lz) == (64, F(1,250), F(101,100), F(167,100)), "source smoothing constants")
    require(6*eps < delta/2 and 3*eps < F(1,16), "endpoint and tubular geometry")
    h = 3*eps
    require((1+h)/(2-2*h-h*h) < LD, "denominator derivative bound")
    require(1/(F(3,5)-h) < Lz, "z logarithm derivative bound")
    endpoint_allowance = F(17,16)*8*eps/(delta*delta)
    smooth_allowance = 6*F(5,4)*LD + 3*F(5,4)*Lz
    require(smooth_allowance == F(printed["A0"]) == F(1107, 80), "printed A_0")
    require(endpoint_allowance + smooth_allowance < int(printed["allowance"]), "full endpoint/smooth allowance")
    require(str(F(13,2)) == printed["self_coeff"], "self-circle coefficient")
    pi = arb.pi()
    T = A(F(9,5)) + A(F(3,5))*(pi*pi + arb(2).log()**2).sqrt()
    require(A(printed["Tlo"]) < T < A(printed["Thi"]), "printed path variation interval")
    B = A(printed["B"])
    cstar = (B*(1+(A(m0*m0)/3).log())/(A(m0)*A(m0).log())
             + A(F(13,2)) + T.log()/A(m0).log()
             + A(smooth_allowance)/(A(m0)*A(m0).log()))
    require(cstar < A(printed["Cstar_upper"]) < A(printed["error_coeff"]), "uniform 870 coefficient")
    # Derivative of (1+log(m^2/3))/(m log m) has the negative
    # numerator -2*(log m)^2 +(log 3-1)*(log m+1).
    ell = A(m0).log()
    require(-2*ell*ell + (arb(3).log()-1)*(ell+1) < 0, "start of decreasing coefficient")
    require(4*ell > arb(3).log()-1, "derivative numerator remains decreasing")
    for key, expected in [("n_x", 50), ("n_v", 52), ("n_zreal", 62), ("n_poisson",157), ("series_N",420)]:
        require(int(printed[key]) == expected, (key, "printed source count"))
    tail = 2*A(F(1,2))**421/421
    require(tail < A(F(printed["series_tail_mantissa"])*F(10)**int(printed["series_tail_power"])), "printed per-pair series remainder")
    tex_hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    depth_source = files[0].parent / "anc" / "analytic" / "depth_sequences.tex"
    require(depth_source.is_file(), "printed depth table source is present")
    tex_hashes[str(depth_source)] = hashlib.sha256(depth_source.read_bytes()).hexdigest()
    report = dict(ok=True, seconds=time.monotonic()-started, python=sys.version,
                  python_flint=flint.__version__, platform=platform.platform(),
                  precision_dps=90, endpoint_precision_bits=100,
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.suffix in (".py", ".json")},
                  tex_sha256=tex_hashes,
                  scope="Finite interval certificate plus exact constants; uniform lemmas are proved in Section 4.",
                  masses=masses, Q=str(q), Q_parts={k:str(v) for k,v in parts.items()},
                  alpha=str(alpha), alpha_from_Q_lower=str(alpha_lo), supports=reports,
                  endpoint=ep, endpoint_allowance=str(endpoint_allowance),
                  smooth_allowance=str(smooth_allowance), path_variation=str(T), C_star=str(cstar))
    (args.out / "result.json").write_text(json.dumps(report, indent=2)+"\n")
    # Keep the complete accepted partition, not just its counts.
    with (args.out / "leaves.json.gz").open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as dst:
            dst.write((json.dumps(leaves, separators=(",", ":"))+"\n").encode())
    print("Endpoint disk and exact density constants: PASS")
    print("alpha =", alpha)
    print("C_* =", cstar)
    print("Printed TeX tables: PASS")
    print("Uniform smoothing formulas: exact constants PASS; their proof is in Section 4.")
    print(f"PASS ({report['seconds']:.3f} seconds)")


if __name__ == "__main__":
    main()
