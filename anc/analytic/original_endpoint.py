#!/usr/bin/env python3
"""Validated full-disk endpoint bounds for the merged q=2 C_joint certificate.

The script is intentionally independent of the source cap verifier.  It only
uses the frozen rational reference and 100-bit Arb interval arithmetic.  The
complex-disk step is analytic: triangle inequalities reduce every supremum or
infimum over |q|<=delta to a real endpoint or to the Poisson radius bound.
"""
from __future__ import annotations
import hashlib, json, math, platform, sys
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, ctx

ctx.prec = 100
ROOT = Path(__file__).resolve().parent
# The frozen bundle stores the reference JSON files beside this verifier.
# Earlier drafts used an ``input/`` subdirectory, which made a clean replay
# fail even though the certificate data were present.
IN = ROOT
DELTA = arb(1) / arb(250)  # 0.004 exactly
RADIUS = arb(3) / arb(5)
M0 = 64
PI_LOWER = Fraction(333, 106)

# A short exact certificate for the only pi bound used below.  Machin's
# identity and the alternating arctan series give
# pi > 16*(1/5 - 1/(3*5^3) + 1/(5*5^5) - 1/(7*5^7)) - 4/239
#     = 333/106 + 3418213/41563593750.
_machin_lower = (16 * (Fraction(1, 5) - Fraction(1, 3 * 5**3)
                       + Fraction(1, 5 * 5**5) - Fraction(1, 7 * 5**7))
                 - Fraction(4, 239))
assert _machin_lower > PI_LOWER

def rat(s: str) -> arb:
    f = Fraction(s)
    return arb(f.numerator) / arb(f.denominator)

def parse_interval(row):
    return rat(row["a"]), rat(row["b"]), rat(row["w"])

def A(u: acb) -> acb:
    # Real part is the primitive of log|q-t|; branch choices do not affect it.
    if abs(u) < arb("1e-80"):
        return acb(0)
    return u * u.log() - u

def interval_potential(q, arr):
    q = acb(q)
    out = arb(0)
    for a, b, w in arr:
        out += w / (b - a) * (A(q-a) - A(q-b)).real
    return out

def exact_mass(rows):
    return sum((Fraction(x["w"]) for x in rows), Fraction(0))

def ceil_fraction(q: Fraction) -> int:
    """Exact ceiling for a nonnegative rational."""
    if q < 0:
        return -((-q.numerator) // q.denominator)
    return (q.numerator + q.denominator - 1) // q.denominator

def ceil_sqrt_integer(n: int) -> int:
    """Smallest integer s with s^2 >= n."""
    r = math.isqrt(n)
    return r if r * r == n else r + 1

def sqrt_upper(q: Fraction) -> Fraction:
    """Rational upper bound for sqrt(q), exact and strict when needed."""
    if q < 0:
        raise ValueError("sqrt_upper expects a nonnegative rational")
    if q == 0:
        return Fraction(0)
    # sqrt(num/den) = sqrt(num*den)/den.
    return Fraction(ceil_sqrt_integer(q.numerator * q.denominator), q.denominator)

def density_budget(data):
    def dmax(rows):
        return max(Fraction(x["w"]) / (Fraction(x["b"]) - Fraction(x["a"])) for x in rows)
    mx_exact = dmax(data["x"])
    mv_exact = dmax(data["y"])
    mzr_exact = dmax(data["z_real"])
    mx = ceil_fraction(mx_exact)
    mv = ceil_fraction(mv_exact)
    mzr = ceil_fraction(mzr_exact)
    r = Fraction(data["radius"])
    mzp_upper = Fraction(0)
    for x in data["z_poisson"]:
        re = Fraction(x["re"])
        im = Fraction(x["im"])
        w = Fraction(x["w"])
        c_hi = sqrt_upper(re * re + im * im)
        if not c_hi < r:
            raise AssertionError(("Poisson center bound", c_hi, r))
        # The exact rational bound 333/106 < pi gives
        # 2*pi > 333/53.  With c<=c_hi, the displayed rational expression is
        # an outward upper bound for the arclength density.  All operations
        # below are exact Fraction arithmetic.
        mzp_upper += w * (r + c_hi) / (2 * PI_LOWER * r * (r - c_hi))
    mzp_ceiling = ceil_fraction(mzp_upper)
    mz = mzr + mzp_ceiling
    Kmod = 8
    # The three fields carry different source masses.  The x field sees only
    # its self interaction; v sees self plus v-z; z sees self plus z-v.
    # Keep these as separate integer coefficients before applying masses.
    cx = 2 * Kmod * mx
    cv = 2 * Kmod * mv + Kmod * mz
    cz = 2 * Kmod * mz + Kmod * mv
    weighted = Fraction(cx) + Fraction(5, 4) * cv + Fraction(cz)
    if weighted.denominator != 1:
        raise AssertionError(("weighted density coefficient is not integral", weighted))
    # Retain the old maximum only as a diagnostic for comparison with older
    # reports; the proof uses the mass-weighted coefficient below.
    ccap = max(cx + Kmod * mz, cv, cz)
    delta_coeff = 3 * weighted + 200
    if delta_coeff.denominator != 1:
        raise AssertionError(("delta coefficient is not integral", delta_coeff))
    return {
        "Mx_exact": str(mx_exact), "Mv_exact": str(mv_exact), "Mzreal_exact": str(mzr_exact),
        "Mx": mx, "Mv": mv, "Mzreal": mzr,
        "Mz_circle_upper": str(mzp_upper), "Mz_circle_ceiling": mzp_ceiling, "Mz": mz,
        "Kmod": Kmod, "Ccap_legacy_max": ccap, "delta_coeff": delta_coeff.numerator,
        "Cx": cx, "Cv": cv, "Cz": cz,
        "weighted_density_coeff": weighted.numerator,
    }

def cstar(B):
    T = 3.7302903070270859171495181793204702746102129224222
    m = 64.0
    LD = 1.01; LZ = 1.67
    return (B*(1+math.log(m*m/3))/(m*math.log(m))
            + 6.5 + math.log(T)/math.log(m)
            + (6*1.25*LD + 3*1.25*LZ)/(m*math.log(m)))

def run(tag: str, filename: str, caps):
    data = json.loads((IN / filename).read_text())
    X = [parse_interval(x) for x in data["x"]]
    V = [parse_interval(x) for x in data["y"]]
    ZR = [parse_interval(x) for x in data["z_real"]]
    ZP = [(rat(x["re"]), rat(x["im"]), rat(x["w"])) for x in data["z_poisson"]]
    radius = rat(data["radius"])
    radius_frac = Fraction(data["radius"])

    # Shape/mass checks are exact rational assertions, not sampled checks.
    assert exact_mass(data["x"]) == Fraction(1)
    assert exact_mass(data["y"]) == Fraction(5, 4)
    assert all(Fraction(x["a"]) >= 0 and Fraction(x["b"]) <= 1 for x in data["x"])
    assert all(Fraction(x["a"]) >= 0 and Fraction(x["b"]) <= 1 for x in data["y"])
    # The frozen z_real reference is the negative real path only; the positive
    # real chain is an evaluation support, not a reference-density component.
    assert all(Fraction(x["a"]) >= -2 and Fraction(x["b"]) <= -Fraction(3,5) for x in data["z_real"])
    assert all(Fraction(x["a"]) < -Fraction(3,5) and Fraction(x["b"]) <= -Fraction(3,5) for x in data["z_real"])
    for re, im, _w in ZP:
        assert re*re + im*im < radius*radius
    assert radius_frac == Fraction(3, 5)

    mx = arb(exact_mass(data["x"]).numerator) / arb(exact_mass(data["x"]).denominator)
    mv = arb(exact_mass(data["y"]).numerator) / arb(exact_mass(data["y"]).denominator)
    mzp_frac = sum((Fraction(x["w"]) for x in data["z_poisson"]), Fraction(0))
    mzp = arb(mzp_frac.numerator) / arb(mzp_frac.denominator)

    # Endpoint primitives.  The signs match the disk extremizers:
    # x,v supports [0,1] -> q=-delta; negative z path -> q=+delta.
    ux_star = interval_potential(-DELTA, X)
    uv_star = interval_potential(-DELTA, V)
    uz_real_star = interval_potential(DELTA, ZR)
    # For |q|<RADIUS, each Poisson component is
    # log|R^2-q*conj(c)| - log R.  Since |c|<R and |q|<=delta,
    # this is at most log(R+delta), by the triangle inequality.
    uz_star = uz_real_star + mzp * (RADIUS + DELTA).log()

    # Denominator lower bounds use |2+qt| >= 2-|q|t >= 2-delta.
    # The reference masses are exact, so these are interval-certified.
    d_v_lower = mv * (arb(2) - DELTA).log()
    d_x_lower = mx * (arb(2) - DELTA).log()
    log_delta = DELTA.log()
    fx_upper = 2*ux_star - d_v_lower + arb(3)/4 * log_delta
    fv_upper = 2*uv_star + uz_star - d_x_lower + arb(1)/4 * log_delta

    kx, kv, kz = map(arb, caps)
    slack_x = kx - fx_upper
    slack_v = kv - fv_upper
    assert fx_upper < kx, (tag, "Fx cap failure", fx_upper, kx)
    assert fv_upper < kv, (tag, "Fv cap failure", fv_upper, kv)

    # Smoothing geometry for epsilon=m^-2.  A translate by 2i epsilon and a
    # radius-epsilon circle move every center by at most 3 epsilon.
    eps64 = arb(1) / arb(M0*M0)
    three_eps = 3*eps64
    six_eps = 6*eps64
    assert three_eps < DELTA/2
    assert six_eps < DELTA/2

    db = density_budget(data)
    return {
        "tag": tag,
        "reference": filename,
        "reference_sha256": hashlib.sha256((IN/filename).read_bytes()).hexdigest(),
        "delta": str(DELTA),
        "radius": str(RADIUS),
        "m0": M0,
        "masses": {"x": str(mx), "v": str(mv), "z_poisson": str(mzp)},
        "support_assertions": {
            "x_v_in_0_1": True,
            "z_real_in_neg_path": True,
            "poisson_centers_strictly_inside": True,
        },
        "endpoint_values": {
            "Ux_minus_delta": str(ux_star),
            "Uv_minus_delta": str(uv_star),
            "Uzreal_plus_delta": str(uz_real_star),
            "Uz_disk_upper": str(uz_star),
            "Dv_disk_lower": str(d_v_lower),
            "Dx_disk_lower": str(d_x_lower),
            "Fx_disk_upper": str(fx_upper),
            "Fv_disk_upper": str(fv_upper),
            "Fx_slack": str(slack_x),
            "Fv_slack": str(slack_v),
        },
        "smearing_geometry_m0": {
            "epsilon_64": str(eps64),
            "3epsilon": str(three_eps),
            "6epsilon": str(six_eps),
            "delta_over_2": str(DELTA/2),
            "endpoint_outside_log_loss_bound": "6*epsilon/delta",
        },
        "caps": {"Kx": caps[0], "Kv": caps[1], "Kz": caps[2]},
        "density_budget": db,
        "C_star_check": {"B_frozen_coarse": 26000 if tag == "cap-tightened" else 27800,
                          "C_star_at_B": cstar(26000 if tag == "cap-tightened" else 27800),
                          "C_star_at_direct_delta_coeff": cstar(db["delta_coeff"])},
    }

def main():
    out = {
        "script": "verify_endpoint_disk.py",
        "date": "2026-10-05",
        "precision_bits": ctx.prec,
        "python": sys.version,
        "flint_version": __import__("flint").__version__,
        "platform": platform.platform(),
        "analytic_claim": {
            "disk": "|q| <= delta",
            "delta": "1/250",
            "Ux_Uv": "|q-t| <= t+delta = |-delta-t| for t in [0,1]",
            "denominator": "|2+q*t| >= 2-delta*t >= 2-delta",
            "z_real": "z_real subset [-2,-3/5], |q-t| <= |t|+delta = |delta-t|",
            "z_poisson": "log|R^2-q*conj(c)|-log R <= log(R+delta), |c|<R",
            "log_field": "log|q| <= log(delta) for q != 0; q=0 gives -infinity",
        },
        "runs": [
            run("cap-tightened", "endpoint-tight_reference.json", ("-4.6277", "-4.4816", "-0.092")),
            run("optimized", "endpoint-optimized_reference.json", ("-4.629", "-4.485", "-0.092")),
        ],
    }
    (ROOT / "endpoint-RESULTS.json").write_text(json.dumps(out, indent=2) + "\n")
    for r in out["runs"]:
        e = r["endpoint_values"]
        print(r["tag"])
        for k in ("Ux_minus_delta", "Uv_minus_delta", "Uzreal_plus_delta", "Uz_disk_upper", "Dv_disk_lower", "Dx_disk_lower", "Fx_disk_upper", "Fv_disk_upper", "Fx_slack", "Fv_slack"):
            print(f"  {k} = {e[k]}")
        print("  m0 geometry:", r["smearing_geometry_m0"])
        print("  density budget:", r["density_budget"])
        print("  C* check:", r["C_star_check"])

if __name__ == "__main__":
    main()
