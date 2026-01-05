from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Optional, dict, List
import sympy as sp

@dataclass
class Finding:
    level: str
    message: str
    cause: Optional[str] = None
    fix: Optional[str] = None
    data: Optional[dict] = None

@dataclass
class PatternHit:
    name: str
    confidence: float
    note: str
    data: dict

x, y, z, t = sp.symbols("x y z t")
nu = sp.Symbol("nu", nonnegative=True)

def _as_matrix(v: Any) -> sp.Matrix:
    if isinstance(v, sp.Matrix):
        return v
    raise TypeError(f"Expected sympy.Matrix, got {type(v)}")

def is_zero(expr: Any) -> bool:
    try:
        expr = sp.together(expr)
        expr = sp.cancel(expr)
        expr = sp.simplify(expr)
        return expr == 0
    except Exception:
        return False

def div(u: sp.Matrix, coords: List[sp.Symbol]) -> sp.Expr:
    return sum(sp.diff(u[i], coords[i]) for i in range(len(coords)))

def grad(p: sp.Expr, coords: List[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix([sp.diff(p, c) for c in coords])

def lap_scalar(f: sp.Expr, coords: List[sp.Symbol]) -> sp.Expr:
    return sum(sp.diff(f, c, 2) for c in coords)

def lap_vec(u: sp.Matrix, coords: List[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix([lap_scalar(u[i], coords) for i in range(u.rows)])

def dt_vec(u: sp.Matrix, t: sp.Symbol) -> sp.Matrix:
    return sp.Matrix([sp.diff(u[i], t) for i in range(u.rows)])

def convection(u: sp.Matrix, coords: List[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix([
        sum(u[j] * sp.diff(u[i], coords[j]) for j in range(len(coords)))
        for i in range(u.rows)
    ])

def curl3d(u: sp.Matrix, coords: List[sp.Symbol]) -> sp.Matrix:
    x, y, z = coords
    return sp.Matrix([
        sp.diff(u[2], y) - sp.diff(u[1], z),
        sp.diff(u[0], z) - sp.diff(u[2], x),
        sp.diff(u[1], x) - sp.diff(u[0], y),
    ])

def curl2d_scalar(u: sp.Matrix, coords: List[sp.Symbol]) -> sp.Expr:
    x, y = coords
    return sp.diff(u[1], x) - sp.diff(u[0], y)

def divergence_of_vector(F: sp.Matrix, coords: List[sp.Symbol]) -> sp.Expr:
    return sum(sp.diff(F[i], coords[i]) for i in range(len(coords)))

def poisson_pressure_residual(u, p, f, coords, t, steady):
    ut = sp.zeros(u.rows, 1) if steady else dt_vec(u, t)
    return lap_scalar(p, coords) + divergence_of_vector(ut + convection(u, coords) - f, coords)

def is_conservative_2d(F, coords):
    curlF = sp.diff(F[1], coords[0]) - sp.diff(F[0], coords[1])
    return is_zero(curlF), curlF

def integrate_grad_2d(F, coords):
    x, y = coords
    try:
        p_x = sp.integrate(F[0], x)
        delta = F[1] - sp.diff(p_x, y)
        p = p_x + sp.integrate(delta, y)
    except Exception:
        return None, "Integration failure"
    if is_zero(sp.diff(p, x) - F[0]) and is_zero(sp.diff(p, y) - F[1]):
        return p, "OK"
    return None, "Gradient mismatch"

def detect_patterns(u, coords):
    hits = []
    if len(coords) == 2:
        if is_zero(curl2d_scalar(u, coords)):
            hits.append(PatternHit("PotentialFlow2D", 0.6, "ω=0", {}))
    if len(coords) == 3:
        w = curl3d(u, coords)
        if all(is_zero(w[i]) for i in range(3)):
            hits.append(PatternHit("PotentialFlow3D", 0.6, "curl(u)=0", {}))
    return hits

class GeRoMV61:
    def __init__(self):
        self.findings = []
        self.patterns = []

    def run(self, ctx):
        self.findings.clear()
        self.patterns.clear()

        u = _as_matrix(ctx["u"])
        coords = ctx["coords"]
        f = ctx.get("f", sp.zeros(u.rows, 1))
        nu = ctx["nu"]
        t = ctx["t"]
        ass = ctx.get("assumptions", {})
        steady = bool(ass.get("steady", False))

        if u.rows != len(coords):
            return self._fatal("Dimension mismatch u/coords")

        if f.rows != u.rows:
            return self._fatal("Dimension mismatch f/u")

        self.patterns = detect_patterns(u, coords)

        if not is_zero(div(u, coords)):
            return self._fatal("div(u) != 0")

        if ctx.get("p") is None and len(coords) == 2:
            gp = - (sp.zeros(u.rows,1) if steady else dt_vec(u,t)) \
                 - convection(u, coords) + nu * lap_vec(u, coords) + f

            ok, curl = is_conservative_2d(gp, coords)
            if not ok:
                return self._fatal("grad(p) non conservative")

            p, _ = integrate_grad_2d(gp, coords)
            if p is None:
                return self._fatal("Pressure integration failed")

            ctx["p"] = sp.simplify(p)

        p = ctx.get("p")
        res = (sp.zeros(u.rows,1) if steady else dt_vec(u,t)) \
              + convection(u, coords) + grad(p, coords) - nu * lap_vec(u, coords) - f

        if not all(is_zero(res[i]) for i in range(res.rows)):
            return self._fatal("Momentum residual non-zero")

        if isinstance(nu, sp.Expr) and nu.free_symbols & set(coords):
            self.findings.append(Finding("WARN", "nu(x) detected: Poisson check indicative only"))

        else:
            r = poisson_pressure_residual(u, p, f, coords, t, steady)
            if not is_zero(r):
                self.findings.append(Finding("WARN", "Poisson residual non-zero"))

        if ass.get("energy_inequality", False):
            self.findings.append(Finding("INFO", "Energy inequality declared"))
        else:
            self.findings.append(Finding("WARN", "Energy inequality not declared"))

        return self._report(ctx)

    def _fatal(self, msg):
        self.findings.append(Finding("ERROR", msg))
        return self._report({})

    def _report(self, ctx):
        return {
            "ok": not any(f.level == "ERROR" for f in self.findings),
            "findings": [f.__dict__ for f in self.findings],
            "patterns": [p.__dict__ for p in self.patterns],
            "context_echo": {
                "u": str(ctx.get("u")),
                "p": str(ctx.get("p")),
                "nu": str(ctx.get("nu")),
            }
        }

# Demo
if __name__ == "__main__":
    ctx = {
        "u": sp.Matrix([-y, x]),
        "nu": nu,
        "coords": [x, y],
        "t": t,
        "assumptions": {"steady": True, "energy_inequality": True},
        "f": sp.Matrix([0, 0])
    }
    gerom = GeRoMV61()
    report = gerom.run(ctx)
    print(json.dumps(report, default=str, indent=2))
