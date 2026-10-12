"""Generate reference results from calfem-python for the calfem-cpp tests.

Every case calls a calfem-python function with fixed (seeded) inputs and stores the inputs
and the returned values in calfem_python_reference.json. test_calfem_python.cpp runs the
same inputs through the C++ functions and compares the results.

Usage (needs numpy and scipy):

    python generate_calfem_python.py <path to calfem-python/src> [output.json]

Optional arguments (eq, nep, ...) are only stored when they were passed to the Python
function, so the C++ test can call the matching overload.
"""

import contextlib
import io
import json
import os
import subprocess
import sys

import numpy as np


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    src = os.path.abspath(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "calfem_python_reference.json")

    sys.path.insert(0, src)
    import calfem
    import calfem.core as cfc
    import calfem.utils as cfu

    gen = Generator(cfc, cfu)
    gen.generate_all()

    try:
        commit = subprocess.run(["git", "-C", src, "rev-parse", "HEAD"], capture_output=True, text=True,
                                check=True).stdout.strip()
    except Exception:
        commit = "unknown"

    data = {
        "generator": "tests/reference_data/generate_calfem_python.py",
        "calfem_python_version": calfem.__version__,
        "calfem_python_commit": commit,
        "numpy_version": np.__version__,
        "cases": gen.cases,
    }

    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)
        f.write("\n")

    counts = {}
    for case in gen.cases:
        counts[case["function"]] = counts.get(case["function"], 0) + 1

    print(f"Wrote {len(gen.cases)} cases to {out}")
    for name, count in counts.items():
        print(f"  {name}: {count}")
    for name, reason in gen.skipped:
        print(f"  skipped {name}: {reason}")


def to_json(value):
    """Converts a Python result to a 2D list (scalars become [[x]], column vectors [[a], [b]])."""
    a = np.asarray(value, dtype=float)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(-1, 1)
    return a.tolist()


class Generator:
    def __init__(self, cfc, cfu=None):
        self.cfc = cfc
        self.cfu = cfu
        self.rng = np.random.default_rng(20261012)
        self.cases = []
        self.skipped = []

    # ---- helpers ----

    def call(self, func, inputs, order, outputs, array_inputs=("ed", "ev", "ea", "eo", "D")):
        """Calls cfc.<func> with the inputs in the given order and records the case."""
        args = []
        for name in order:
            value = inputs.get(name)
            if name in array_inputs and value is not None and not isinstance(value, np.ndarray):
                value = np.array(value, dtype=float)
            args.append(value)

        # Drop trailing None arguments so that defaults are used, as a user would call it
        while args and args[-1] is None:
            args.pop()

        try:
            with contextlib.redirect_stdout(io.StringIO()):
                result = getattr(self.cfc, func)(*args)
        except Exception as e:  # e.g. np.arange giving the wrong number of points
            self.skipped.append((func, f"{type(e).__name__}: {e}"))
            return

        if not isinstance(result, tuple):
            result = (result,)

        if len(result) != len(outputs):
            raise RuntimeError(f"{func}: expected {len(outputs)} outputs, got {len(result)}")

        self.cases.append({
            "function": func,
            "inputs": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in inputs.items()
                       if v is not None},
            "outputs": [{"name": n, "value": to_json(r)} for n, r in zip(outputs, result)],
        })

    def uniform(self, lo, hi, n=None):
        v = self.rng.uniform(lo, hi, n)
        return v.tolist() if n is not None else float(v)

    def ex1(self):
        x1 = self.uniform(-3.0, 3.0)
        return [x1, x1 + self.uniform(1.0, 5.0)]

    def exy2(self):
        x1, y1 = self.uniform(-3.0, 3.0), self.uniform(-3.0, 3.0)
        L, angle = self.uniform(1.0, 5.0), self.uniform(-np.pi, np.pi)
        return [x1, x1 + L * np.cos(angle)], [y1, y1 + L * np.sin(angle)]

    def disp(self, n, scale=1e-3):
        return self.uniform(-scale, scale, n)

    # Variants of the optional arguments: (with eq, nep)
    VARIANTS = [(False, None), (True, None), (False, 4), (True, 2), (True, 5), (True, 7)]

    # ---- elements ----

    def generate_all(self):
        self.springs()
        self.bars_1d()
        self.bars_2d()
        self.beams_1d()
        self.beams_2d()
        self.beams_2d_second_order()
        self.beams_2d_dynamic()
        self.elements_3d()
        self.solvers()
        self.time_integration()
        self.field_elements()
        self.continuum()
        self.plates()
        self.utilities()

    def springs(self):
        for _ in range(3):
            k = self.uniform(1e3, 1e6)
            self.call("spring1e", {"ep": k}, ["ep"], ["Ke"])
            self.call("spring1s", {"ep": k, "ed": self.disp(2)}, ["ep", "ed"], ["es"])

    def stiffness_and_section(self, e_func, s_func, make_inputs, e_order, s_order, load):
        for with_eq, nep in self.VARIANTS:
            inputs = make_inputs()
            eq = load() if with_eq else None

            e_inputs = {k: inputs[k] for k in e_order if k in inputs}
            e_inputs["eq"] = eq
            self.call(e_func, e_inputs, e_order + ["eq"], ["Ke", "fe"] if with_eq else ["Ke"])

            s_inputs = {k: inputs[k] for k in s_order if k in inputs}
            s_inputs["eq"] = eq
            s_inputs["nep"] = nep
            outputs = ["es", "edi", "eci"] if nep is not None else ["es"]
            self.call(s_func, s_inputs, s_order + ["eq", "nep"], outputs)

    def bars_1d(self):
        self.stiffness_and_section(
            "bar1e", "bar1s",
            lambda: {"ex": self.ex1(), "ep": [self.uniform(1e9, 3e11), self.uniform(1e-4, 1e-2)], "ed": self.disp(2)},
            ["ex", "ep"], ["ex", "ep", "ed"], lambda: [self.uniform(-1e4, 1e4)])

        self.stiffness_and_section(
            "bar1we", "bar1ws",
            lambda: {"ex": self.ex1(), "ep": [self.uniform(1e9, 3e11), self.uniform(1e-4, 1e-2),
                                              self.uniform(1e6, 1e8)], "ed": self.disp(2)},
            ["ex", "ep"], ["ex", "ep", "ed"], lambda: [self.uniform(-1e4, 1e4)])

    def bars_2d(self):
        def inputs():
            ex, ey = self.exy2()
            return {"ex": ex, "ey": ey, "ep": [self.uniform(1e9, 3e11), self.uniform(1e-4, 1e-2)], "ed": self.disp(4)}

        self.stiffness_and_section("bar2e", "bar2s", inputs, ["ex", "ey", "ep"], ["ex", "ey", "ep", "ed"],
                                   lambda: [self.uniform(-1e4, 1e4)])

        for nep in [None, None, 3, 6]:
            d = inputs()
            QX = self.uniform(-1e5, 1e5)
            self.call("bar2ge", {"ex": d["ex"], "ey": d["ey"], "ep": d["ep"], "QX": QX}, ["ex", "ey", "ep", "QX"],
                      ["Ke"])
            self.call("bar2gs", {"ex": d["ex"], "ey": d["ey"], "ep": d["ep"], "ed": d["ed"], "nep": nep},
                      ["ex", "ey", "ep", "ed", "nep"], ["es", "QX", "edi", "eci"] if nep else ["es", "QX"])

    def beams_1d(self):
        def ed():
            return self.disp(4)

        self.stiffness_and_section(
            "beam1e", "beam1s",
            lambda: {"ex": self.ex1(), "ep": [self.uniform(1e9, 3e11), self.uniform(1e-6, 1e-4)], "ed": ed()},
            ["ex", "ep"], ["ex", "ep", "ed"], lambda: [self.uniform(-1e4, 1e4)])

        self.stiffness_and_section(
            "beam1we", "beam1ws",
            lambda: {"ex": self.ex1(), "ep": [self.uniform(1e9, 3e11), self.uniform(1e-6, 1e-4),
                                              self.uniform(1e5, 1e7)], "ed": ed()},
            ["ex", "ep"], ["ex", "ep", "ed"], lambda: [self.uniform(-1e4, 1e4)])

    def beams_2d(self):
        def load():
            return [self.uniform(-1e4, 1e4), self.uniform(-1e4, 1e4)]

        def inputs(ep):
            def make():
                ex, ey = self.exy2()
                return {"ex": ex, "ey": ey, "ep": ep(), "ed": self.disp(6)}
            return make

        order_e = ["ex", "ey", "ep"]
        order_s = ["ex", "ey", "ep", "ed"]

        self.stiffness_and_section(
            "beam2e", "beam2s",
            inputs(lambda: [self.uniform(1e9, 3e11), self.uniform(1e-3, 1e-1), self.uniform(1e-6, 1e-4)]),
            order_e, order_s, load)

        self.stiffness_and_section(
            "beam2we", "beam2ws",
            inputs(lambda: [self.uniform(1e9, 3e11), self.uniform(1e-3, 1e-1), self.uniform(1e-6, 1e-4),
                            self.uniform(1e5, 1e7), self.uniform(1e5, 1e7)]),
            order_e, order_s, load)

        self.stiffness_and_section(
            "beam2te", "beam2ts",
            inputs(lambda: [self.uniform(1e9, 3e11), self.uniform(5e8, 1e11), self.uniform(1e-3, 1e-1),
                            self.uniform(1e-6, 1e-4), self.uniform(0.5, 1.0)]),
            order_e, order_s, load)

    def beams_2d_second_order(self):
        # Compression, tension and (nearly) zero axial force
        for e_func, s_func in [("beam2ge", "beam2gs"), ("beam2gxe", "beam2gxs")]:
            for QX_range in [(-2e5, -1e3), (1e3, 2e5), (0.0, 0.0)]:
                for with_eq, nep in self.VARIANTS:
                    ex, ey = self.exy2()
                    ep = [self.uniform(1e10, 3e11), self.uniform(1e-3, 1e-1), self.uniform(1e-5, 1e-4)]
                    QX = self.uniform(*QX_range) if QX_range[0] != QX_range[1] else 0.0
                    eq = [self.uniform(-1e4, 1e4)] if with_eq else None

                    self.call(e_func, {"ex": ex, "ey": ey, "ep": ep, "QX": QX, "eq": eq},
                              ["ex", "ey", "ep", "QX", "eq"], ["Ke", "fe"] if with_eq else ["Ke"])

                    self.call(s_func, {"ex": ex, "ey": ey, "ep": ep, "ed": self.disp(6), "QX": QX, "eq": eq,
                                       "nep": nep},
                              ["ex", "ey", "ep", "ed", "QX", "eq", "nep"],
                              ["es", "QX", "edi", "eci"] if nep is not None else ["es", "QX"])

    def beams_2d_dynamic(self):
        for damping in [False, True, False, True]:
            ex, ey = self.exy2()
            ep = [self.uniform(1e9, 3e11), self.uniform(1e-3, 1e-1), self.uniform(1e-6, 1e-4),
                  self.uniform(10.0, 500.0)]
            if damping:
                ep += [self.uniform(0.0, 0.1), self.uniform(0.0, 1e-3)]

            self.call("beam2de", {"ex": ex, "ey": ey, "ep": ep}, ["ex", "ey", "ep"],
                      ["Ke", "Me", "Ce"] if damping else ["Ke", "Me"])

            nt = 3
            self.call("beam2ds", {"ex": ex, "ey": ey, "ep": ep,
                                  "ed": [self.disp(6) for _ in range(nt)],
                                  "ev": [self.disp(6, 1e-2) for _ in range(nt)],
                                  "ea": [self.disp(6, 1e-1) for _ in range(nt)]},
                      ["ex", "ey", "ep", "ed", "ev", "ea"], ["es"])

    def elements_3d(self):
        # The 3D elements predate this harness (they come from ObjectiveFrame); compare them too.
        def ez3():
            x1, y1, z1 = self.uniform(-3, 3, 3)
            d = self.rng.normal(size=3)
            d = d / np.linalg.norm(d) * self.uniform(1.0, 5.0)
            return [x1, x1 + d[0]], [y1, y1 + d[1]], [z1, z1 + d[2]], d

        for with_eq, nep in self.VARIANTS:
            ex, ey, ez, d = ez3()
            ep = [self.uniform(1e9, 3e11), self.uniform(1e-4, 1e-2)]
            eq = [self.uniform(-1e4, 1e4)] if with_eq else None
            self.call("bar3e", {"ex": ex, "ey": ey, "ez": ez, "ep": ep, "eq": eq}, ["ex", "ey", "ez", "ep", "eq"],
                      ["Ke", "fe"] if with_eq else ["Ke"])
            self.call("bar3s", {"ex": ex, "ey": ey, "ez": ez, "ep": ep, "ed": self.disp(6), "eq": eq, "nep": nep},
                      ["ex", "ey", "ez", "ep", "ed", "eq", "nep"],
                      ["es", "edi", "eci"] if nep is not None else ["es"])

        for with_eq, nep in self.VARIANTS:
            ex, ey, ez, d = ez3()
            eo = np.cross(d, self.rng.normal(size=3)).tolist()
            ep = [self.uniform(1e9, 3e11), self.uniform(5e8, 1e11), self.uniform(1e-3, 1e-1),
                  self.uniform(1e-6, 1e-4), self.uniform(1e-6, 1e-4), self.uniform(1e-6, 1e-4)]
            eq = self.uniform(-1e4, 1e4, 4) if with_eq else None
            self.call("beam3e", {"ex": ex, "ey": ey, "ez": ez, "eo": eo, "ep": ep, "eq": eq},
                      ["ex", "ey", "ez", "eo", "ep", "eq"], ["Ke", "fe"] if with_eq else ["Ke"])
            self.call("beam3s", {"ex": ex, "ey": ey, "ez": ez, "eo": eo, "ep": ep, "ed": self.disp(12), "eq": eq,
                                 "nep": nep},
                      ["ex", "ey", "ez", "eo", "ep", "ed", "eq", "nep"],
                      ["es", "edi", "eci"] if nep is not None else ["es"])


    # ---- solvers ----

    def spd(self, n, scale=1.0):
        """Random symmetric positive definite matrix with distinct eigenvalues."""
        Q, _ = np.linalg.qr(self.rng.normal(size=(n, n)))
        d = np.sort(self.rng.uniform(1.0, 10.0, n)) + np.arange(n)
        return scale * (Q @ np.diag(d) @ Q.T)

    def solvers(self):
        A = self.rng.normal(size=(6, 6))
        self.call("red", {"A": A, "b": np.array([5, 2])}, ["A", "b"], ["B"])
        self.call("red", {"A": A, "bc": np.array([[2, 0], [6, 0]])}, ["A", "bc"], ["B"])

        for cd in [np.array([5, 2]), np.array([1]), np.array([6, 3, 4])]:
            self.call("statcon", {"K": self.spd(6, 1e3), "f": self.rng.normal(size=6), "cd": cd},
                      ["K", "f", "cd"], ["K1", "f1"])

        for b in [None, np.array([1, 4]), np.array([6, 2, 3])]:
            self.call("eigen", {"K": self.spd(6, 1e4), "M": self.spd(6), "b": b}, ["K", "M", "b"], ["L", "X"])

        for G, dt in [(np.array([[0.0, 0.0], [0.5, 1.0], [1.5, -1.0], [2.0, 0.0]]), 0.25),
                      (np.array([[0.0, 1.0], [1.0, 3.0], [4.0, 2.0]]), 0.5),
                      (np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 1.0]]), 0.3)]:
            ti, _ = self.cfc.gfunc(G, dt)
            expected = int(np.floor((G[-1, 0] - G[0, 0]) / dt + 1e-9)) + 1
            if len(ti) != expected:
                self.skipped.append(("gfunc", f"np.arange gave {len(ti)} points instead of {expected}"))
                continue
            self.call("gfunc", {"G": G, "dt": dt}, ["G", "dt"], ["ti", "g"])

    def time_integration(self):
        n = 4
        K = np.array([[2.0, -1.0, 0.0, 0.0], [-1.0, 2.0, -1.0, 0.0], [0.0, -1.0, 2.0, -1.0],
                      [0.0, 0.0, -1.0, 1.0]]) * 100.0
        C = np.diag([1.0, 2.0, 1.5, 0.5])
        M = np.diag([2.0, 1.0, 3.0, 1.5])
        a0 = np.array([[0.0], [0.01], [-0.02], [0.03]])
        da0 = np.array([[0.1], [0.0], [-0.1], [0.2]])
        empty = np.array([])

        dt, T = 0.05, 1.0
        nstep = int(T / dt)
        f_const = np.array([[0.0], [1.0], [0.0], [-2.0]])
        f_hist = np.vstack([np.zeros(nstep + 1), np.sin(np.linspace(0, 3, nstep + 1)),
                            np.zeros(nstep + 1), np.linspace(0, 1, nstep + 1)])
        bc_fixed = np.array([[1, 0.0]])
        bc_hist = np.hstack([np.array([[1.0]]), 0.01 * np.sin(np.linspace(0, 2, nstep + 1)).reshape(1, -1)])
        times = np.array([0.0, 0.3, 0.6, 0.9])
        dofs = np.array([2, 4])

        # (f, bc, times, dofs). calfem-python records dofhist only when times are given.
        variants = [(f_const, empty, empty, dofs), (f_hist, bc_fixed, times, dofs), (empty, bc_hist, times, dofs),
                    (f_const, bc_fixed, empty, empty), (f_hist, empty, times, dofs)]

        def outputs(names, times, dofs):
            out = ["modelhist_" + k for k in names]
            if times.size and dofs.size:
                out += ["dofhist_" + k for k in names]
            return out

        def unpack(func, *args):
            modelhist, dofhist = getattr(self.cfc, func)(*args)
            keys = ["a", "da"] + (["d2a"] if func == "step2" else [])
            result = tuple(modelhist[k] for k in keys)
            if dofhist:
                result += tuple(dofhist[k] for k in keys)
            return result

        for alpha in [0.5, 1.0]:
            for f, bc, t, d in variants:
                inputs = {"K": K, "C": C, "f": f, "a0": a0, "bc": bc, "ip": [dt, T, alpha], "times": t, "dofs": d}
                order = ["K", "C", "f", "a0", "bc", "ip", "times", "dofs"]
                self.call_custom("step1", inputs, order, outputs(["a", "da"], t, d), unpack)

        for damping in [True, False]:
            for f, bc, t, d in variants:
                inputs = {"K": K, "C": C if damping else empty, "M": M, "f": f, "a0": a0, "da0": da0, "bc": bc,
                          "ip": [dt, T, 0.25, 0.5], "times": t, "dofs": d}
                order = ["K", "C", "M", "f", "a0", "da0", "bc", "ip", "times", "dofs"]
                self.call_custom("step2", inputs, order, outputs(["a", "da", "d2a"], t, d), unpack)

    # ---- field problems ----

    def quad(self, x0=0.0, y0=0.0):
        """Convex quadrilateral: perturbed rectangle, counter-clockwise."""
        w, h = self.uniform(0.5, 2.0), self.uniform(0.5, 2.0)
        ex = np.array([0, w, w, 0]) + x0 + self.rng.uniform(-0.1, 0.1, 4) * w
        ey = np.array([0, 0, h, h]) + y0 + self.rng.uniform(-0.1, 0.1, 4) * h
        return ex.tolist(), ey.tolist()

    def quad8(self):
        ex, ey = self.quad()
        mx = [(ex[i] + ex[(i + 1) % 4]) / 2 for i in range(4)]
        my = [(ey[i] + ey[(i + 1) % 4]) / 2 for i in range(4)]
        mx = (np.array(mx) + self.rng.uniform(-0.03, 0.03, 4)).tolist()
        my = (np.array(my) + self.rng.uniform(-0.03, 0.03, 4)).tolist()
        return ex + mx, ey + my

    def brick(self):
        a, b, c = self.uniform(0.5, 2.0), self.uniform(0.5, 2.0), self.uniform(0.5, 2.0)
        ex = np.array([0, a, a, 0, 0, a, a, 0]) + self.rng.uniform(-0.05, 0.05, 8)
        ey = np.array([0, 0, b, b, 0, 0, b, b]) + self.rng.uniform(-0.05, 0.05, 8)
        ez = np.array([0, 0, 0, 0, c, c, c, c]) + self.rng.uniform(-0.05, 0.05, 8)
        return ex.tolist(), ey.tolist(), ez.tolist()

    def triangle(self):
        ex = [0.0, self.uniform(1.0, 2.0), self.uniform(-0.5, 0.5)]
        ey = [0.0, self.uniform(-0.3, 0.3), self.uniform(1.0, 2.0)]
        return ex, ey

    def field_elements(self):
        def at2(func, *args):
            return tuple(np.atleast_2d(r) for r in getattr(self.cfc, func)(*args))

        for with_eq in [False, True, False, True]:
            ex, ey = self.triangle()
            D, ep = self.spd(2), [self.uniform(0.1, 1.0)]
            eq = self.uniform(1.0, 100.0) if with_eq else None
            self.call("flw2te", {"ex": ex, "ey": ey, "ep": ep, "D": D, "eq": eq}, ["ex", "ey", "ep", "D", "eq"],
                      ["Ke", "fe"] if with_eq else ["Ke"])

        for nel in [1, 3]:
            tris = [self.triangle() for _ in range(nel)]
            ex = np.array([t[0] for t in tris])
            ey = np.array([t[1] for t in tris])
            inputs = {"ex": ex, "ey": ey, "D": self.spd(2), "ed": self.rng.normal(size=(nel, 3))}
            self.call_custom("flw2ts", inputs, ["ex", "ey", "D", "ed"], ["es", "et"], at2)

        for with_eq in [False, True, False, True]:
            ex, ey = self.quad()
            D, ep = self.spd(2), [self.uniform(0.1, 1.0)]
            eq = self.uniform(1.0, 100.0) if with_eq else None
            self.call("flw2qe", {"ex": ex, "ey": ey, "ep": ep, "D": D, "eq": eq}, ["ex", "ey", "ep", "D", "eq"],
                      ["Ke", "fe"] if with_eq else ["Ke"])
            # calfem-python's flw2qs fails for more than one row in ed and returns 1D arrays for one row
            inputs = {"ex": ex, "ey": ey, "ep": ep, "D": D, "ed": self.rng.normal(size=(1, 4))}
            order = ["ex", "ey", "ep", "D", "ed"]
            if with_eq:
                inputs["eq"] = eq
                order.append("eq")
            self.call_custom("flw2qs", inputs, order, ["es", "et"], at2)

        for e_func, s_func, geometry, nnodes, rules in [("flw2i4e", "flw2i4s", self.quad, 4, [1, 2, 3]),
                                                       ("flw2i8e", "flw2i8s", self.quad8, 8, [2, 3])]:
            for ir in rules:
                for with_eq in [False, True]:
                    ex, ey = geometry()
                    D, ep = self.spd(2), [self.uniform(0.1, 1.0), ir]
                    eq = self.uniform(1.0, 100.0) if with_eq else None
                    self.call(e_func, {"ex": ex, "ey": ey, "ep": ep, "D": D, "eq": eq},
                              ["ex", "ey", "ep", "D", "eq"], ["Ke", "fe"] if with_eq else ["Ke"])
                    self.call(s_func, {"ex": ex, "ey": ey, "ep": ep, "D": D,
                                       "ed": self.rng.normal(size=(1 if with_eq else 2, nnodes))},
                              ["ex", "ey", "ep", "D", "ed"], ["es", "et", "eci"])

        # calfem-python's ir = 3 rule for flw3i8 is incorrect (see flow.h), only ir = 2 is compared
        for with_eq in [False, True, False]:
            ex, ey, ez = self.brick()
            D, ep = self.spd(3), [2]
            eq = self.uniform(1.0, 100.0) if with_eq else None
            self.call("flw3i8e", {"ex": ex, "ey": ey, "ez": ez, "ep": ep, "D": D, "eq": eq},
                      ["ex", "ey", "ez", "ep", "D", "eq"], ["Ke", "fe"] if with_eq else ["Ke"])
            self.call("flw3i8s", {"ex": ex, "ey": ey, "ez": ez, "ep": ep, "D": D,
                                  "ed": self.rng.normal(size=(2, 8))},
                      ["ex", "ey", "ez", "ep", "D", "ed"], ["es", "et", "eci"])

    # ---- continuum ----

    def continuum(self):
        cfc = self.cfc

        def material(ptype):
            return np.asarray(cfc.hooke(ptype, self.uniform(1e9, 3e11), self.uniform(0.1, 0.45)))

        def rows(func, *args):
            return tuple(np.atleast_2d(r) for r in getattr(cfc, func)(*args))

        for ptype in [1, 2, 3, 4]:
            self.call("hooke", {"ptype": ptype, "E": self.uniform(1e9, 3e11), "v": self.uniform(0.0, 0.45)},
                      ["ptype", "E", "v"], ["D"])

        for ptype, ncomps in [(1, 3), (1, 4), (2, 4), (3, 4), (4, 6)]:
            self.call("effmises", {"es": self.rng.normal(size=(5, ncomps)) * 1e6, "ptype": ptype}, ["es", "ptype"],
                      ["eseff"])

        # Four 4-node elements in a 2 x 2 mesh, 2 dofs per node
        nodes = np.arange(1, 10).reshape(3, 3)
        topo = [[nodes[i, j], nodes[i, j + 1], nodes[i + 1, j + 1], nodes[i + 1, j]] for i in range(2) for j in range(2)]
        edof = np.array([[2 * n - 1 if k % 2 == 0 else 2 * n for n in t for k in range(2)] for t in topo])
        self.call("stress2nodal", {"eseff": self.rng.uniform(1e6, 1e8, 4), "edof": edof}, ["eseff", "edof"], ["ev"])

        # Plane elements: ptype 1/2 with 3 x 3 and 4 x 4 D
        for ptype, dtype in [(1, 1), (1, 2), (2, 2), (2, 1)]:
            for with_eq in [False, True]:
                D = material(dtype)
                t = self.uniform(0.01, 0.2)
                eq = self.uniform(-1e4, 1e4, 2) if with_eq else None
                ex, ey = self.triangle()

                inputs = {"ex": ex, "ey": ey, "ep": [ptype, t], "D": D}
                order = ["ex", "ey", "ep", "D"]
                if with_eq:
                    inputs["eq"] = eq
                    order.append("eq")
                # calfem-python returns fe as a row
                self.call_custom("plante", inputs, order, ["Ke", "fe"] if with_eq else ["Ke"],
                                 lambda f, *a: (lambda r: (r,) if not isinstance(r, tuple) else
                                                (r[0], np.asarray(r[1]).reshape(-1, 1)))(getattr(cfc, f)(*a)))

                nel = 3 if with_eq else 1
                tris = [self.triangle() for _ in range(nel)]
                self.call_custom("plants", {"ex": np.array([tr[0] for tr in tris]), "ey": np.array([tr[1] for tr in tris]),
                                            "ep": [ptype, t], "D": D, "ed": self.rng.normal(size=(nel, 6)) * 1e-3},
                                 ["ex", "ey", "ep", "D", "ed"], ["es", "et"], rows)

                self.call("plantf", {"ex": ex, "ey": ey, "ep": [ptype, t],
                                     "es": self.rng.normal(size=D.shape[0] if D.shape[0] < 6 else 3) * 1e6},
                          ["ex", "ey", "ep", "es"], ["ef"], array_inputs=("es",))

                qx, qy = self.quad()
                inputs = {"ex": qx, "ey": qy, "ep": [ptype, t], "D": D}
                order = ["ex", "ey", "ep", "D"]
                if with_eq:
                    inputs["eq"] = eq
                    order.append("eq")
                self.call("planqe", inputs, order, ["Ke", "fe"] if with_eq else ["Ke"])

                inputs = {"ex": qx, "ey": qy, "ep": [ptype, t], "D": D, "ed": self.rng.normal(size=8) * 1e-3}
                order = ["ex", "ey", "ep", "D", "ed"]
                if with_eq:
                    inputs["eq"] = eq
                    order.append("eq")
                self.call_custom("planqs", inputs, order, ["es", "et"], rows)

                for ir in [1, 2, 3]:
                    inputs = {"ex": qx, "ey": qy, "ep": [ptype, t, ir], "D": D}
                    order = ["ex", "ey", "ep", "D"]
                    if with_eq:
                        inputs["eq"] = eq
                        order.append("eq")
                    # calfem-python always returns (Ke, fe); without eq only Ke is compared
                    self.call_custom("plani4e", inputs, order, ["Ke", "fe"] if with_eq else ["Ke"],
                                     lambda f, *a: getattr(cfc, f)(*a))
                    self.call_custom("plani4s", {"ex": qx, "ey": qy, "ep": [ptype, t, ir], "D": D,
                                                 "ed": self.rng.normal(size=8) * 1e-3},
                                     ["ex", "ey", "ep", "D", "ed"], ["es", "et"], lambda f, *a: getattr(cfc, f)(*a))

        # Solid element
        for ir in [1, 2, 3]:
            for with_eq in [False, True]:
                ex, ey, ez = self.brick()
                D = material(4)
                inputs = {"ex": ex, "ey": ey, "ez": ez, "ep": [ir], "D": D}
                order = ["ex", "ey", "ez", "ep", "D"]
                if with_eq:
                    inputs["eq"] = self.uniform(-1e4, 1e4, 3)
                    order.append("eq")
                self.call("soli8e", inputs, order, ["Ke", "fe"] if with_eq else ["Ke"])
                self.call("soli8s", {"ex": ex, "ey": ey, "ez": ez, "ep": [ir], "D": D,
                                     "ed": self.rng.normal(size=24) * 1e-3},
                          ["ex", "ey", "ez", "ep", "D", "ed"], ["et", "es", "eci"])

    # ---- plates ----

    def plates(self):
        cfc = self.cfc
        for with_eq in [False, True, False, True]:
            x0, y0 = self.uniform(-2, 2), self.uniform(-2, 2)
            lx, ly = self.uniform(0.5, 2.0), self.uniform(0.5, 2.0)
            ex = np.array([x0, x0 + lx, x0 + lx, x0])
            ey = np.array([y0, y0, y0 + ly, y0 + ly])
            D = np.asarray(cfc.hooke(1, self.uniform(1e9, 3e11), self.uniform(0.1, 0.45)))
            inputs = {"ex": ex, "ey": ey, "ep": [self.uniform(0.01, 0.3)], "D": D}
            order = ["ex", "ey", "ep", "D"]
            if with_eq:
                inputs["eq"] = self.uniform(-1e4, 1e4)
                order.append("eq")
            # calfem-python returns fe as a row
            self.call_custom("platre", inputs, order, ["Ke", "fe"] if with_eq else ["Ke"],
                             lambda f, *a: (lambda r: (r,) if not isinstance(r, tuple) else
                                            (r[0], np.asarray(r[1]).reshape(-1, 1)))(getattr(cfc, f)(*a)))

    # ---- dof utilities, boundary conditions and loads ----

    def mesh(self, nx, ny, dofs_per_node, dim=2):
        """Structured mesh of 4-node elements: coords, dofs, edof (nodes numbered row by row)."""
        coords = []
        for j in range(ny + 1):
            for i in range(nx + 1):
                c = [i * 1.5 + self.uniform(-0.1, 0.1), j * 0.8 + self.uniform(-0.1, 0.1)]
                if dim == 3:
                    c.append(self.uniform(-0.5, 0.5))
                coords.append(c)
        coords = np.array(coords)
        dofs = self.cfc.create_dofs(len(coords), dofs_per_node)
        node = lambda i, j: j * (nx + 1) + i
        edof = np.array([np.hstack([dofs[node(i, j)], dofs[node(i + 1, j)], dofs[node(i + 1, j + 1)],
                                    dofs[node(i, j + 1)]]) for j in range(ny) for i in range(nx)])
        return coords, dofs, edof

    def utilities(self):
        cfc, cfu = self.cfc, self.cfu

        for n, d in [(4, 1), (5, 2), (3, 3)]:
            self.call("create_dofs", {"nCoords": n, "nDof": d}, ["nCoords", "nDof"], ["dofs"])

        for dpn, dim in [(1, 2), (2, 2), (3, 3)]:
            coords, dofs, edof = self.mesh(3, 2, dpn, dim)
            outputs = ["ex", "ey"] if dim == 2 else ["ex", "ey", "ez"]
            self.call("coord_extract", {"edof": edof, "coords": coords, "dofs": dofs},
                      ["edof", "coords", "dofs"], outputs, array_inputs=())

        # 2D frame: three beam elements
        ex = np.array([[0.0, 0.0], [0.0, 3.0], [3.0, 3.5]])
        ey = np.array([[0.0, 2.5], [2.5, 2.8], [2.8, 0.0]])
        ed = self.rng.normal(size=(3, 6)) * 1e-3
        for mag in [1.0, 50.0]:
            self.call("beam2crd", {"ex": ex, "ey": ey, "ed": ed, "mag": mag}, ["ex", "ey", "ed", "mag"],
                      ["excd", "eycd"], array_inputs=())

        coords, dofs, edof = self.mesh(3, 2, 2)
        ex2, ey2 = cfc.coord_extract(edof, coords, dofs)
        ed2 = self.rng.normal(size=edof.shape) * 1e-3
        for rat in [0.2, 0.1]:
            self.call_custom("scalfact2", {"ex": ex2, "ey": ey2, "ed": ed2, "rat": rat}, ["ex", "ey", "ed", "rat"],
                             ["sfac"], lambda f, *a: (getattr(cfu, f)(*a),))

        coords3, dofs3, edof3 = self.mesh(2, 2, 3, 3)
        ex3, ey3, ez3 = cfc.coord_extract(edof3, coords3, dofs3)
        for ignore_first, e in [(False, edof3), (True, np.hstack([np.arange(1, len(edof3) + 1).reshape(-1, 1),
                                                                    edof3]))]:
            self.call_custom("convert_to_node_topo",
                             {"edof": e, "ex": ex3, "ey": ey3, "ez": ez3, "n_dofs_per_node": 3,
                              "ignore_first": ignore_first},
                             ["edof", "ex", "ey", "ez", "n_dofs_per_node", "ignore_first"],
                             ["coords", "topo", "node_dofs"], lambda f, *a: getattr(cfu, f)(*a))

        # Boundary dofs as produced by the mesh generator: marker -> dofs (node by node)
        bdofs2 = {10: [1, 2, 3, 4, 5, 6], 20: [7, 8, 13, 14], 30: [9, 10]}
        bdofs3 = {1: [1, 2, 3, 4, 5, 6], 2: [10, 11, 12]}
        bc0, bcv0 = np.array([3, 30]), np.array([0.5, -0.2])

        def unpack_bc(f, *a):
            return getattr(cfu, f)(*a)

        for marker, value, dimension in [(10, 0.0, 0), (10, 1.5, 1), (10, -2.0, 2), (20, 0.3, 0)]:
            self.call_custom("apply_bc", {"bdofs": bdofs2, "bcPrescr": bc0, "bcVal": bcv0, "marker": marker,
                                          "value": value, "dimension": dimension},
                             ["bdofs", "bcPrescr", "bcVal", "marker", "value", "dimension"], ["bc", "bcVal"], unpack_bc)

        for marker, value, dimension in [(1, 0.0, 0), (1, 1.0, 3), (2, 2.0, 2)]:
            self.call_custom("apply_bc_3d", {"bdofs": bdofs3, "bcPrescr": bc0, "bcVal": bcv0, "marker": marker,
                                             "value": value, "dimension": dimension},
                             ["bdofs", "bcPrescr", "bcVal", "marker", "value", "dimension"], ["bc", "bcVal"], unpack_bc)

        dofs_n = cfc.create_dofs(4, 3)
        for node, value, dimension in [(0, 0.0, 0), (2, 1.5, 2), (3, -1.0, 3)]:
            self.call_custom("apply_bc_node", {"nodeIdx": node, "dofs": dofs_n, "bcPrescr": bc0, "bcVal": bcv0,
                                               "value": value, "dimension": dimension},
                             ["nodeIdx", "dofs", "bcPrescr", "bcVal", "value", "dimension"], ["bc", "bcVal"],
                             unpack_bc)

        # Assembly, static solvers and element displacements
        import scipy.sparse as sp

        for with_f in [False, True]:
            K = self.rng.normal(size=(6, 6))
            K = K + K.T
            f = self.rng.normal(size=(6, 1))
            Ke = self.rng.normal(size=(3, 3))
            fe = self.rng.normal(size=(3, 1))
            edof_e = np.array([5, 2, 3])
            inputs = {"edof": edof_e, "K": K, "Ke": Ke}
            order = ["edof", "K", "Ke"]
            if with_f:
                inputs.update({"f": f, "fe": fe})
                order += ["f", "fe"]
            self.call_custom("assem", inputs, order, ["K", "f"] if with_f else ["K"],
                             lambda fn, *a: (lambda r: r if isinstance(r, tuple) else (r,))(
                                 getattr(cfc, fn)(*[x.copy() if isinstance(x, np.ndarray) else x for x in a])))

        for bc_vals in [None, np.array([0.01, -0.02])]:
            K = self.spd(6, 1e3)
            f = self.rng.normal(size=(6, 1))
            inputs = {"K": K, "f": f, "bcPrescr": np.array([5, 2])}
            order = ["K", "f", "bcPrescr"]
            if bc_vals is not None:
                inputs["bcVal"] = bc_vals
                order.append("bcVal")
            self.call_custom("solveq", inputs, order, ["a", "r"], lambda fn, *a: getattr(cfc, fn)(*a))
            self.call_custom("spsolveq", inputs, order, ["a", "r"],
                             lambda fn, K, *a: getattr(cfc, fn)(sp.csr_matrix(K), *a))

        edof_d = np.array([[1, 2, 5, 6], [3, 4, 5, 6]])
        self.call_custom("extract_eldisp", {"edof": edof_d, "a": self.rng.normal(size=(6, 1))}, ["edof", "a"],
                         ["ed"], lambda fn, *a: (getattr(cfc, fn)(*a),))

        # Loads modify f in place: record f after the call
        def unpack_force(func, bdofs, f, *rest):
            f = f.copy()
            getattr(cfu, func)(bdofs, f, *rest)
            return (f,)

        f0 = self.rng.normal(size=30)
        for func, bdofs, cases in [("apply_force", bdofs2, [(10, 5.0, 0), (20, -3.0, 1), (10, 2.0, 2)]),
                                   ("apply_force_total", bdofs2, [(10, 6.0, 0), (20, -4.0, 2)]),
                                   ("apply_force_3d", bdofs3, [(1, 5.0, 0), (1, 2.0, 3)]),
                                   ("apply_force_total_3d", bdofs3, [(1, 9.0, 0), (2, 3.0, 1)])]:
            for marker, value, dimension in cases:
                self.call_custom(func, {"bdofs": bdofs, "f": f0, "marker": marker, "value": value,
                                        "dimension": dimension},
                                 ["bdofs", "f", "marker", "value", "dimension"], ["f"], unpack_force)

    def call_custom(self, func, inputs, order, outputs, unpack):
        """Like call(), for functions whose results need unpacking (dicts)."""
        args = [inputs[name] for name in order]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                result = unpack(func, *args)
        except Exception as e:
            self.skipped.append((func, f"{type(e).__name__}: {e}"))
            return

        # Only keep the outputs that are recorded (dofhist is buggy without times, see header)
        result = result[:len(outputs)]

        self.cases.append({
            "function": func,
            "inputs": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in inputs.items()},
            "outputs": [{"name": n, "value": to_json(r)} for n, r in zip(outputs, result)],
        })


if __name__ == "__main__":
    main()
