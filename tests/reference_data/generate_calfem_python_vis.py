"""Generate reference drawings from calfem-python's vis_mpl for the calfem-cpp tests.

Every case calls a calfem-python drawing function (eldraw2, eldisp2, ...) with fixed inputs on
an empty matplotlib figure (Agg backend) and records what was drawn: lines (Line2D and polygon
collections), node markers (scatter collections) and text, plus the returned scale factor.
test_vis_reference.cpp draws the same inputs with calfem::vis and compares the figures.

Usage (needs numpy, scipy and matplotlib):

    python generate_calfem_python_vis.py <path to calfem-python/src> [output.json]
"""

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
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__),
                                                             "calfem_python_vis_reference.json")

    import matplotlib
    matplotlib.use("Agg")

    sys.path.insert(0, src)
    import calfem
    import calfem.vis_mpl as cfv

    gen = Generator(cfv)
    gen.generate_all()

    try:
        commit = subprocess.run(["git", "-C", src, "rev-parse", "HEAD"], capture_output=True, text=True,
                                check=True).stdout.strip()
    except Exception:
        commit = "unknown"

    data = {
        "generator": "tests/reference_data/generate_calfem_python_vis.py",
        "calfem_python_version": calfem.__version__,
        "calfem_python_commit": commit,
        "matplotlib_version": matplotlib.__version__,
        "cases": gen.cases,
    }

    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)
        f.write("\n")

    counts = {}
    for case in gen.cases:
        counts[case["function"]] = counts.get(case["function"], 0) + 1
    print(f"Wrote {len(gen.cases)} cases to {out}: {counts}")


def to_json(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [to_json(v) for v in value]
    return value


def rgb(color):
    import matplotlib.colors as mcolors
    return [float(c) for c in mcolors.to_rgb(color)]


def line_style(style):
    """Line style name from a matplotlib linestyle or collection dash tuple."""
    if isinstance(style, tuple):
        offset, dashes = style
        if dashes is None or len(dashes) == 0:
            return "solid"
        return "dashed" if abs(dashes[0] - dashes[1]) < 1e-9 else "dotted"
    return {"-": "solid", "solid": "solid", "--": "dashed", "dashed": "dashed", ":": "dotted",
            "dotted": "dotted"}.get(style, str(style))


def marker_name(path):
    """Name of the scatter marker whose (transformed) path is path."""
    from matplotlib.markers import MarkerStyle
    for name, label in (("o", "circle"), ("x", "cross"), ("*", "star"), (".", "dot")):
        m = MarkerStyle(name)
        p = m.get_path().transformed(m.get_transform())
        if p.vertices.shape == path.vertices.shape and np.allclose(p.vertices, path.vertices):
            return label
    return "unknown"


def harvest(fig):
    """Lines, markers and texts drawn in the axes of fig."""
    import matplotlib.collections as mcoll

    lines, markers, texts = [], [], []
    for ax in fig.axes:
        for line in ax.lines:
            lines.append({
                "points": line.get_xydata().tolist(),
                "closed": False,
                "color": rgb(line.get_color()),
                "width": float(line.get_linewidth()),
                "style": line_style(line.get_linestyle()),
            })
        for coll in ax.collections:
            if isinstance(coll, mcoll.PathCollection):  # scatter
                edge = coll.get_edgecolor()
                marker = marker_name(coll.get_paths()[0])
                # matplotlib ignores the edge colour of the unfilled "x" marker, so calfem-python's
                # crosses are invisible (calfem::vis draws them in the line colour): no colour stored.
                color = None if marker == "cross" or len(edge) == 0 else rgb(edge[0])
                markers.append({
                    "points": np.asarray(coll.get_offsets()).tolist(),
                    "marker": marker,
                    "color": color,
                })
            elif isinstance(coll, mcoll.PolyCollection):
                edge = coll.get_edgecolor()
                for path in coll.get_paths():
                    verts = path.vertices
                    closed = path.codes is not None and path.codes[-1] == path.CLOSEPOLY
                    if closed:
                        verts = verts[:-1]
                    lines.append({
                        "points": verts.tolist(),
                        "closed": bool(closed),
                        "color": rgb(edge[0]),
                        "width": float(coll.get_linewidth()[0]),
                        "style": line_style(coll.get_linestyle()[0]),
                    })
        for text in ax.texts:
            x, y = text.get_position()
            texts.append({
                "position": [float(x), float(y)],
                "text": text.get_text(),
                "color": rgb(text.get_color()),
                "halign": text.get_horizontalalignment(),
                "valign": text.get_verticalalignment(),
            })
    return {"lines": lines, "markers": markers, "texts": texts}


class Generator:
    def __init__(self, cfv):
        self.cfv = cfv
        self.cases = []
        self.rng = np.random.default_rng(20260412)

    def call(self, name, inputs, **kwargs):
        """Draws with cfv.<name>(**inputs, **kwargs) on a new figure and stores the case.

        Optional arguments are only passed (and stored) when given in inputs.
        """
        import matplotlib.pyplot as plt

        func = getattr(self.cfv, name)
        fig = plt.figure()
        result = func(**inputs)
        drawing = harvest(fig)
        plt.close(fig)

        case = {"function": name, "inputs": {k: to_json(v) for k, v in inputs.items()}}
        if result is not None:
            case["sfac"] = float(result)
        case.update(drawing)
        case.update(kwargs)
        self.cases.append(case)

    # ---- Element coordinates ----

    def bars(self):
        ex = np.array([[0.0, 3.0], [3.0, 3.0], [0.0, 3.0]])
        ey = np.array([[0.0, 0.0], [0.0, 4.0], [0.0, 4.0]])
        return ex, ey

    def quads(self):
        # 2 x 2 distorted quadrilaterals
        x = np.array([0.0, 1.1, 2.0])
        y = np.array([0.0, 0.9, 2.0])
        ex, ey = [], []
        for j in range(2):
            for i in range(2):
                ex.append([x[i], x[i + 1], x[i + 1] + 0.1 * j, x[i] + 0.1 * j])
                ey.append([y[j], y[j] + 0.05 * i, y[j + 1], y[j + 1]])
        return np.array(ex), np.array(ey)

    def triangles(self):
        ex, ey = self.quads()
        tx = np.vstack([ex[:, [0, 1, 2]], ex[:, [0, 2, 3]]])
        ty = np.vstack([ey[:, [0, 1, 2]], ey[:, [0, 2, 3]]])
        return tx, ty

    def quads8(self):
        ex, ey = self.quads()
        mx = 0.5 * (ex + np.roll(ex, -1, axis=1))
        my = 0.5 * (ey + np.roll(ey, -1, axis=1))
        # Curved midside nodes
        mx = mx + self.rng.uniform(-0.05, 0.05, mx.shape)
        my = my + self.rng.uniform(-0.05, 0.05, my.shape)
        return np.hstack([ex, mx]), np.hstack([ey, my])

    def disp(self, rows, cols, scale=1e-3):
        return self.rng.uniform(-1.0, 1.0, (rows, cols)) * scale

    # ---- Cases ----

    def generate_all(self):
        self.eldraw2()
        self.eldisp2()
        self.dispbeam2()
        self.secforce2()
        self.scalgraph2()

    def eldraw2(self):
        ex, ey = self.bars()
        self.call("eldraw2", {"ex": ex, "ey": ey})
        self.call("eldraw2", {"ex": ex, "ey": ey, "plotpar": [2, 4, 2], "elnum": [1, 2, 3]})
        self.call("eldraw2", {"ex": ex[0], "ey": ey[0], "plotpar": [1, 1, 0]})

        tx, ty = self.triangles()
        self.call("eldraw2", {"ex": tx, "ey": ty, "plotpar": [3, 3, 0]})

        qx, qy = self.quads()
        self.call("eldraw2", {"ex": qx, "ey": qy, "plotpar": [1, 2, 1], "elnum": [5, 6, 7, 8]})

        q8x, q8y = self.quads8()
        self.call("eldraw2", {"ex": q8x, "ey": q8y, "plotpar": [1, 1, 1]})

    def eldisp2(self):
        ex, ey = self.bars()
        self.call("eldisp2", {"ex": ex, "ey": ey, "ed": self.disp(3, 4)})
        self.call("eldisp2", {"ex": ex, "ey": ey, "ed": self.disp(3, 6), "plotpar": [1, 4, 2], "sfac": 250.0})
        self.call("eldisp2", {"ex": ex[1], "ey": ey[1], "ed": self.disp(1, 6)[0], "plotpar": [2, 2, 0]})

        tx, ty = self.triangles()
        self.call("eldisp2", {"ex": tx, "ey": ty, "ed": self.disp(len(tx), 6), "plotpar": [1, 3, 3]})

        qx, qy = self.quads()
        self.call("eldisp2", {"ex": qx, "ey": qy, "ed": self.disp(len(qx), 8), "plotpar": [3, 2, 1],
                              "sfac": 100.0})

        q8x, q8y = self.quads8()
        self.call("eldisp2", {"ex": q8x, "ey": q8y, "ed": self.disp(len(q8x), 16)})

    def dispbeam2(self):
        # calfem-python draws dispbeam2 with a solid line whatever the line type in plotpar
        # (calfem::vis uses it), so the line style is not compared.
        ex, ey = np.array([0.0, 3.0]), np.array([0.0, 4.0])
        self.call("dispbeam2", {"ex": ex, "ey": ey, "edi": self.disp(11, 2)}, compare_style=False)
        self.call("dispbeam2", {"ex": ex, "ey": ey, "edi": self.disp(7, 2), "plotpar": [1, 4, 0], "sfac": 50.0},
                  compare_style=False)
        self.call("dispbeam2", {"ex": np.array([6.0, 0.0]), "ey": np.array([4.0, 4.0]), "edi": self.disp(5, 2),
                                "plotpar": [2, 3, 2]}, compare_style=False)

    def secforce2(self):
        ex, ey = np.array([0.0, 3.0]), np.array([0.0, 4.0])
        es = self.rng.uniform(-1.0, 1.0, 11) * 1e4
        self.call("secforce2", {"ex": ex, "ey": ey, "es": es})
        self.call("secforce2", {"ex": np.array([6.0, 0.0]), "ey": np.array([4.0, 4.0]), "es": es[:5],
                                "plotpar": [4, 1], "sfac": 2e-4})
        eci = np.array([0.0, 0.5, 2.0, 4.5, 5.0])
        self.call("secforce2", {"ex": ex, "ey": ey, "es": es[:5], "plotpar": [3, 2], "sfac": 1e-4, "eci": eci})

    def scalgraph2(self):
        # calfem-python fails for magnitude = [N] (list times float), so only [N, x, y] is used.
        self.call("scalgraph2", {"sfac": 2e-4, "magnitude": [3e4, 1.5, 0.0]})
        self.call("scalgraph2", {"sfac": 25.0, "magnitude": [1e-2, 0.5, -0.25], "plotpar": 1})
        self.call("scalgraph2", {"sfac": 1e4, "magnitude": [2.5e-5, -1.0, 2.0], "plotpar": 4})


if __name__ == "__main__":
    main()
