"""Generate reference results from calfem-python's vis_mpl mesh and field drawing functions.

Every case calls a calfem-python function (draw_mesh, draw_element_values,
draw_nodal_values_contourf, ..., eliso2, elval2, elflux2) with fixed inputs on an empty
matplotlib figure (Agg backend) and records what test_vis_fields_reference.cpp compares:

- lines (mesh edges) and node markers, as generate_calfem_python_vis.py,
- filled polygons with their face and edge colours (element values),
- filled contours: the levels and, per level interval, the colour and the filled area
  (matplotlib merges the regions, calfem::vis draws one polygon per triangle),
- contour lines: the levels and, per level, the colour and the total line length,
- colour limits of shaded drawings (tripcolor),
- arrows (quiver): centres, vectors and colours,
- the returned scale factor or levels.

Usage (needs numpy and matplotlib):

    python generate_calfem_python_vis_fields.py <path to calfem-python/src> [output.json]
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
                                                             "calfem_python_vis_fields_reference.json")

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
        "generator": "tests/reference_data/generate_calfem_python_vis_fields.py",
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


def polygon_area(v):
    x, y = v[:, 0], v[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def path_area(path):
    """Filled area of a contour path (holes have the opposite orientation)."""
    return abs(sum(polygon_area(p) for p in path.to_polygons(closed_only=False) if len(p) >= 3))


def path_length(path):
    """Total length of the line segments of a path."""
    total = 0.0
    for p in path.to_polygons(closed_only=False):
        if len(p) >= 2:
            total += float(np.sum(np.hypot(np.diff(p[:, 0]), np.diff(p[:, 1]))))
    return total


def harvest(fig):
    """Lines, markers, polygons, contours, colour limits and arrows drawn in fig."""
    import matplotlib.collections as mcoll
    import matplotlib.contour as mcontour
    import matplotlib.quiver as mquiver

    fig.canvas.draw()  # updates the colour mapping of the collections

    result = {"lines": [], "markers": [], "polygons": [], "contourf": [], "contour": [], "clim": None,
              "arrows": []}
    for ax in fig.axes:
        for line in ax.lines:
            xy = line.get_xydata()
            if line.get_linestyle() in ("None", "", " ") and line.get_marker() not in ("None", None, ""):
                result["markers"].append({"points": xy.tolist(), "marker": "circle",
                                          "color": rgb(line.get_markerfacecolor())})
            else:
                result["lines"].append({"points": xy.tolist(), "closed": False, "color": rgb(line.get_color()),
                                        "width": float(line.get_linewidth())})

        for coll in ax.collections:
            if isinstance(coll, mcontour.ContourSet):
                levels = [float(v) for v in coll.levels]
                if coll.filled:
                    colors = coll.get_facecolor()
                    bands = [{"color": rgb(colors[i]), "area": path_area(p)}
                             for i, p in enumerate(coll.get_paths())]
                    result["contourf"].append({"levels": levels, "bands": bands})
                else:
                    colors = coll.get_edgecolor()
                    lines = [{"color": rgb(colors[i if len(colors) > 1 else 0]), "length": path_length(p)}
                             for i, p in enumerate(coll.get_paths())]
                    result["contour"].append({"levels": levels, "lines": lines})
            elif isinstance(coll, mquiver.Quiver):
                offsets = np.asarray(coll.get_offsets())
                colors = coll.get_facecolor()
                for i in range(len(offsets)):
                    result["arrows"].append({
                        "center": offsets[i].tolist(),
                        "vector": [float(coll.U[i]), float(coll.V[i])],
                        "color": rgb(colors[i if len(colors) > 1 else 0]),
                    })
            elif isinstance(coll, mcoll.TriMesh):  # tripcolor with Gouraud shading
                result["clim"] = [float(v) for v in coll.get_clim()]
            elif isinstance(coll, mcoll.PathCollection):  # scatter
                face = coll.get_facecolor()
                edge = coll.get_edgecolor()
                filled = len(face) > 0 and face[0][3] > 0
                result["markers"].append({
                    "points": np.asarray(coll.get_offsets()).tolist(),
                    "marker": "circle",
                    "color": rgb(face[0] if filled else edge[0]),
                    "filled": bool(filled),
                })
            elif isinstance(coll, mcoll.PolyCollection):
                face = coll.get_facecolor()
                edge = coll.get_edgecolor()
                width = float(coll.get_linewidth()[0])
                mapped = coll.get_array() is not None
                if mapped:
                    result["clim"] = [float(v) for v in coll.get_clim()]
                for i, path in enumerate(coll.get_paths()):
                    verts = path.vertices
                    closed = path.codes is not None and path.codes[-1] == path.CLOSEPOLY
                    if closed:
                        verts = verts[:-1]
                    has_face = len(face) > 0 and face[i % len(face)][3] > 0
                    has_edge = len(edge) > 0 and edge[i % len(edge)][3] > 0
                    if not has_face:
                        result["lines"].append({"points": verts.tolist(), "closed": bool(closed),
                                                "color": rgb(edge[i % len(edge)]), "width": width})
                    else:
                        result["polygons"].append({
                            "points": verts.tolist(),
                            "fill": rgb(face[i % len(face)]),
                            "edge": rgb(edge[i % len(edge)]) if has_edge else None,
                            "width": width,
                        })
    return result


class Generator:
    def __init__(self, cfv):
        self.cfv = cfv
        self.cases = []
        self.rng = np.random.default_rng(20261012)

    def call(self, name, inputs, **kwargs):
        """Draws with cfv.<name>(**inputs) on a new figure and stores the case."""
        import matplotlib.pyplot as plt

        func = getattr(self.cfv, name)
        fig = plt.figure()
        result = func(**inputs)
        drawing = harvest(fig)
        plt.close(fig)

        case = {"function": name, "inputs": {k: to_json(v) for k, v in inputs.items()}}
        if isinstance(result, (float, int, np.floating)):
            case["result"] = float(result)
        case.update(drawing)
        case.update(kwargs)
        self.cases.append(case)

    # ---- Meshes ----

    def grid(self, nx=4, ny=3):
        """Distorted structured grid: node coordinates and quadrilateral topology (1-based)."""
        coords = []
        for j in range(ny + 1):
            for i in range(nx + 1):
                x = 2.0 * i / nx + 0.06 * np.sin(1.7 * j + i)
                y = 1.5 * j / ny + 0.05 * np.cos(2.3 * i - j)
                coords.append([x, y])
        topo = []
        for j in range(ny):
            for i in range(nx):
                n = j * (nx + 1) + i + 1
                topo.append([n, n + 1, n + nx + 2, n + nx + 1])
        return np.array(coords), np.array(topo)

    def triangles(self, topo):
        return np.vstack([topo[:, [0, 1, 2]], topo[:, [0, 2, 3]]])

    def edof(self, topo, dofs_per_node):
        cols = []
        for n in topo.T:
            for d in range(dofs_per_node):
                cols.append((n - 1) * dofs_per_node + d + 1)
        return np.array(cols).T

    def field(self, coords):
        x, y = coords[:, 0], coords[:, 1]
        return 3.7 * x ** 2 - 2.1 * y + 0.4 * x * y + 1.3

    def element_coords(self, coords, topo):
        return coords[topo - 1, 0], coords[topo - 1, 1]

    def quads8(self, coords, topo):
        ex, ey = self.element_coords(coords, topo)
        mx = 0.5 * (ex + np.roll(ex, -1, axis=1))
        my = 0.5 * (ey + np.roll(ey, -1, axis=1))
        return np.hstack([ex, mx]), np.hstack([ey, my])

    # ---- Cases ----

    def generate_all(self):
        self.mesh_cases()
        self.element_value_cases()
        self.flux_cases()
        self.nodal_value_cases()
        self.eliso2_cases()

    def mesh_cases(self):
        coords, topo = self.grid()
        tri = self.triangles(topo)

        self.call("draw_mesh", {"coords": coords, "edof": self.edof(topo, 2), "dofs_per_node": 2, "el_type": 3})
        self.call("draw_mesh", {"coords": coords, "edof": self.edof(tri, 1), "dofs_per_node": 1, "el_type": 2,
                                "color": (0.0, 0.0, 1.0), "filled": True, "show_nodes": True,
                                "node_color": (1.0, 0.0, 0.0)})

        ex, ey = self.element_coords(coords, topo)
        self.call("draw_elements", {"ex": ex, "ey": ey, "color": (1.0, 0.0, 0.0)})
        self.call("draw_elements", {"ex": ex[:3], "ey": ey[:3], "filled": True, "face_color": (0.2, 0.6, 0.4)})
        self.call("draw_node_circles", {"ex": ex[:2], "ey": ey[:2]})
        self.call("draw_node_circles", {"ex": ex[:2], "ey": ey[:2], "filled": True, "color": (0.0, 0.0, 1.0)})

        a = self.rng.uniform(-1.0, 1.0, 2 * len(coords)) * 1e-2
        self.call("draw_displacements", {"a": a, "coords": coords, "edof": self.edof(topo, 2), "dofs_per_node": 2,
                                         "el_type": 3, "magnfac": 5.0, "color": (0.3, 0.3, 0.3)})
        # Automatic magnification factor (calfem-python ignores color and always draws (0.3, 0.3, 0.3))
        self.call("draw_displacements", {"a": a.reshape(-1, 2), "coords": coords, "edof": self.edof(tri, 2),
                                         "dofs_per_node": 2, "el_type": 2, "draw_undisplaced_mesh": True},
                  magnfac=float(0.1 * np.max(np.ptp(coords, axis=0)) /
                                np.max(np.linalg.norm(a.reshape(-1, 2), axis=1))))

    def element_value_cases(self):
        coords, topo = self.grid()
        tri = self.triangles(topo)
        ev = self.rng.uniform(-3.0, 5.0, len(topo))

        self.call("draw_element_values", {"values": ev, "coords": coords, "edof": self.edof(topo, 2),
                                          "dofs_per_node": 2, "el_type": 3})
        a = self.rng.uniform(-1.0, 1.0, 2 * len(coords)) * 1e-2
        self.call("draw_element_values", {"values": self.rng.uniform(0.0, 1.0, len(tri)), "coords": coords,
                                          "edof": self.edof(tri, 2), "dofs_per_node": 2, "el_type": 2,
                                          "displacements": a, "magnfac": 3.0, "draw_elements": False,
                                          "draw_undisplaced_mesh": True})

        ex, ey = self.element_coords(coords, topo)
        self.call("elval2", {"ex": ex, "ey": ey, "ev": ev})
        self.call("elval2", {"ex": ex, "ey": ey, "ev": ev, "showMesh": True})

    def flux_cases(self):
        coords, topo = self.grid()
        tri = self.triangles(topo)
        flux = self.rng.uniform(-2.0, 2.0, (len(topo), 2))

        self.call("draw_element_flux", {"flux": flux, "coords": coords, "edof": self.edof(topo, 1),
                                        "dofs_per_node": 1, "el_type": 3})
        self.call("draw_element_flux", {"flux": self.rng.uniform(-2.0, 2.0, (len(tri), 2)), "coords": coords,
                                        "edof": self.edof(tri, 1), "dofs_per_node": 1, "el_type": 2, "scale": 0.05,
                                        "color_by_magnitude": True, "draw_elements": True})

        ex, ey = self.element_coords(coords, topo)
        self.call("elflux2", {"ex": ex, "ey": ey, "es": flux})
        tx, ty = self.element_coords(coords, tri)
        self.call("elflux2", {"ex": tx, "ey": ty, "es": self.rng.uniform(-2.0, 2.0, (len(tri), 2)),
                              "plotcolor": [4], "sfac": 0.04})

    def nodal_value_cases(self):
        coords, topo = self.grid()
        tri = self.triangles(topo)
        v = self.field(coords)

        self.call("draw_nodal_values_contourf", {"values": v, "coords": coords, "edof": topo})
        self.call("draw_nodal_values_contourf", {"values": v, "coords": coords, "edof": tri, "levels": 5})
        self.call("draw_nodal_values_contourf", {"values": v, "coords": coords, "edof": topo,
                                                 "levels": [0.0, 2.5, 5.0, 10.0, 14.0]})
        self.call("draw_nodal_values_contourf", {"values": v * 1e-4 + 2.0, "coords": coords, "edof": topo,
                                                 "levels": 7, "dofs_per_node": 1, "el_type": 3,
                                                 "draw_elements": True})

        self.call("draw_nodal_values_contour", {"values": v, "coords": coords, "edof": topo})
        self.call("draw_nodal_values_contour", {"values": -v, "coords": coords, "edof": tri, "levels": 4})
        self.call("draw_nodal_values_contour", {"values": v, "coords": coords, "edof": topo,
                                                "levels": [1.0, 3.3, 8.0]})

        self.call("draw_nodal_values_shaded", {"values": v, "coords": coords, "edof": topo})
        self.call("draw_nodal_values_shaded", {"values": v, "coords": coords, "edof": tri, "dofs_per_node": 1,
                                               "el_type": 2, "draw_elements": True})

        self.call("draw_nodal_values", {"node_vals": v, "coords": coords, "edof": self.edof(topo, 2),
                                        "dofs_per_node": 2, "el_type": 3})
        self.call("draw_nodal_values", {"node_vals": v, "coords": coords, "edof": self.edof(tri, 1),
                                        "dofs_per_node": 1, "el_type": 2, "clim": (2.0, 9.0),
                                        "draw_elements": False})

    def eliso2_cases(self):
        coords, topo = self.grid()
        tri = self.triangles(topo)

        ex, ey = self.element_coords(coords, topo)
        ed = self.field(np.column_stack((ex.ravel(), ey.ravel()))).reshape(ex.shape)
        self.call("eliso2", {"ex": ex, "ey": ey, "ed": ed})

        tx, ty = self.element_coords(coords, tri)
        td = self.field(np.column_stack((tx.ravel(), ty.ravel()))).reshape(tx.shape)
        self.call("eliso2", {"ex": tx, "ey": ty, "ed": td, "isov": 6, "plotpar": [2, 4]})

        q8x, q8y = self.quads8(coords, topo)
        q8d = self.field(np.column_stack((q8x.ravel(), q8y.ravel()))).reshape(q8x.shape)
        self.call("eliso2", {"ex": q8x, "ey": q8y, "ed": q8d, "isov": [2.0, 4.0, 8.0, 12.0]})


if __name__ == "__main__":
    main()
