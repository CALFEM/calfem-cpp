"""Generate geometry/mesh reference results from calfem-python for the calfem-cpp tests.

For each case a geometry is described as a list of operations (point, spline, circle, ...),
built with calfem.geometry, written as a Gmsh .geo text by calfem-python's mesh generator and
meshed with calfem.mesh.create_mesh (Gmsh Python module). The results are stored in
calfem_python_mesh_reference.json; test_mesh_reference.cpp replays the operations with the C++
Geometry and compares the .geo text and (when calfem is built with Gmsh) the mesh.

Usage (needs numpy, scipy and gmsh):

    python generate_calfem_python_mesh.py <path to calfem-python/src> [output.json]
"""

import io
import json
import math
import os
import sys

import numpy as np


def square_tri():
    ops = [{"op": "point", "coord": c} for c in ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0])]
    ops += [{"op": "spline", "points": [i, (i + 1) % 4], "marker": 10 * (i + 1)} for i in range(4)]
    ops += [{"op": "surface", "outer_loop": [0, 1, 2, 3], "marker": 1}]
    return {"name": "square_tri", "ops": ops, "el_type": 2, "el_size_factor": 0.3, "dofs_per_node": 1}


def plate_with_hole_quads():
    w, h, r = 2.0, 1.0, 0.25
    ops = [{"op": "point", "coord": c, "el_size": 0.2} for c in ([0.0, 0.0], [w, 0.0], [w, h], [0.0, h])]
    ops += [{"op": "point", "coord": c, "el_size": 0.08} for c in
            ([1.0, 0.5], [1.0 + r, 0.5], [1.0, 0.5 + r], [1.0 - r, 0.5], [1.0, 0.5 - r])]
    ops += [{"op": "spline", "points": [i, (i + 1) % 4], "marker": [10, 20, 30, 40][i]} for i in range(4)]
    ops += [{"op": "circle", "points": [5 + i, 4, 5 + (i + 1) % 4], "marker": 50} for i in range(4)]
    ops += [{"op": "surface", "outer_loop": [0, 1, 2, 3], "holes": [[4, 5, 6, 7]], "marker": 100}]
    return {"name": "plate_with_hole_quads", "ops": ops, "el_type": 3, "el_size_factor": 1.0, "dofs_per_node": 2}


def structured():
    ops = [{"op": "point", "coord": c} for c in ([0.0, 0.0], [3.0, 0.0], [3.0, 1.0], [0.0, 1.0])]
    ops += [{"op": "spline", "points": [0, 1], "marker": 10, "el_on_curve": 8, "el_distrib_type": "progression",
             "el_distrib_val": 1.2},
            {"op": "spline", "points": [1, 2], "marker": 20, "el_on_curve": 4},
            {"op": "spline", "points": [2, 3], "marker": 30, "el_on_curve": 8, "el_distrib_type": "progression",
             "el_distrib_val": 1 / 1.2},
            {"op": "spline", "points": [3, 0], "marker": 40, "el_on_curve": 4}]
    ops += [{"op": "struct_surface", "outer_loop": [0, 1, 2, 3], "marker": 7}]
    return {"name": "structured", "ops": ops, "el_type": 3, "el_size_factor": 1.0, "dofs_per_node": 2}


def curves():
    ops = [{"op": "point", "coord": c, "el_size": 0.15} for c in
           ([0.0, 0.0], [2.0, 0.0], [0.0, 1.0], [-1.0, 1.5], [-1.5, 0.5])]
    ops[1]["marker"] = 5  # marked point
    ops += [{"op": "ellipse", "points": [1, 0, 1, 2], "marker": 10},
            {"op": "bspline", "points": [2, 3, 4, 0], "marker": 0},
            {"op": "spline", "points": [0, 1], "marker": 0}]
    ops += [{"op": "surface", "outer_loop": [0, 1, 2], "marker": 3}]
    return {"name": "curves", "ops": ops, "el_type": 2, "el_size_factor": 1.0, "dofs_per_node": 1}


def two_surfaces():
    ops = [{"op": "point", "coord": c, "el_size": 0.25} for c in
           ([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0])]
    ops += [{"op": "spline", "points": p, "marker": m} for p, m in
            ([[0, 1], 1], [[1, 2], 2], [[2, 3], 1], [[3, 0], 1], [[1, 4], 1], [[4, 5], 3], [[5, 2], 1])]
    ops += [{"op": "surface", "outer_loop": [0, 1, 2, 3], "marker": 11},
            {"op": "surface", "outer_loop": [4, 5, 6, 1], "marker": 12}]
    return {"name": "two_surfaces", "ops": ops, "el_type": 2, "el_size_factor": 1.0, "dofs_per_node": 1}


def box_hex():
    a, b, c, n = 2.0, 1.0, 1.0, [4, 2, 2]
    pts = [[0, 0, 0], [a, 0, 0], [a, b, 0], [0, b, 0], [0, 0, c], [a, 0, c], [a, b, c], [0, b, c]]
    ops = [{"op": "point", "coord": [float(x) for x in p]} for p in pts]
    edges = [([0, 1], 0), ([1, 2], 1), ([2, 3], 0), ([3, 0], 1),
             ([4, 5], 0), ([5, 6], 1), ([6, 7], 0), ([7, 4], 1),
             ([0, 4], 2), ([1, 5], 2), ([2, 6], 2), ([3, 7], 2)]
    ops += [{"op": "spline", "points": p, "el_on_curve": n[d]} for p, d in edges]
    faces = [([0, 1, 2, 3], 10), ([4, 5, 6, 7], 20), ([0, 9, 4, 8], 30), ([2, 11, 6, 10], 40),
             ([1, 10, 5, 9], 50), ([3, 8, 7, 11], 60)]
    ops += [{"op": "struct_surface", "outer_loop": l, "marker": m} for l, m in faces]
    ops += [{"op": "struct_volume", "outer_surfaces": [0, 1, 2, 3, 4, 5], "marker": 90}]
    return {"name": "box_hex", "ops": ops, "el_type": 5, "el_size_factor": 1.0, "dofs_per_node": 3}


CASES = [square_tri, plate_with_hole_quads, structured, curves, two_surfaces, box_hex]


def build(cfg, ops):
    g = cfg.Geometry()
    for op in ops:
        kind = op["op"]
        if kind == "point":
            g.point(list(op["coord"]), marker=op.get("marker", 0), el_size=op.get("el_size", 1))
        elif kind in ("spline", "bspline", "circle", "ellipse"):
            getattr(g, kind)(op["points"], marker=op.get("marker", 0), el_on_curve=op.get("el_on_curve"),
                             el_distrib_type=op.get("el_distrib_type"), el_distrib_val=op.get("el_distrib_val"))
        elif kind == "surface":
            g.surface(op["outer_loop"], holes=op.get("holes", []), marker=op.get("marker", 0))
        elif kind == "struct_surface":
            g.struct_surface(op["outer_loop"], marker=op.get("marker", 0))
        elif kind == "struct_volume":
            g.struct_volume(op["outer_surfaces"], marker=op.get("marker", 0))
        else:
            raise ValueError(kind)
    return g


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    src = os.path.abspath(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__),
                                                             "calfem_python_mesh_reference.json")
    sys.path.insert(0, src)
    import calfem
    import calfem.geometry as cfg
    import calfem.mesh as cfm
    import gmsh

    cases = []
    for make in CASES:
        case = make()
        g = build(cfg, case["ops"])

        gen = cfm.GmshMeshGenerator(g, el_type=case["el_type"], el_size_factor=case["el_size_factor"],
                                    dofs_per_node=case["dofs_per_node"])
        gen.geofile = io.StringIO()
        gen._writeGeoFile()
        case["geo"] = gen.geofile.getvalue()

        coords, edof, dofs, bdofs, markers = cfm.create_mesh(g, el_type=case["el_type"],
                                                             el_size_factor=case["el_size_factor"],
                                                             dofs_per_node=case["dofs_per_node"])
        case["coords"] = np.asarray(coords, dtype=float).tolist()
        case["edof"] = np.asarray(edof, dtype=int).tolist()
        case["dofs"] = np.asarray(dofs, dtype=int).tolist()
        case["bdofs"] = {str(k): sorted(int(d) for d in v) for k, v in bdofs.items()}
        case["element_markers"] = [int(m) for m in markers]
        cases.append(case)
        print(f"  {case['name']}: {len(coords)} nodes, {len(edof)} elements")

    data = {"generator": "tests/reference_data/generate_calfem_python_mesh.py",
            "calfem_python_version": calfem.__version__, "gmsh_version": gmsh.__version__, "cases": cases}

    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)
        f.write("\n")
    print(f"Wrote {len(cases)} cases to {out}")


if __name__ == "__main__":
    main()
