# CALFEM for C++

A C++ implementation of [CALFEM](https://github.com/CALFEM/calfem-python), the finite element toolbox
from the Division of Structural Mechanics, Lund University. The library is built on
[Eigen](https://eigen.tuxfamily.org) and is used as the solver core of
[ObjectiveFrame](https://github.com/jonaslindemann/objectiveframe).

Currently available:

- **Springs**: `spring1e/s`
- **Bars**: `bar1e/s`, `bar1we/ws` (elastic foundation), `bar2e/s`, `bar2ge/gs` (second order), `bar3e/s`
- **Beams**: `beam1e/s`, `beam1we/ws`, `beam2e/s`, `beam2we/ws`, `beam2te/ts` (Timoshenko), `beam2ge/gs` and
  `beam2gxe/gxs` (second order), `beam2de/ds` (dynamics), `beam3e/s`
- **Field problems** (heat flow etc.): `flw2te/ts`, `flw2qe/qs`, `flw2i4e/s`, `flw2i8e/s`, `flw3i8e/s`
- **Plane stress/strain**: `plante/s/f`, `planqe/s`, `plani4e/s`; **solids**: `soli8e/s`; **plates**: `platre`
- **Materials**: `hooke`, `effmises`, `stress2nodal`
- **Solvers**: `solveq`, `spsolveq`, `SparseSolver`, `red`, `statcon`, `eigen`
- **Time integration**: `step1` (first order systems), `step2` (Newmark), `gfunc` (time functions)
- **Boundary conditions and loads**: `apply_bc`, `apply_bc_3d`, `apply_bc_node`, `apply_force` (`_3d`, `_node`,
  `_total`, `_total_3d`)
- **DOF utilities**: `create_dofs`, `coord_extract`, `extract_edof`, `extract_eldisp`, `convert_to_node_topo`,
  `beam2crd`, `scalfact2`
- Dense and sparse assembly

- **Geometry and mesh generation**: `Geometry` (points, splines, B-splines, circle and ellipse arcs, surfaces
  with holes, structured surfaces and volumes, markers), `create_mesh` / `MeshGenerator` with a built-in 2D
  mesher and an optional Gmsh backend

All element and solver routines of calfem-python's `core.py` are available.

## Requirements

- CMake 3.24 or later
- A C++17 compiler (MSVC 2022+, GCC 9+, Clang 10+)
- [vcpkg](https://github.com/microsoft/vcpkg) (dependencies: Eigen 3.4, Catch2 for tests)

## Building

The presets use Ninja and locate vcpkg through the `VCPKG_ROOT` environment variable. All presets share a
single vcpkg tree, `vcpkg_installed/`, so the dependencies are installed once for all build trees. Since vcpkg
removes packages a configuration does not request, every preset requests the same manifest features (including
`qt`); the Qt presets differ only in enabling `CALFEM_WITH_QT`.

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

On Windows, run these from a Developer Command Prompt/PowerShell (or let VS Code's CMake Tools select the MSVC
kit) so that Ninja finds the compiler.

To get a Visual Studio solution instead, use the `vs` preset (Windows only, Visual Studio 2026). It generates a
multi-config solution in `build-vs/calfem.slnx` that does not need a developer environment:

```bash
cmake --preset vs
cmake --build --preset vs-release   # or vs-debug
ctest --preset vs-release           # or vs-debug
```

Open `build-vs/calfem.slnx` in Visual Studio to build and debug from the IDE. To use another Visual Studio
version, override the generator in a `CMakeUserPresets.json`, e.g. `"generator": "Visual Studio 17 2022"`.

To use a different vcpkg installation than the one in `VCPKG_ROOT`, create a `CMakeUserPresets.json` (ignored by
git):

```json
{
  "version": 6,
  "configurePresets": [
    { "name": "my-debug", "inherits": "debug", "environment": { "VCPKG_ROOT": "E:/vcpkg" } }
  ],
  "buildPresets": [ { "name": "my-debug", "configurePreset": "my-debug" } ],
  "testPresets": [ { "name": "my-debug", "configurePreset": "my-debug" } ]
}
```

### CMake options

| Option                  | Default (top-level / subproject) | Description                          |
|-------------------------|----------------------------------|--------------------------------------|
| `CALFEM_BUILD_TESTS`    | `ON` / `OFF`                     | Build the Catch2 unit tests          |
| `CALFEM_BUILD_EXAMPLES` | `ON` / `OFF`                     | Build the examples in `src/examples` |
| `CALFEM_INSTALL`        | `ON` / `OFF`                     | Generate install/export rules        |
| `CALFEM_WITH_MESH`      | `ON` / `OFF`                     | Geometry and mesh generation (CDT)   |
| `CALFEM_WITH_GMSH`      | `OFF`                            | Gmsh mesh backend (GPL, see below)   |
| `CALFEM_WITH_QT`        | `OFF`                            | Qt windows for `calfem::vis` (LGPL)  |

## Headers

`#include <calfem/calfem.h>` includes the whole library. Individual modules can also be included on their own:

| Header              | Contents                                                                                         |
|---------------------|--------------------------------------------------------------------------------------------------|
| `calfem/types.h`    | `Matrix`, `RowVec`, `ColVec`, `SpMatrix`, `IntMatrix`, ...                                       |
| `calfem/elements.h` | all elements (`elements/spring.h`, `bar.h`, `beam.h`, `flow.h`, `plane.h`, `solid.h`, `plate.h`) |
| `calfem/material.h` | `hooke`, `effmises`                                                                              |
| `calfem/dofs.h`     | `create_dofs`, `coord_extract`, `extract_edof`, `extract_eldisp`, `stress2nodal`, ...            |
| `calfem/assembly.h` | `assem`, `spassem`                                                                               |
| `calfem/boundary.h` | `BoundaryDofs`, `apply_bc`, `apply_force`, ...                                                   |
| `calfem/solvers.h`  | `solveq`, `spsolveq`, `SparseSolver`, `red`, `statcon`, `eigen`, `SolverError`                   |
| `calfem/dynamics.h` | `step1`, `step2`, `gfunc`, `TimeHistory`                                                         |
| `calfem/geometry.h` | `Geometry` (with `CALFEM_WITH_MESH`)                                                             |
| `calfem/mesh.h`     | `create_mesh`, `MeshGenerator`, `Mesh`, `geo_string`, `apply_traction_linear_element`            |
| `calfem/log.h`      | `setLogFunction`                                                                                 |
| `calfem/utils.h`    | `print`, `spy`, `scalfact2`                                                                      |
| `calfem/vis.h`      | `Figure`, `figure`, `eldraw2`, `eldisp2`, `dispbeam2`, `secforce2`, `scalgraph2`, `save_figure`, `show` |
| `calfem/vis_qt.h`   | Qt backend: `figure_widget` (target `calfem::vis_qt`, with `CALFEM_WITH_QT`)                     |

## Examples

Examples live in [src/examples](src/examples) and are built to `bin/debug` or `bin/release`:

- `truss3d` - three-bar 3D truss: topology, assembly, sparse solve and normal forces.
- `frame2d` - 2D portal frame with beam elements, distributed load, dense solve and section forces.
- `plane2d` - plane stress cantilever with 4-node elements: tip deflection versus beam theory, von Mises stress.
- `plate2d` - simply supported plate (`platre`, `apply_bc_3d`) compared with Navier's solution.
- `heat2d` - steady heat conduction in a plate with 4-node isoparametric elements, compared with the exact
  solution.
- `mesh2d` - plate with a hole: geometry, mesh generation, plane stress analysis and the stress concentration
  compared with Kirsch/Heywood. `mesh2d mesh.svg 2` writes the mesh as SVG and uses triangles.
- `vibration2d` - natural frequencies of a cantilever (`beam2de`, `eigen`) and its response to a load pulse
  (`gfunc`, `step2`).

Each example has a visualisation version `<example>_vis` that solves the same problem and plots the results
with `calfem::vis`. The figures are shown in Qt windows when built with `CALFEM_WITH_QT`, otherwise (or with
`--svg`) written as `<example>_<figure>.svg`:

- `truss3d_vis` - undeformed and deformed truss in a 3D projection, supports, load and normal forces.
- `frame2d_vis` - deformed frame and normal force, shear force and moment diagrams.
- `plane2d_vis` - deformed mesh, the centre line deflection versus beam theory, and the effective stress
  (`draw_element_values`).
- `plate2d_vis` - deflected plate as a 3D wireframe, the centre line deflection versus Navier's solution, and
  deflection contours (`draw_nodal_values_contourf`/`_contour`).
- `heat2d_vis` - temperature contours with heat flux arrows (`draw_element_flux`), and the temperature versus
  the exact solution.
- `mesh2d_vis` - mesh, deformed mesh, `sigma_xx` above the hole versus Kirsch's solution, and `sigma_xx`
  contours (`mesh2d_vis 2` uses triangles).
- `vibration2d_vis` - the first three bending modes, and the tip displacement history.

## Visualisation

`calfem/vis.h` ports calfem-python's drawing functions (`calfem.vis_mpl`). They draw into figures, which are
written to SVG without further dependencies, or shown in interactive windows with the Qt backend:

```cpp
#include <calfem/calfem.h>
#include <calfem/vis_qt.h>              // optional: Qt windows (link calfem::vis_qt)

vis::figure(1);
vis::eldraw2(ex, ey, {1, 2, 1});        // undeformed mesh
double sfac = vis::eldisp2(ex, ey, ed); // deformed mesh, returns the scale factor
vis::scalgraph2(sfac, {1e-3, 0.5, 0.0});
vis::title("Displacements");

vis::figure(2);
for (Eigen::Index i = 0; i < ex.rows(); i++)
    vis::secforce2(ex.row(i), ey.row(i), es[i].col(2), {2, 1}, scale);

vis::save_figure("moment.svg");         // .png/.jpg/.pdf need the Qt backend
vis::show();                            // Qt windows; without Qt: SVG opened in the default viewer
```

Mesh and field drawings follow calfem-python's `draw_*` functions and use matplotlib's colour maps (viridis by
default, `vis::colormap()` selects jet or gray) and contour levels:

```cpp
vis::figure(3);
vis::draw_nodal_values_contourf(a, coords, edof, 12, "Temperature", 1, 3, true); // filled contours + mesh
vis::colorbar();
vis::draw_element_flux(q, coords, edof, 1, 3);                                    // flux arrows

vis::figure(4);
vis::draw_element_values(eseff, coords, edof, 2, 3, a, false, true, 10.0);        // on the displaced mesh
vis::colorbar();
```

Also available: `draw_mesh`, `draw_elements`, `draw_node_circles`, `draw_displacements`,
`draw_nodal_values_contour`, `draw_nodal_values_shaded`, `draw_nodal_values` and the element-coordinate forms
`eliso2`, `elval2` and `elflux2`. Deviations from calfem-python (e.g. Gouraud shading approximated by 64
contour intervals) are listed in `vis.h`.

The Qt windows pan with the mouse, zoom with the wheel, fit with Home or a double-click, toggle the grid with G
and save the current view as SVG, PNG or PDF. `show()` blocks until the windows are closed (set
`CALFEM_NO_BLOCK` to return immediately); `figure_widget(fig, parent)` embeds a figure in a Qt application.

The Qt backend is built with `-DCALFEM_WITH_QT=ON` and the vcpkg feature `qt` (presets `vs-qt`, `debug-qt`,
`release-qt`). vcpkg builds Qt from source the first time, which takes a while; the binary cache makes later
configures fast. Qt is LGPL licensed and linked dynamically.

If your user profile path contains a space (e.g. `C:\Users\First Last`), vcpkg's Qt build fails with
`Feature "pkg_config": Forcing to "ON" breaks its condition`. Point `PKG_CONFIG` to vcpkg's `pkg-config.exe`
through a path without spaces (e.g. its 8.3 short name) for the first configure:

```powershell
$env:PKG_CONFIG = "C:/Users/FIRSTL~1/AppData/Local/vcpkg/downloads/tools/msys2/<hash>/mingw64/bin/pkg-config.exe"
$env:VCPKG_KEEP_ENV_VARS = "PKG_CONFIG"
cmake --preset vs-qt
```

## Geometry and mesh generation

Geometries are described as in calfem-python and meshed with `create_mesh` or `MeshGenerator`:

```cpp
Geometry g;
g.point({0.0, 0.0});                 // points 0..3
g.point({2.0, 0.0});
g.point({2.0, 1.0});
g.point({0.0, 1.0});
g.spline({0, 1}, -1, 10);            // curve 0, marker 10
g.spline({1, 2}, -1, 20);
g.spline({2, 3}, -1, 30);
g.spline({3, 0}, -1, 40);
g.surface({0, 1, 2, 3});

auto [coords, edof, dofs, bdofs, element_markers] = create_mesh(g, 3, 0.1, 2); // quads, size, dofs/node

MeshGenerator mesher(g);             // more options, boundary elements
mesher.el_type = 2;
mesher.el_size_factor = 0.05;
Mesh m = mesher.create();
```

Two mesh generators are available (`MeshGenerator::backend`):

- **Built-in** (always available, MIT/MPL): 2D geometries with triangles (`el_type` 2) or quadrilaterals
  (`el_type` 3). Curves are divided according to `el_on_curve` (with progression/bump distributions) or the
  point element sizes, structured surfaces get a mapped mesh, other surfaces a constrained Delaunay
  triangulation ([CDT](https://github.com/artem-ogre/CDT)) with graded interior points and smoothing.
  Quadrilaterals on unstructured surfaces are made by splitting triangles (needs an even number of elements on
  the curves of such surfaces).
- **Gmsh** (optional): 3D, higher order elements, Gmsh's algorithms and options (`min_size`, `max_size`,
  `meshing_algorithm`, `gmsh_options`). The geometry is written as a `.geo` description in the same way as
  calfem-python and meshed with the Gmsh library in-process (no executable), giving the same mesh as
  calfem-python.

`MeshBackend::Auto` (default) uses Gmsh when calfem is built with it, otherwise the built-in mesher.

### Building with Gmsh

**License:** Gmsh is licensed under the GPL (version 2 or later). calfem-cpp itself stays MIT, but programs
that link calfem built with `CALFEM_WITH_GMSH` must comply with the GPL when they are distributed. The option
is therefore off by default.

Gmsh support uses the Gmsh SDK (C API with the `gmsh.h_cwrap` C++ wrapper, which works with any compiler).
The vcpkg `gmsh` port is built without the mesh module and cannot be used. Download the SDK from
[gmsh.info](https://gmsh.info) (or `pip install gmsh`, which contains the same files) and point `GMSH_ROOT` to it:

```bash
cmake --preset release -DCALFEM_WITH_GMSH=ON -DGMSH_ROOT=/path/to/gmsh-sdk
```

The Gmsh shared library is copied next to the test and example executables.

## Testing against calfem-python

The element routines are checked against results computed with calfem-python.
[tests/reference_data/generate_calfem_python.py](tests/reference_data/generate_calfem_python.py) calls the
Python functions with fixed inputs and stores inputs and results in `calfem_python_reference.json`, which
`test_calfem_python.cpp` compares with the C++ results. To regenerate it (requires numpy and scipy):

```bash
python tests/reference_data/generate_calfem_python.py <path to calfem-python>/src
python tests/reference_data/generate_calfem_python_mesh.py <path to calfem-python>/src   # also needs gmsh
```

The mesh reference stores calfem-python's `.geo` text (compared with `geo_string` in every build) and its Gmsh
mesh (compared with the Gmsh backend when available).

## Using the library

As an installed package:

```cmake
find_package(calfem CONFIG REQUIRED)
target_link_libraries(myapp PRIVATE calfem::calfem)
```

As a subproject (git submodule or `FetchContent`):

```cmake
add_subdirectory(external/calfem-cpp)
target_link_libraries(myapp PRIVATE calfem::calfem)
```

```cpp
#include <calfem/calfem.h>

using namespace calfem;

RowVec ex(2), ey(2), ez(2), ep(2);
ex << 0.0, 2.0;
ey << 0.0, 0.0;
ez << 0.0, 0.0;
ep << 210e9, 0.01; // E, A
double eq = 1e3;   // distributed axial load

Matrix Ke = bar3e(ex, ey, ez, ep);          // Python: Ke = bar3e(ex, ey, ez, ep)
auto [Ke2, fe] = bar3e(ex, ey, ez, ep, eq); // Python: Ke, fe = bar3e(ex, ey, ez, ep, eq)
```

The element routines follow calfem-python: results are returned, as a single value or as a `std::tuple` for
use with structured bindings. As in Python, the load vector is only returned when `eq` is given, and the `*s`
routines only return `edi` and `eci` when the number of evaluation points `nep` is given:

```cpp
ColVec es = bar3s(ex, ey, ez, ep, ed);                    // es at the element ends
auto [es, edi, eci] = bar3s(ex, ey, ez, ep, ed, eq, nep); // es, edi, eci at nep points
```

The same style applies to the rest of the API:

```cpp
auto [edof, ex, ey, ez] = extract_edof(dofs, coords, topo);
Matrix ed = extract_eldisp(edof, a);          // Python: ed = extract_eldisp(edof, a)

assem(edof.row(i), K, Ke);                    // Python: assem(edof, K, Ke)
assem(edof.row(i), K, Ke, f, fe);             // Python: assem(edof, K, Ke, f, fe)

auto [a, r] = solveq(K, f, bcDofs);           // Python: a, r = solveq(K, f, bc)
auto [a, r] = spsolveq(K, f, bcDofs, bcVals); // Python: a, r = spsolveq(K, f, bc, bcVal)

SparseSolver solver;
solver.setup(K, bcDofs, bcVals);
auto [a, r] = solver.solve(f);
```

`extract_edof` and `extract_eldisp` are aliases for `extractEdof` and `extractElDisp`, so either naming
style can be used. `assem`/`spassem` update `K` in place (copying the global matrix for every element would be far too slow).
The value-returning solvers throw `calfem::SolverError` (a `std::runtime_error`) on invalid input or a failed
factorization.

The original output-argument form, used by ObjectiveFrame, is still available for all of these. There the
solvers return `false` on failure instead of throwing:

```cpp
Matrix Ke;
ColVec fe;
bar3e(ex, ey, ez, ep, eq, Ke, fe);

ColVec a, r;
bool ok = spsolveq(K, f, bcDofs, bcVals, a, r);
```

Degrees of freedom, topology and boundary condition arrays use 1-based numbering, as in CALFEM for MATLAB and
Python. Boundary condition dofs can be given in any order.

The solvers return `false` on invalid input or a failed factorization. Diagnostic messages (problem size,
the reason for a failure, a warning when K is nearly singular because of missing supports) are passed to an
optional log function; the library prints nothing by default:

```cpp
calfem::setLogFunction([](const std::string &message) { std::cerr << message << "\n"; });
```

## License

MIT, see [LICENSE](LICENSE).
