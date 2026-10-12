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

The presets use Ninja and locate vcpkg through the `VCPKG_ROOT` environment variable.

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

On Windows, run these from a Developer Command Prompt/PowerShell (or let VS Code's CMake Tools select the MSVC
kit) so that Ninja finds the compiler.

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
