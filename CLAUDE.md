# CLAUDE.md

C++ port of CALFEM (Eigen-based), extracted from ObjectiveFrame (`objectiveframe/src/ofsolve/calfem.*`) to become a
standalone library. Reference implementation for new routines: calfem-python
(<https://github.com/CALFEM/calfem-python>, `src/calfem/core.py`, `geometry.py`, `mesh.py`).

## Build and test

```bash
cmake --preset debug && cmake --build --preset debug && ctest --preset debug
```

Requires `VCPKG_ROOT` and, on Windows, an MSVC developer environment (Ninja generator). The vcpkg
`builtin-baseline` matches ObjectiveFrame's so both resolve the same Eigen version; keep them (and the commit in
`.github/workflows/ci.yml`) in sync. Build trees are `build-debug`/`build-release`; executables go to
`bin/debug`/`bin/release`. On Windows the `vs` preset generates a multi-config Visual Studio 2026 solution in
`build-vs` (build/test presets `vs-debug`/`vs-release`), sharing the same `bin/<config>` output directories. vcpkg packages are shared by all build trees in a
single `vcpkg_installed/`; all presets request the same `VCPKG_MANIFEST_FEATURES` (vcpkg removes unrequested
packages), and the Qt presets only add `CALFEM_WITH_QT=ON`. Keep the feature lists identical.

## Layout

- `include/calfem/` public headers, `src/calfem/` sources (same structure), target `calfem::calfem` (static by
  default). `calfem.h` is an umbrella header including everything; modules are `types.h`, `log.h`, `utils.h`,
  `material.h`, `dofs.h`, `assembly.h`, `boundary.h`, `solvers.h`, `dynamics.h`, `vis.h` and `elements/*.h` (one header per element family, collected by
  `elements.h`). New element families get `include/calfem/elements/<family>.h` + `src/calfem/elements/<family>.cpp`
  and an include in `elements.h`; future modules (materials, geometry, mesh) follow the same pattern.
- `src/calfem/detail/` private headers (not installed), e.g. `detail/log.h` with `logMessage`/`logError`.
- `src/calfem/mesh/` geometry meshing (only with `CALFEM_WITH_MESH`): `curves.*` curve evaluation/division,
  `builtin_mesher.cpp` (2D: transfinite + CDT), `geo_writer.cpp` (port of calfem-python's `_writeGeoFile`; keep it
  identical, it is compared with Python's output), `gmsh_mesher.cpp` (Gmsh API via `gmsh.h_cwrap`, only with
  `CALFEM_WITH_GMSH`), `mesh.cpp` (front end, dofs/bdofs). Gmsh is GPL: never make it a default dependency.
- `src/calfem/vis/` visualisation (port of calfem-python's `vis_mpl`): `figure.cpp` (Figure, current-figure
  management, layout incl. colour bar, backend dispatch), `svg.cpp` (SVG writer), `draw2.cpp` (eldraw2,
  eldisp2, ...), `draw_fields.cpp` (draw_mesh, draw_element_values, draw_nodal_values_contourf/contour/shaded,
  eliso2, elval2, elflux2: contouring on triangles with matplotlib's MaxNLocator levels and colour
  normalisation), `colormap.cpp` (matplotlib's viridis/jet/gray tables; colours must match matplotlib exactly).
  `src/calfem/detail/vis_codes.h` holds the plotpar colour/line codes. Drawing functions only add items (Line,
  Markers, Text, Polygon) to a `Figure` and set its colour scale; renderers use `vis::layout()` and must handle
  every item type and the colour bar. `src/calfem_qt/` is the optional Qt
  backend (target `calfem::vis_qt`, `CALFEM_WITH_QT`, vcpkg feature `qt`): `painter.cpp` (QPainter renderer,
  PNG/PDF), `figure_widget.cpp` (widget/window), `qt_backend.cpp` (`show`, installed by including
  `vis_qt.h`). `vis_qt.h` only forward-declares `QWidget`, so it passes the header check without Qt. The core
  library never depends on Qt. Reference drawings: `generate_calfem_python_vis.py` (needs matplotlib) and
  `test_vis_reference.cpp` (item by item); `generate_calfem_python_vis_fields.py` and
  `test_vis_fields_reference.cpp` for the field functions (filled contours compared by level, colour and area
  per level interval, contour lines by colour and length, since matplotlib merges regions); Qt tests are in `tests/qt/` (`minimal` platform; vcpkg has no `offscreen` plugin on Windows; `calfem_deploy_qt_plugins()` copies the plugins next to executables).
- Every public header must compile on its own (checked by the `calfem_header_check` target in `tests/`).
- `src/calfem/detail/element_utils.h` shared element helpers (evaluation points, 2D transformations,
  axial/bending section results); `detail/gauss.h` Gauss rules (CALFEM constants and point order);
  `detail/shape_functions.h` isoparametric shape functions (quad4, quad8, brick8).
- `tests/reference_data/` calfem-python reference results. When porting a routine: add cases to
  `generate_calfem_python.py`, regenerate the JSON (needs numpy/scipy and a calfem-python checkout), and add a
  handler in `test_calfem_python.cpp` (it fails if a handler has no cases). Meshing has its own generator
  (`generate_calfem_python_mesh.py`, needs the `gmsh` Python package) and test (`test_mesh_reference.cpp`).
- Gmsh build for testing: `cmake --preset release -B build-gmsh -DCALFEM_WITH_GMSH=ON -DGMSH_ROOT=<sdk>
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=<dir>` (the `gmsh` pip package is a usable SDK).
- `tests/` Catch2 tests, one `test_*.cpp` per area, picked up by glob.
- `tests/reference/` frozen copy of the calfem code ObjectiveFrame was built against (namespace `calfem_ref`).
  `test_objectiveframe_compat.cpp` checks signatures and compares results against it. Never edit the reference.
- `src/examples/` example programs (listed in `src/examples/CMakeLists.txt`). Every example `<name>` has a
  visualisation version `<name>_vis.cpp` (added automatically by the CMake list) using the shared helpers in
  `example_vis.h` (graphs, arrows, 3D projection, `finish()`: Qt windows or `--svg`).
- `cmake/calfemConfig.cmake.in` package config for `find_package(calfem)`.

## Constraints

- **ObjectiveFrame compatibility**: do not change the signatures or semantics of `bar3e`, `bar3s`, `beam3e`,
  `beam3s`, `spassem`, `extractAllDofs`, `toZeroIndex`, `SparseSolver` (`setup`/`solve`/`recompute`) or the
  type aliases in `calfem.h`. Add new APIs as overloads instead.
- **API style**: functions follow calfem-python and return their results (a single value, or a `std::tuple`
  when Python returns a tuple), using overloads by argument count for Python's optional arguments (`eq`, `nep`,
  `bcVal`, `f`/`fe`). Value-returning solvers throw `SolverError`; report failures with `detail::logError()`
  (`src/calfem/detail/log.h`) so the exception carries the message. The output-argument form exists only for the original
  functions; new functions get the returning form only. Overloads must not share an argument count with an
  existing output-argument form. Exception: `assem`/`spassem` update `K` in place.
- **Naming**: use calfem-python names for new functions. Older camelCase names get inline snake_case aliases
  next to the originals (`extract_eldisp` for `extractElDisp`, `extract_edof` for `extractEdof`) covering all overloads.
- **Porting from calfem-python**: keep the numerics of `core.py`; deviations must be documented in the header
  (e.g. evaluation points computed as `i * L/(nep-1)` instead of `np.arange`, which can produce an extra point;
  `beam2de` always returning `Ce`; Python debug prints dropped; `step1`/`step2` recording `dofhist` also
  without `times`; `eigen` requiring symmetric K and SPD M; correct 27-point rule in `flw3i8e/s`; `flw2qs`
  accepting several rows in `ed`; `*e` routines called with `eq` always returning `(Ke, fe)` and without `eq` only `Ke`, also for `plani4e`;
  `plants` evaluating every row of `ed` for a single element). Return order follows calfem-python even where
  it differs from MATLAB (`soli8s` returns `et, es, eci`). `apply_force_node` loads the selected dof
  (calfem-python is off by one); `coord_extract` always returns `ex, ey, ez`; `convert_to_node_topo` returns
  zero-based `topo` like calfem-python. Not ported: `dof_hash`/`c_mul` (Python hashing helpers);
  `apply_traction_linear_element` belongs with the mesh module (needs boundary elements). calfem-python has
  bugs; verify suspicious behaviour numerically before copying it, and leave such cases out of the reference
  data. Single-value loads are `double`, multi-component
  loads `RowVec`.
- DOF, topology and BC arrays are 1-based, as in CALFEM.
- The library never writes to stdout/stderr (except `print`/`spy`); diagnostics go through `setLogFunction`.
- Library targets C++17; do not require newer standard features in public headers.
- Code style follows `.clang-format` (Microsoft-based, 4-space indent).
