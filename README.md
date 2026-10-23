# GeoPartitioning

Spatial and spatiotemporal domain decomposition for geostatistical workflows.
Extracted verbatim from the `bstm` package's `src/partitioning.jl` engine
(v1.0.0) so that it can be shared by `bstm` and `GeoMotion` — and used
standalone — without any dependency between the three packages in the other
direction.

## What it provides

- **Areal-unit assignment** (`assign_spatial_units`,
  `assign_spatial_units_inferred`): Voronoi-centered (`:cvt`), K-means-seeded
  (`:kvt`), quasi-optimal (`:qvt`), balanced (`:bvt`), adaptive-capacity
  (`:avt`), hierarchical (`:hvt`), regular lattice, and hexagonal
  tessellations, with exact-unit control, area/merge constraints, and
  polygon/cityblock metrics.
- **Graph topology** (`spatial_weights_matrix`, `spatial_knn_graph`,
  `spatial_radius_graph`, `scaling_factor_bym2`, and the internal lattice
  adjacency builder): weighting-matrix styles (binary/row-standardized),
  connectivity enforcement, and BYM2 precision scaling.
- **Temporal & spatiotemporal discretization** (`discretize_data`,
  `assign_time_units`, `assign_spatiotemporal_units`): quantile-regular,
  Jenks-natural, k-means, and pan-biological time binning.
- **Cross-validation** (`spatial_block_cv`): k-means and grid spatial block
  folds.
- **Hexagonal meshing** (`build_hex_mesh_planar`, `map_point_to_units`,
  `lonlat_to_xy_km`, `xy_km_to_lonlat`): tangent-plane / CRS-aware planar hex
  meshes with nearest-centroid point mapping and KDTree lookup.
- **HSI field utilities** (`_interpolate_hsi`,
  `build_monthly_hsi_matrix`, `load_hsi_jld2`, `reshard_hsi_field`): habitat
  suitability index interpolation and field resharding carried over from the
  original engine.
- **Column vocabulary** (`STANDARD_SPATIAL_COORDINATE_PAIRS`,
  `STANDARD_TEMPORAL_CANDIDATES`, plus the internal `_detect_xy_columns` /
  `_detect_time_column` resolvers).

Full methodology, mathematics, and worked examples: see
[`docs/spatial_partitioning.md`](docs/spatial_partitioning.md) — ported from
`bstm`'s partitioning guide.

## Consumers

- **`bstm`** (`../bstm`): imports the public API and re-exports it unchanged;
  `bstm`'s `:partitioning` segment keeps testing the original behaviour
  through this package.
- **`GeoMotion`** (`../GeoMotion`): imports the geometric/hex-mesh primitives;
  keeps its domain-specific `assign_spatial_units` (hexagonal-only),
  `map_to_units`, and `build_hex_mesh_clipped` extension, and adds a
  land/water-masked `map_point_to_units` method in `movement.jl`.

## Setup

```sh
julia --project=. scripts/setup_env.jl   # instantiate + resolve
julia --project=. test/runtests.jl       # run the test suite
julia --project=. scripts/parse_check.jl # parse-only sanity check
```

## Provenance & deliberate fixes

Per-file provenance headers map every `src/*.jl` file back to the exact line
range of the original `bstm/src/partitioning.jl` (a 3,039-line file split at
its section banners; reassembly is byte-identical — see the fidelity check in
the repo history). Two latent defects found in the original were fixed here
rather than propagated:

1. `_HAS_COORDSYS` was referenced (`hexmesh.jl`, from `lonlat_to_xy_km` /
   `xy_km_to_lonlat`) but never defined anywhere in `bstm`, so any call with
   `crs ≠ nothing` raised `UndefVarError`. `CoordRefSystems` is a hard
   dependency here, so the flag is defined as `true` and declared before the
   includes, preserving the historical branch shape.
2. `load_shapefile_to_libgeos` calls `Shapefile.Table`, but `Shapefile` was
   never declared in `bstm`'s project (guaranteed `UndefVarError`). It is a
   hard dependency here.
