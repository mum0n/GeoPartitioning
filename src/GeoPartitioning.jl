"""
    GeoPartitioning

Spatial and spatiotemporal domain decomposition for geostatistical workflows:
Voronoi/CVT/KVT/QVT/BVT/HVT/AVT tessellation, hexagonal gridding, areal-unit
assignment, spatial graph construction and weighting matrices, temporal and
spatiotemporal discretization, and spatial block cross-validation.

Extracted from the `bstm` package's `src/partitioning.jl` (v1.0.0); see the
per-file provenance headers in `src/` and `docs/spatial_partitioning.md` for
the full methodology documentation.
"""
module GeoPartitioning

# External dependencies used by the extracted files. Everything the extracted
# code needs is listed here explicitly: no `bstm` dependency is permitted in
# either direction from this package (bstm consumes *this* package).
using
    Clustering,
    CoordRefSystems,
    DataFrames,
    Dates,
    DelaunayTriangulation,
    Graphs,
    JLD2,
    LibGEOS,
    LinearAlgebra,
    NearestNeighbors,
    Shapefile,
    SparseArrays,
    StaticArrays,
    Statistics,
    StatsBase,
    Unitful

# `CoordRefSystems` is a hard dependency, so projection through it is always
# available. This flag is declared before the includes so `lonlat_to_xy_km` /
# `xy_km_to_lonlat` (hexmesh.jl) keep the historical branching shape they had
# in bstm, where it was referenced-but-never-defined (any `crs ≠ nothing` call
# raised `UndefVarError` there). See the notes at the bottom of the module.
const _HAS_COORDSYS = true

include("columns.jl")         # STANDARD_* coordinate/time vocabulary + detectors
include("geometry.jl")        # hulls, polygon primitives, KDE seeding
include("tessellation.jl")    # CVT/KVT/QVT/BVT/HVT/AVT/lattice/hexagonal centroids
include("spatial_units.jl")   # assign_spatial_units dispatcher + granular sizing
include("spatial_graphs.jl")  # weights matrices, knn/radius graphs, BYM2 scaling
include("temporal.jl")        # discretize_data, time/spatiotemporal units, block CV
include("hexmesh.jl")         # lonlat <-> x/y km projection + hex mesh + point mapping
include("hsi.jl")             # HSI per-unit interpolation and field resharding

# Public API. Underscore-prefixed helpers (`_detect_*`, `_jenks_breaks`,
# `_interpolate_hsi`, ...) remain internal-to-the-module, matching the bstm
# conventions they were extracted from.
export
    # Column vocabulary & resolvers (previously exported by bstm)
    STANDARD_SPATIAL_COORDINATE_PAIRS,
    STANDARD_TEMPORAL_CANDIDATES,

    # Spatial unit assignment
    assign_spatial_units,
    assign_spatial_units_inferred,
    assign_time_units,
    assign_spatiotemporal_units,
    discretize_data,
    spatial_block_cv,

    # Geometric primitives & meshes
    expand_hull,
    build_hex_mesh_planar,
    map_point_to_units,
    lonlat_to_xy_km,
    xy_km_to_lonlat,

    # Graph topology & weighting matrices
    spatial_weights_matrix,
    spatial_knn_graph,
    spatial_radius_graph,
    scaling_factor_bym2

end # module GeoPartitioning
