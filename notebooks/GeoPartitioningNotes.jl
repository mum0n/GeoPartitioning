### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
begin
    using Pkg
    Pkg.activate(@__DIR__)
    using GeoPartitioning
    using DataFrames
    using Random
    using Statistics
    using Plots
    default(size=(700,450), legend=false)
    "Environment active & GeoPartitioning loaded."
end

# ╔═╡ 00000000-0000-0000-0000-000000000002
md"""
# GeoPartitioning — spatial & spatiotemporal domain decomposition

This notebook walks through the main workflows of `GeoPartitioning`:

1. **Spatial partitioning** of point data into areal units (hexagonal, CVT, K-means, …)
2. **Spatial graph topology** (k-NN, radius graphs, weighting matrices)
3. **Temporal discretisation** of a time variable
4. **Spatiotemporal synchronisation** (`st = (t-1)S + s`)
5. **Hexagonal meshing** on a tangent plane

Each section generates a synthetic dataset, runs a partitioning method, and visualises the result.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000003
md"""
## 1. Hexagonal spatial partitioning

Generate random point observations and tessellate the domain into hexagonal areal units.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000004
begin
    rng = MersenneTwister(42)
    N = 400
    df = DataFrame(
        lon = rand(rng, N) .* 100.0,
        lat = rand(rng, N) .* 50.0,
    )
    au = assign_spatial_units(df; x=:lon, y=:lat, area_method=:hexagonal, target_units=16, exact_units=true)
    df.s_idx = au.assignments
    (n_units = au.n_units, n_points = N, cv_area = au.metrics.cv_area)
end

# ╔═╡ 00000000-0000-0000-0000-000000000005
function plot_partition(au, df; title="")
    p = scatter(df.lon, df.lat, group=df.s_idx, ms=3, alpha=0.7, title=title)
    for poly in au.polygons
        xs = [pt[1] for pt in poly]
        ys = [pt[2] for pt in poly]
        push!(xs, poly[1][1]); push!(ys, poly[1][2])
        plot!(p, xs, ys, color=:black, lw=0.8, label="")
    end
    for c in au.centroids
        scatter!(p, [c[1]], [c[2]], color=:red, ms=4, label="")
    end
    p
end

plot_partition(au, df; title="Hexagonal partitioning ($(au.n_units) units)")

# ╔═╡ 00000000-0000-0000-0000-000000000006
md"""
## 2. K-means / CVT partitioning

Compare a density-balanced K-means-seeded tessellation (`:kvt`) with a centroidal Voronoi tessellation (`:cvt`).
"""

# ╔═╡ 00000000-0000-0000-0000-000000000007
begin
    au_kvt = assign_spatial_units(df; x=:lon, y=:lat, area_method=:kvt, target_units=16)
    au_cvt = assign_spatial_units(df; x=:lon, y=:lat, area_method=:cvt, target_units=16)
    p1 = plot_partition(au_kvt, df; title="K-means seeded (:kvt)")
    p2 = plot_partition(au_cvt, df; title="Centroidal Voronoi (:cvt)")
    plot(p1, p2, layout=(1,2), size=(900,400))
end

# ╔═╡ 00000000-0000-0000-0000-000000000008
md"""
## 3. Spatial graph topology

Build a k-nearest-neighbour graph and a radius graph from the unit centroids, and inspect a row-standardised weights matrix.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000009
begin
    coords = au.centroids
    g_knn, W_knn = spatial_knn_graph(coords, 4)
    g_rad, W_rad = spatial_radius_graph(coords, 15.0)
    W_std = spatial_weights_matrix(W_knn; style=:row_standardized)
    (knn_edges = count(!iszero, W_knn) ÷ 2, radius_edges = count(!iszero, W_rad) ÷ 2,
     row_sums = round.(sum(W_std, dims=2)[1:3], digits=3))
end

# ╔═╡ 00000000-0000-0000-0000-000000000010
md"""
## 4. Temporal discretisation

Bin a continuous time variable into `t_N` regular intervals.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000011
begin
    df2 = DataFrame(
        lon = df.lon[1:N],
        lat = df.lat[1:N],
        year = rand(rng, 2018:2024, N),
    )
    au_t = assign_time_units(df2.year; time_method="unique", t_N=6)
    df2.t_idx = au_t.idx
    (n_time_bins = au_t.N_cat, first_bins = first(df2.t_idx, 10))
end

# ╔═╡ 00000000-0000-0000-0000-000000000012
md"""
## 5. Spatiotemporal synchronisation

Combine spatial and temporal units into a single index `st = (t-1)S + s`.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000013
begin
    st = assign_spatiotemporal_units(df2;
        space_x=:lon, space_y=:lat, time_var=:year,
        area_method=:hexagonal, target_units=12, time_method="unique", t_N=4)
    (S = length(st.au_spatial.centroids), T = st.au_temporal.N_cat,
     n_st_units = maximum(st.st_idx), first_st = first(st.st_idx, 10))
end

# ╔═╡ 00000000-0000-0000-0000-000000000014
md"""
## 6. Hexagonal mesh on a tangent plane

Build a planar hexagonal mesh and map points to the nearest hexagon centroid.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000015
begin
    hex = build_hex_mesh_planar(df.lon, df.lat; radius_km=1000.0)
    p = plot(legend=false, title="Planar hex mesh ($(hex.n_units) hexes)")
    for poly in hex.polygons_lonlat
        xs = [pt[1] for pt in poly]; ys = [pt[2] for pt in poly]
        push!(xs, poly[1][1]); push!(ys, poly[1][2])
        plot!(p, xs, ys, color=:black, lw=0.8, label="")
    end
    scatter!(p, [c[1] for c in hex.centroids_lonlat], [c[2] for c in hex.centroids_lonlat], color=:red, ms=3)
    p
end
