# -----------------------------------------------------------------------------
# geometry.jl -- (GeoPartitioning) extracted verbatim from bstm src/partitioning.jl
# lines 1-282 (Section 1: Geometric Primitives and Boundary Utilities).
# ------------------------------------------------------------------------------
"""
    partitioning.jl

Spatial and spatiotemporal domain decomposition, Voronoi/hexagonal tessellation,
areal unit assignment, graph connectivity, and spatial block cross-validation for BSTM.

Version: v1.0.0
"""

# -----------------------------------------------------------------------------
# Section 1: Geometric Primitives & Boundary Utilities
# -----------------------------------------------------------------------------

"""
    expand_hull(s_x::AbstractVector{<:Real}, s_y::AbstractVector{<:Real}, buffer_dist::Real)

Computes the convex hull of 2D points `(s_x, s_y)` and expands it by `buffer_dist`
using `LibGEOS.buffer`. Returns a `LibGEOS` polygon geometry.
"""
function expand_hull(s_x::AbstractVector{<:Real}, s_y::AbstractVector{<:Real}, buffer_dist::Real)
    if isempty(s_x) || isempty(s_y)
        return LibGEOS.Polygon([[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]])
    end
    coords_vec = [[Float64(s_x[i]), Float64(s_y[i])] for i in 1:length(s_x)]
    points_geom = LibGEOS.MultiPoint(coords_vec)
    hull = LibGEOS.convexhull(points_geom)
    if buffer_dist > 0.0
        return LibGEOS.buffer(hull, Float64(buffer_dist))
    else
        return hull
    end
end

expand_hull(coords::Vector{<:Tuple{Real, Real}}, buffer_dist::Real) =
    expand_hull([p[1] for p in coords], [p[2] for p in coords], buffer_dist)


"""
    get_polygon_area(poly_coords::AbstractVector)

Calculates the 2D surface area of a polygon using the Shoelace formula.
Coordinates are expected as a vector of `(x, y)` tuples or vectors.
"""
function get_polygon_area(poly_coords::AbstractVector)
    valid_pts = Tuple{Float64, Float64}[]
    for p in poly_coords
        x, y = float(p[1]), float(p[2])
        if !isnan(x) && !isinf(x) && !isnan(y) && !isinf(y)
            push!(valid_pts, (x, y))
        end
    end

    if length(valid_pts) > 1 && valid_pts[1] == valid_pts[end]
        valid_pts = valid_pts[1:end-1]
    end

    if length(valid_pts) < 3
        return 0.0
    end

    x = [p[1] for p in valid_pts]
    y = [p[2] for p in valid_pts]

    # Shoelace formula
    return 0.5 * abs(dot(x, circshift(y, 1)) - dot(y, circshift(x, 1)))
end

get_polygon_area(s_x::AbstractVector, s_y::AbstractVector) = get_polygon_area(tuple.(s_x, s_y))


"""
    is_valid_polygon_coords(poly_coords)

Validates whether a collection of vertices contains at least 3 non-NaN/non-Inf points.
"""
function is_valid_polygon_coords(poly_coords)
    valid_pts = [
        p for p in poly_coords
        if !isnan(p[1]) && !isinf(p[1]) && !isnan(p[2]) && !isinf(p[2])
    ]
    return length(valid_pts) >= 3
end


"""
    get_coords_from_geom(geom)

Extracts coordinate tuples `Vector{Tuple{Float64, Float64}}` from `LibGEOS` geometries
(`Point`, `Polygon`, `MultiPolygon`, `LineString`, `LinearRing`).
"""
function get_coords_from_geom(geom)
    coords = Tuple{Float64, Float64}[]
    if geom === nothing
        return coords
    end
    local type_id = -1
    try
        type_id = LibGEOS.geomTypeId(geom)
        if type_id == LibGEOS.GEOS_POINT
            seq = LibGEOS.getCoordSeq(geom)
            push!(coords, (LibGEOS.getX(seq, 1), LibGEOS.getY(seq, 1)))
            return coords
        elseif type_id == LibGEOS.GEOS_POLYGON
            ring = LibGEOS.exteriorRing(geom)
            n = LibGEOS.numPoints(ring)
            for i in 1:n
                p = LibGEOS.getPoint(ring, i)
                seq = LibGEOS.getCoordSeq(p)
                push!(coords, (LibGEOS.getX(seq, 1), LibGEOS.getY(seq, 1)))
            end
        elseif type_id == LibGEOS.GEOS_MULTIPOLYGON
            n_geoms = LibGEOS.numGeometries(geom)
            for i in 1:n_geoms
                poly = LibGEOS.getGeometryN(geom, i)
                ring = LibGEOS.exteriorRing(poly)
                n = LibGEOS.numPoints(ring)
                for j in 1:n
                    p = LibGEOS.getPoint(ring, j)
                    seq = LibGEOS.getCoordSeq(p)
                    push!(coords, (LibGEOS.getX(seq, 1), LibGEOS.getY(seq, 1)))
                end
                if i < n_geoms
                    push!(coords, (NaN, NaN))
                end
            end
        elseif type_id in [LibGEOS.GEOS_LINESTRING, LibGEOS.GEOS_LINEARRING]
            n = LibGEOS.numPoints(geom)
            for i in 1:n
                p = LibGEOS.getPoint(geom, i)
                seq = LibGEOS.getCoordSeq(p)
                push!(coords, (LibGEOS.getX(seq, 1), LibGEOS.getY(seq, 1)))
            end
        end
    catch e
        @warn "Coordinate extraction failed for LibGEOS type $type_id: $e"
    end
    return coords
end


"""
    get_kde_seeds(coords, target_u)

Fast density-weighted seeding for spatial tessellation. Uses KDTree to estimate local
point densities in O(N log N) without allocating large pairwise distance matrices.
"""
function get_kde_seeds(coords::Vector{<:Tuple{Real, Real}}, target_u::Integer)
    u_pts = unique(coords)
    n = length(u_pts)
    if n == 0
        return Tuple{Float64, Float64}[]
    elseif n <= target_u
        return [(Float64(p[1]), Float64(p[2])) for p in u_pts]
    end

    # Fast KD-Tree 5-nearest-neighbors density proxy
    pts_mat = hcat([[Float64(p[1]), Float64(p[2])] for p in u_pts]...)
    tree = KDTree(pts_mat)
    k_eval = min(6, n)
    _, dists = knn(tree, pts_mat, k_eval)

    # Inverse mean distance as density weight
    weights = [1.0 / (mean(d[2:end]) + 1e-6) for d in dists]
    idx = StatsBase.sample(1:n, Weights(weights), min(target_u, n), replace=false)
    return [(Float64(u_pts[i][1]), Float64(u_pts[i][2])) for i in idx]
end


"""
    get_voronoi_polygons_and_edges(centroids, hull_geom, tol=1e-7)

Generates Voronoi polygons for `centroids`, clips them to `hull_geom` via `LibGEOS`,
and extracts contiguous neighbor edges. Handles edge cases for 0, 1, and 2 centroids.
"""
function get_voronoi_polygons_and_edges(
    centroids::Vector{<:Tuple{Real, Real}}, hull_geom, tol::Real=1e-7
)
    n_c = length(centroids)
    if n_c == 0
        return Tuple{Float64, Float64}[][],
               Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}[]
    elseif n_c == 1
        return [get_coords_from_geom(hull_geom)],
               Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}[]
    elseif n_c == 2
        p1 = (Float64(centroids[1][1]), Float64(centroids[1][2]))
        p2 = (Float64(centroids[2][1]), Float64(centroids[2][2]))
        mid = ((p1[1] + p2[1]) / 2.0, (p1[2] + p2[2]) / 2.0)
        dx, dy = p2[1] - p1[1], p2[2] - p1[2]
        px, py = -dy, dx

        env = LibGEOS.envelope(hull_geom)
        hull_bbox_coords = get_coords_from_geom(env)
        L = if !isempty(hull_bbox_coords)
            min_x, max_x = extrema(c[1] for c in hull_bbox_coords)
            min_y, max_y = extrema(c[2] for c in hull_bbox_coords)
            2.0 * sqrt((max_x - min_x)^2 + (max_y - min_y)^2) + 1.0
        else
            1e7
        end

        pt1 = (mid[1] + L*px, mid[2] + L*py)
        pt2 = (mid[1] - L*px, mid[2] - L*py)
        side1_pts = [pt1, pt2, (pt2[1] - L*dx, pt2[2] - L*dy), (pt1[1] - L*dx, pt1[2] - L*dy), pt1]
        poly1_box = LibGEOS.Polygon([[[p[1], p[2]] for p in side1_pts]])
        side2_pts = [pt1, pt2, (pt2[1] + L*dx, pt2[2] + L*dy), (pt1[1] + L*dx, pt1[2] + L*dy), pt1]
        poly2_box = LibGEOS.Polygon([[[p[1], p[2]] for p in side2_pts]])
        res1 = LibGEOS.intersection(hull_geom, poly1_box)
        res2 = LibGEOS.intersection(hull_geom, poly2_box)
        return [get_coords_from_geom(res1), get_coords_from_geom(res2)], [(p1, p2)]
    end

    # Deduplicate centroids
    u_centroids = unique(centroids)
    if length(u_centroids) < n_c
        u_polys, u_edges = get_voronoi_polygons_and_edges(u_centroids, hull_geom, tol)
        return [u_polys[findfirst(==(c), u_centroids)] for c in centroids], u_edges
    end

    # 3+ centroids: Delaunay -> Voronoi
    pts_dt = [(Float64(c[1]), Float64(c[2])) for c in centroids]
    tri = triangulate(pts_dt)
    hull_coords = get_coords_from_geom(hull_geom)
    xs = [p[1] for p in hull_coords if !isnan(p[1])]
    ys = [p[2] for p in hull_coords if !isnan(p[2])]
    if isempty(xs) || isempty(ys)
        return [Tuple{Float64, Float64}[] for _ in 1:length(centroids)], []
    end
    
    bbox = (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
    vorn = voronoi(tri)
    final_coords = [Tuple{Float64, Float64}[] for _ in 1:length(centroids)]
    valid_geoms = Dict{Int, Any}()

    for i in each_generator(vorn)
        if i < 1 || i > length(centroids)
            continue
        end
        vertices = get_polygon_coordinates(vorn, i, bbox)
        if !isempty(vertices)
            poly_pts = [[v[1], v[2]] for v in vertices]
            if poly_pts[1] != poly_pts[end]
                push!(poly_pts, poly_pts[1])
            end
            try
                lg_poly = LibGEOS.Polygon([poly_pts])
                clipped = LibGEOS.intersection(lg_poly, hull_geom)
                if !LibGEOS.isEmpty(clipped) &&
   LibGEOS.geomTypeId(clipped) in [LibGEOS.GEOS_POLYGON, LibGEOS.GEOS_MULTIPOLYGON]
                    final_coords[i] = get_coords_from_geom(clipped)
                    valid_geoms[i] = clipped
                end
            catch
            end
        end
    end

    v_edges = Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}[]
    active_ids = sort(collect(keys(valid_geoms)))
    for idx in 1:length(active_ids)
        i = active_ids[idx]
        g1 = valid_geoms[i]
        for jdx in idx+1:length(active_ids)
            j = active_ids[jdx]
            g2 = valid_geoms[j]
            if LibGEOS.touches(g1, g2)
                push!(v_edges, ((Float64(centroids[i][1]), Float64(centroids[i][2])),
                                (Float64(centroids[j][1]), Float64(centroids[j][2]))))
            else
                g1_b = LibGEOS.buffer(g1, Float64(tol))
                if LibGEOS.intersects(g1_b, g2)
                    push!(v_edges, ((Float64(centroids[i][1]), Float64(centroids[i][2])),
                                    (Float64(centroids[j][1]), Float64(centroids[j][2]))))
                end
            end
        end
    end

    return final_coords, v_edges
end


