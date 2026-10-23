# Environment setup for GeoPartitioning: resolves the [deps] of Project.toml
# against the installed registries and builds Manifest.toml in-place.
# Usage: julia --project=. scripts/setup_env.jl
# Deliberately minimal: no Pkg.add-by-name, no test deps pulled in here
# ([extras] entries are handled by `Pkg.test()` / `julia --project=. test/runtests.jl`).

using Pkg

Pkg.instantiate()
println("GeoPartitioning environment instantiated successfully.")
