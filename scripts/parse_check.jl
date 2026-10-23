# Parse-only sanity check for every source/test file in GeoPartitioning.
# Usage: julia --project=. scripts/parse_check.jl
# Exits non-zero on the first malformed file; prints file + line + reason.

dirs = ["src", "test"]
ok = true
for d in dirs
    for f in sort(readdir(joinpath(@__DIR__, "..", d); join = true))
        endswith(f, ".jl") || continue
        code = read(f, String)
        result = Meta.parseall(code; filename = f)
        for ex in result.args
            if ex isa Expr && ex.head == :incomplete
                println(stderr, "PARSE FAIL: ", f)
                println(stderr, "    ", ex.args[1])
                global ok = false
                break
            elseif ex isa Expr && ex.head == :error
                println(stderr, "PARSE FAIL: ", f)
                println(stderr, "    ", ex.args[1])
                global ok = false
                break
            elseif ex isa Expr && ex.head == :toplevel
                for sub in ex.args
                    if sub isa Expr && sub.head in (:incomplete, :error)
                        println(stderr, "PARSE FAIL: ", f)
                        println(stderr, "    ", sub.args[1])
                        global ok = false
                        break
                    end
                end
            end
        end
        ok || exit(1)
        println("parse ok: ", f)
    end
end
exit(ok ? 0 : 1)
