# Markdown tables from benchmark/run.jl output:  julia benchmark/summary.jl OUTDIR
using TOML

rows = Dict{String,Any}[]; cpu = ""
for f in sort(filter(endswith(".toml"), readdir(ARGS[1]; join=true)))
    t = TOML.parsefile(f); global cpu = t["cpu"]; append!(rows, t["rows"])
end
genomes = unique(r["genome"] for r in rows)
tools = unique(r["tool"] for r in rows)
fmt(x) = string(round(x; digits=2))

println("CPU: `", cpu, "`. Median of repeated runs, whole-genome search with all three profiles.\n")
for T in sort(unique(r["threads"] for r in rows))
    println("**", T, T == 1 ? " thread" : " threads", "** (seconds)\n")
    println("| genome | ", join(tools, " | "), " |")
    println("|---", "|---"^length(tools), "|")
    for g in genomes
        cells = map(tools) do tool
            r = findfirst(r -> (r["genome"], r["tool"], r["threads"]) == (g, tool, T), rows)
            r === nothing ? "–" : fmt(rows[r]["seconds"])
        end
        println("| ", g, " | ", join(cells, " | "), " |")
    end
    println()
end

println("**Agreement** with ProfileHMMs' hits (matched / missing / extra, largest endpoint difference)\n")
nh = filter(startswith("nhmmer"), tools)
println("| genome | ", join(nh, " | "), " |")
println("|---", "|---"^length(nh), "|")
for g in genomes
    cells = map(nh) do tool
        r = findfirst(r -> (r["genome"], r["tool"]) == (g, tool), rows)
        r === nothing && return "–"
        a = rows[r]["agreement"]
        "$(a["matched"]) / $(a["missing"]) / $(a["extra"]) ($(a["max_endpoint_error"]) bp)"
    end
    println("| ", g, " | ", join(cells, " | "), " |")
end
println()

println("**Work**: fraction of residues (both strands) passed on by the first filter, mean over genomes\n")
models = sort(collect(keys(first(rows)["first_filter_fraction"])))
println("| tool | ", join(models, " | "), " |")
println("|---", "|---"^length(models), "|")
for tool in unique(replace(t, r" (viterbi|forward)$" => "") for t in tools)
    rs = filter(r -> startswith(r["tool"], tool) && r["threads"] == 1, rows)
    isempty(rs) && continue
    cells = [string(round(100 * sum(r["first_filter_fraction"][m] for r in rs) / length(rs); digits=2), "%") for m in models]
    println("| ", tool, " | ", join(cells, " | "), " |")
end
