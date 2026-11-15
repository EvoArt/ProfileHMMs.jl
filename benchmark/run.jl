# Time ProfileHMMs against nhmmer on the same profiles, genomes, E-value
# cutoff and thread count, and check that they find the same hits.
#
#   julia -t T --project=. benchmark/run.jl OUTDIR MODEL.hmm GENOME.fasta...
#
# Environment: NHMMERS="label=/path/to/nhmmer,label2=/path/to/other"
#              BENCH_REPEATS (default 3)
# nhmmer runs with --cpu T, or --cpu 0 (fully serial) when T = 1.
using ProfileHMMs, TOML
const P = ProfileHMMs

const EVALUE = 1e-5
med(x) = (s = sort(x); n = length(s); isodd(n) ? s[cld(n, 2)] : (s[n÷2] + s[n÷2+1]) / 2)

function parse_tbl(path)
    rows = NamedTuple[]
    for line in eachline(path)
        (startswith(line, "#") || isempty(strip(line))) && continue
        p = split(line); a, b = parse.(Int, p[7:8])
        push!(rows, (target=String(p[1]), model=String(p[3]), strand=only(p[12]),
                     start=min(a, b), stop=max(a, b), bits=parse(Float64, p[14])))
    end
    rows
end

# nhmmer's "Residues passing ... filter" fractions, per query.
function pipeline_stats(path)
    stats = Dict{String,Dict{String,Float64}}(); q = ""
    for line in eachline(path)
        m = match(r"^Query:\s+(\S+)", line)
        m === nothing || (q = m[1]; stats[q] = Dict{String,Float64}())
        m = match(r"^Residues passing (\S+) filter:\s+\d+\s+\(([0-9.e-]+)\)", line)
        m === nothing || (stats[q][m[1]] = parse(Float64, m[2]))
    end
    stats
end

# One-to-one matching on sequence, model and strand; endpoints within 10 bp.
function agreement(hits, ref)
    left = collect(eachindex(hits)); err = 0; missing = 0
    for r in ref
        c = filter(i -> (hits[i].target_name, hits[i].query_name, hits[i].strand) == (r.target, r.model, r.strand), left)
        e = [max(abs(hits[i].target_start - r.start), abs(hits[i].target_end - r.stop)) for i in c]
        if isempty(c) || minimum(e) > 10
            missing += 1
        else
            k = argmin(e); err = max(err, e[k]); filter!(!=(c[k]), left)
        end
    end
    (matched=length(ref) - missing, missing=missing, extra=length(left), max_endpoint_error=err)
end

# Fraction of residues (both strands) inside ProfileHMMs' seeded windows: the
# counterpart of nhmmer's "passing SSV filter".
function seed_coverage(hmm, fasta)
    sp = P._SSVProfile(hmm, default_seed_bits(hmm)); cov = 0; tot = 0
    for rec in P._read_records(fasta), s in ('+', '-')
        xs = P._strand_indices(rec, s)
        seeds = P._ssv!(Tuple{Int,Int}[], sp, view(xs, 1:length(xs)), 0)
        cov += sum(w -> w[2] - w[1] + 1, P._windows(hmm, seeds, length(xs)); init=0)
        tot += length(xs)
    end
    cov / tot
end

function main()
    out, model = ARGS[1], ARGS[2]; genomes = ARGS[3:end]
    mkpath(out)
    T = Threads.nthreads()
    repeats = parse(Int, get(ENV, "BENCH_REPEATS", "3"))
    nhmmers = [Tuple(split(x, "=")) for x in split(get(ENV, "NHMMERS", ""), ",") if !isempty(x)]
    hmms = Plan7HMM[build_plan7(p) for p in parse_hmm_file(model)]
    rows = Dict{String,Any}[]
    for g in genomes
        name = splitext(basename(g))[1]
        first_call = @elapsed scan_fasta(hmms, g; progress=false, evalue_threshold=EVALUE)
        vit = scan_fasta(hmms, g; progress=false, evalue_threshold=EVALUE)
        for mode in (:viterbi, :forward)
            ts = [@elapsed(scan_fasta(hmms, g; progress=false, evalue_threshold=EVALUE, scoring=mode)) for _ in 1:repeats]
            push!(rows, Dict("genome" => name, "threads" => T, "tool" => "ProfileHMMs $mode",
                             "seconds" => med(ts), "first_call_seconds" => mode == :viterbi ? first_call : 0.0,
                             "hits" => length(vit),
                             "first_filter_fraction" => Dict(h.name => seed_coverage(h, g) for h in hmms)))
        end
        for (label, exe) in nhmmers
            tbl = joinpath(out, "$name-$label-t$T.tbl"); txt = joinpath(out, "$name-$label-t$T.txt")
            cmd = `$exe --cpu $(T == 1 ? 0 : T) --noali -E $EVALUE --tblout $tbl $model $g`
            ts = [@elapsed(run(pipeline(cmd; stdout=txt))) for _ in 1:repeats]
            ref = parse_tbl(tbl); st = pipeline_stats(txt)
            push!(rows, Dict("genome" => name, "threads" => T, "tool" => "nhmmer $label",
                             "seconds" => med(ts), "hits" => length(ref),
                             "agreement" => Dict(string(k) => v for (k, v) in pairs(agreement(vit, ref))),
                             "first_filter_fraction" => Dict(q => get(d, "SSV", NaN) for (q, d) in st),
                             "viterbi_filter_fraction" => Dict(q => get(d, "Vit", NaN) for (q, d) in st)))
        end
        println(name, " done (", T, " threads)"); flush(stdout)
    end
    open(joinpath(out, "results-t$T.toml"), "w") do io
        TOML.print(io, Dict("julia" => string(VERSION), "cpu" => Sys.CPU_NAME, "rows" => rows))
    end
end
main()
