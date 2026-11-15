# Search pipeline, after nhmmer: an ungapped SSV pass over each strand finds
# seeds; windows around the seeds (one model length plus room for insertions
# either side) are scored with local Viterbi; the flanks of each accepted
# alignment are searched again so that nearby copies are not lost.

"A nucleotide profile hit. Coordinates are 1-based inclusive on the original strand."
struct SearchHit
    target_name::String
    target_desc::String
    query_name::String
    target_start::Int
    target_end::Int
    strand::Char
    bit_score::Float64
    evalue::Float64
    seq_length::Int
end
Base.show(io::IO, h::SearchHit) = print(io, "SearchHit(", h.query_name, ", ",
    h.target_name, ":", h.target_start, "-", h.target_end, "(", h.strand, "))")

"Viterbi P-value of a bit score (Gumbel tail, the profile's calibration)."
function viterbi_pvalue(hmm::Plan7HMM, bits)
    mu, lambda = hmm.stats_viterbi
    isfinite(mu) && isfinite(lambda) && lambda > 0 ||
        throw(ArgumentError("profile $(hmm.name) lacks Viterbi calibration statistics"))
    -expm1(-exp(-lambda * (bits - mu)))
end

"""
    search_space(hmm, residues)

E-value multiplier for a search of `residues` bases on both strands: the number
of non-overlapping windows of the profile's maximum alignment length, as nhmmer.
"""
search_space(hmm::Plan7HMM, residues) = max(1.0, 2residues / hmm.maxl)

# Room for net insertions on each side of a seed's diagonal: MAXL allows for
# very long inserts, so cap it; a longer insertion splits the call in two.
_flank(hmm::Plan7HMM) = clamp(hmm.maxl - hmm.M, 32, max(100, hmm.M ÷ 4))

# Merge seed windows on one strand into disjoint ranges, each with its seeds.
function _windows(hmm, seeds, L)
    sort!(seeds); unique!(seeds)
    ext = _flank(hmm)
    out = Tuple{Int,Int,Vector{Int}}[]
    for (r, k) in seeds
        lo = max(1, r - k + 1 - ext); hi = min(L, r + hmm.M - k + ext)
        if !isempty(out) && lo <= out[end][2] + 1
            a, b, rs = out[end]
            push!(rs, r); out[end] = (a, max(b, hi), rs)
        else
            push!(out, (lo, hi, [r]))
        end
    end
    out
end

# Best alignment in xs[lo:hi]; accept it, then search flanks that hold seeds.
function _search_window!(hits, hmm, xs, lo, hi, rs, cutoff, min_length, buf)
    stack = [(lo, hi)]
    while !isempty(stack)
        a0, b0 = pop!(stack)
        any(r -> a0 <= r <= b0, rs) || continue
        w = view(xs, a0:b0)
        # Score-only pass first; only accepted windows pay for endpoint tracking.
        bits = _viterbi_bits(hmm, w, buf)
        bits > 0 && viterbi_pvalue(hmm, bits) <= cutoff || continue
        v = _viterbi(hmm, w, buf)
        v.target_start > 0 || continue
        a = a0 + v.target_start - 1; b = a0 + v.target_end - 1
        b - a + 1 >= min_length && push!(hits, (a, b, v.bit_score))
        a > a0 && push!(stack, (a0, a - 1))
        b < b0 && push!(stack, (b + 1, b0))
    end
    hits
end

struct _Record
    name::String
    desc::String
    seq::LongDNA{4}
end

# Jobs on one strand of one record for one model: SSV over chunks that overlap
# by a model length (no diagonal is longer), then Viterbi per window.
function _ssv_jobs(L, M; chunk=1 << 20)
    L <= chunk + M && return [(1, L, L)]
    [(a, min(L, a + chunk + M - 1), min(chunk, L - a + 1)) for a in 1:chunk:L]
end

_strand_indices(rec, strand) = seq_to_indices(strand == '+' ? rec.seq : reverse_complement(rec.seq))

function _run(f, jobs, threaded)
    if threaded && Threads.nthreads() > 1 && length(jobs) > 1
        tasks = [Threads.@spawn f(j) for j in jobs]
        return map(fetch, tasks)
    end
    map(f, jobs)
end

# Search a batch of records with every model; return hits (unfiltered by E).
function _search_batch(hmms, ssvs, records, cutoffs, min_length, threaded, bar)
    strands = [(ri, s) for ri in eachindex(records) for s in ('+', '-')]
    xs = _run(((ri, s),) -> _strand_indices(records[ri], s), strands, threaded)
    idx = Dict(strands[i] => xs[i] for i in eachindex(strands))
    jobs = [(mi, ri, s, a, b, n) for mi in eachindex(hmms) for (ri, s) in strands
            for (a, b, n) in _ssv_jobs(length(records[ri].seq), hmms[mi].M)]
    seedlists = _run(jobs, threaded) do (mi, ri, s, a, b, n)
        out = _ssv!(Tuple{Int,Int}[], ssvs[mi], view(idx[(ri, s)], a:b), a - 1)
        _next!(bar, n * hmms[mi].M)
        out
    end
    seeds = Dict{Tuple{Int,Int,Char},Vector{Tuple{Int,Int}}}()
    for ((mi, ri, s), sl) in zip(jobs, seedlists)
        append!(get!(seeds, (mi, ri, s), Tuple{Int,Int}[]), sl)
    end
    wjobs = [(mi, ri, s, w) for ((mi, ri, s), sl) in seeds
             for w in _windows(hmms[mi], sl, length(records[ri].seq))]
    found = _run(wjobs, threaded) do (mi, ri, s, (lo, hi, rs))
        _search_window!(Tuple{Int,Int,Float64}[], hmms[mi], idx[(ri, s)], lo, hi, rs,
                        cutoffs[mi], min_length, _ViterbiBuffer(hmms[mi].M))
    end
    hits = Tuple{Int,Int,Char,Int,Int,Float64}[]
    for ((mi, ri, s, _), fs) in zip(wjobs, found), (a, b, bits) in fs
        push!(hits, (mi, ri, s, a, b, bits))
    end
    hits
end

function _hit(hmm, rec, s, a, b, bits, Z)
    L = length(rec.seq)
    start, stop = s == '+' ? (a, b) : (L - b + 1, L - a + 1)
    SearchHit(rec.name, rec.desc, hmm.name, start, stop, s, bits, Z * viterbi_pvalue(hmm, bits), L)
end

_hitorder(h) = (h.target_name, h.target_start, h.strand, h.query_name)

"""
    search_sequence(hmm, seq, name=""; evalue_threshold=1e-5, kwargs...)

Search both strands of one sequence. See [`scan_fasta`](@ref) for options.
"""
function search_sequence(hmm::Plan7HMM, seq, name::AbstractString=""; kwargs...)
    seq = seq isa AbstractString ? LongDNA{4}(seq) : LongDNA{4}(seq)
    _scan([hmm], [_Record(String(name), "", seq)]; kwargs...)
end

function _scan(hmms, records; evalue_threshold=1e-5, seed_bits=default_seed_bits,
               db_size=nothing, min_length=1, threaded=true, bar=_progress(0, ""; enabled=false),
               batch=32_000_000)
    evalue_threshold > 0 || throw(ArgumentError("evalue_threshold must be positive"))
    min_length >= 1 || throw(ArgumentError("min_length must be positive"))
    db_size === nothing || db_size >= 1 || throw(ArgumentError("db_size must be at least one"))
    nres = sum(r -> length(r.seq), records; init=0)
    Zs = [db_size === nothing ? search_space(h, nres) : Float64(db_size) for h in hmms]
    cutoffs = evalue_threshold ./ Zs
    ssvs = [_SSVProfile(h, seed_bits isa Real ? seed_bits : seed_bits(h)) for h in hmms]
    hits = SearchHit[]
    i = 1
    while i <= length(records)       # batches bound the memory held as indices
        j = i; n = 0
        while j <= length(records) && (j == i || n + length(records[j].seq) <= batch)
            n += length(records[j].seq); j += 1
        end
        part = records[i:j-1]
        for (mi, ri, s, a, b, bits) in _search_batch(hmms, ssvs, part, cutoffs, min_length, threaded, bar)
            push!(hits, _hit(hmms[mi], part[ri], s, a, b, bits, Zs[mi]))
        end
        i = j
    end
    sort!(hits; by=_hitorder)
end

function _read_records(fasta)
    records = _Record[]
    seen = Set{String}()
    open(FASTA.Reader, fasta) do reader
        for record in reader
            name = String(identifier(record))
            name in seen && throw(ArgumentError("duplicate FASTA identifier: $name"))
            push!(seen, name)
            push!(records, _Record(name, String(description(record)), LongDNA{4}(sequence(record))))
        end
    end
    records
end

"""
    scan_fasta(hmms, fasta; evalue_threshold=1e-5, seed_bits=default_seed_bits,
               db_size=nothing, min_length=1, threaded=true, progress, io=stderr)

Search both strands of every sequence in `fasta` with one profile or a vector
of profiles. Returns `SearchHit`s sorted by sequence, position and strand.

- `evalue_threshold`: report hits with E-value at most this. E-values are
  Viterbi P-values times [`search_space`](@ref) of the whole file, unless
  `db_size` fixes the multiplier.
- `seed_bits`: SSV seed threshold in bits, a number or a function of the profile.
  Lower values trade speed for sensitivity to weak, gapped hits.
- `threaded`: use all Julia threads.
"""
function scan_fasta(hmms::AbstractVector{Plan7HMM}, fasta; progress::Bool=_progress_default(),
                    io::IO=stderr, kwargs...)
    records = _read_records(fasta)
    nres = sum(r -> length(r.seq), records; init=0)
    bar = _progress(cld(nres, 1000), "Searching"; unit="kb", enabled=progress, io=io)
    # Seed-scan work per residue is proportional to model length.
    scaled = _Scaled(bar, Threads.Atomic{Int}(0), 2000 * sum(h -> h.M, hmms; init=0))
    ok = false; hits = SearchHit[]
    try
        hits = _scan(hmms, records; bar=scaled, kwargs...)
        ok = true
    finally
        _finish!(bar; ok=ok, info=_summary(hmms, hits))
    end
    hits
end
scan_fasta(hmm::Plan7HMM, fasta; kwargs...) = scan_fasta([hmm], fasta; kwargs...)

# Progress in kb of input: work units (residue × model state, both strands)
# divided by `per_kb`.
struct _Scaled
    bar::_Progress
    work::Threads.Atomic{Int}
    per_kb::Int
end
function _next!(s::_Scaled, k::Integer)
    s.bar.enabled || return s
    n = Threads.atomic_add!(s.work, Int(k)) + k
    _update!(s.bar, n ÷ s.per_kb)
    s
end

function _summary(hmms, hits)
    counts = [count(h -> h.query_name == m.name, hits) for m in hmms]
    length(hmms) == 1 && return "$(only(hmms).name) · $(_fmt(length(hits))) hits"
    join(("$(m.name) $(c)" for (m, c) in zip(hmms, counts)), " · ")
end

"""
    nhmmscan(hmm_file, fasta; models=nothing, kwargs...)

Search a FASTA file with all (or the named) profiles of a HMMER3 nucleotide
profile file. Keywords as for [`scan_fasta`](@ref).
"""
function nhmmscan(hmm_file, fasta; models=nothing, kwargs...)
    profiles = parse_hmm_file(hmm_file)
    if models !== nothing
        models isa AbstractString && throw(ArgumentError("models must be a collection of names"))
        unknown = setdiff(models, [p.name for p in profiles])
        isempty(unknown) || throw(ArgumentError("unknown models: $(join(unknown, ", "))"))
        filter!(p -> p.name in models, profiles)
    end
    scan_fasta(Plan7HMM[build_plan7(p) for p in profiles], fasta; kwargs...)
end
