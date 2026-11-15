"A nucleotide profile hit. Coordinates are 1-based inclusive on the original strand."
struct SearchHit
    target_name::String
    target_desc::String
    query_name::String
    target_start::Int
    target_end::Int
    strand::Char
    bit_score::Float64
    evalue::Float64 # Approximate Viterbi tail probability * db_size, not nhmmer calibrated.
    seq_length::Int
end
Base.show(io::IO, h::SearchHit) = print(io, "SearchHit(", h.query_name, ", ",
    h.target_name, ":", h.target_start, "-", h.target_end, "(", h.strand, "))")

"""
    search_sequence(hmm, seq, name=""; msv_threshold=0.02, evalue_threshold=1e-5,
                    db_size=1, min_length=1)

Search both strands; iteratively search the flanks of each accepted hit to find
multiple non-overlapping copies. E-values are approximate Viterbi statistics.
Set `msv_threshold=Inf` to disable the heuristic ungapped filter.
"""
function search_sequence(hmm::Plan7HMM, seq, name="";
        msv_threshold=0.02, evalue_threshold=1e-5, db_size=1, min_length=1)
    msv_threshold >= 0 || throw(ArgumentError("msv_threshold must be nonnegative"))
    evalue_threshold > 0 || throw(ArgumentError("evalue_threshold must be positive"))
    db_size >= 1 || throw(ArgumentError("db_size must be at least one"))
    min_length >= 1 || throw(ArgumentError("min_length must be positive"))
    seq = seq isa AbstractString ? LongDNA{4}(seq) : seq
    hits = SearchHit[]
    L = length(seq)
    for (strand_seq, strand) in ((seq, '+'), (reverse_complement(seq), '-'))
        indices = seq_to_indices(strand_seq)
        pending = [(1, L)]
        while !isempty(pending)
            lo, hi = pop!(pending)
            hi - lo + 1 >= min_length || continue
            xs = view(indices, lo:hi)
            if isfinite(msv_threshold)
                msv_evalue(hmm, msv_score(hmm, xs), length(xs)) <= msv_threshold || continue
            end
            vit = viterbi_score(hmm, xs)
            ev = db_size * _tail_probability(vit.bit_score, hmm.stats_viterbi)
            vit.bit_score > 0 && ev <= evalue_threshold && vit.target_start > 0 || continue
            a = lo + vit.target_start - 1; b = lo + vit.target_end - 1
            if b - a + 1 >= min_length
                start, stop = strand == '+' ? (a,b) : (L-b+1,L-a+1)
                push!(hits, SearchHit(String(name), "", hmm.name, start, stop, strand,
                                     vit.bit_score, ev, L))
            end
            a > lo && push!(pending, (lo, a-1))
            b < hi && push!(pending, (b+1, hi))
        end
    end
    sort!(hits; by=h -> (h.target_start, h.strand))
end

# Suppress duplicate calls from overlapping windows without joining independent
# loci or extending a selected alignment's endpoints.
function _merge_hits(hits; reciprocal_overlap=0.5)
    kept = SearchHit[]
    for h in sort(hits; by=h -> (-h.bit_score, h.target_start, h.target_end))
        duplicate = any(kept) do k
            (h.target_name, h.query_name, h.strand) == (k.target_name, k.query_name, k.strand) || return false
            overlap = max(0, min(h.target_end,k.target_end)-max(h.target_start,k.target_start)+1)
            overlap / min(h.target_end-h.target_start+1,k.target_end-k.target_start+1) >= reciprocal_overlap
        end
        duplicate || push!(kept,h)
    end
    sort!(kept; by=h -> (h.target_name, h.target_start, h.strand, h.query_name))
end

"Scan a FASTA in overlapping windows; each window owns its results, safe across thread pools."
function scan_fasta(hmm::Plan7HMM, fasta; window=max(5000, 3hmm.M),
        overlap=max(2000, 2hmm.M), progress::Bool=_progress_default(), io::IO=stderr,
        threaded::Bool=true, kwargs...)
    window > 0 && 0 <= overlap < window || throw(ArgumentError("require window > overlap >= 0"))
    all_hits = SearchHit[]
    bar = _progress(0, "Searching"; unit="windows", enabled=progress, io=io)
    ok = false
    seen = Set{String}()
    try
        open(FASTA.Reader, fasta) do reader
            for record in reader
                seq = LongDNA{4}(sequence(record)); L = length(seq)
                name = String(identifier(record))
                name in seen && throw(ArgumentError("duplicate FASTA identifier: $name"))
                push!(seen,name)
                L == 0 && continue
                starts = L <= window ? [1] : collect(1:(window-overlap):(L-window+1))
                L > window && starts[end] != L-window+1 && push!(starts,L-window+1)
                chunks = Vector{Vector{SearchHit}}(undef,length(starts))
                function work(wi)
                    a = starts[wi]; b = min(a+window-1,L)
                    hits = search_sequence(hmm, seq[a:b], name; kwargs...)
                    chunks[wi] = [SearchHit(h.target_name, String(description(record)),h.query_name,
                        h.target_start+a-1,h.target_end+a-1,h.strand,h.bit_score,h.evalue,L) for h in hits]
                    _next!(bar)
                end
                if threaded && length(starts) > 1
                    Threads.@threads for wi in eachindex(starts)
                        work(wi)
                    end
                else
                    foreach(work, eachindex(starts))
                end
                append!(all_hits, _merge_hits(reduce(vcat,chunks)))
            end
        end
        ok = true
    finally
        _finish!(bar; ok=ok, info="$(length(all_hits)) hits")
    end
    sort!(all_hits; by=h -> (h.target_name,h.target_start,h.strand))
end

"""
    nhmmscan(hmm_file, fasta; models=nothing, target_length=2000, kwargs...)

Search all or named nucleotide profiles. This is a Julia Viterbi pipeline,
not an implementation of HMMER's complete nhmmer/nhmmscan pipeline.
"""
function nhmmscan(hmm_file, fasta; models=nothing, target_length=2000, kwargs...)
    profiles = parse_hmm_file(hmm_file)
    if models !== nothing
        models isa AbstractString && throw(ArgumentError("models must be a collection of names"))
        missing = setdiff(models, [p.name for p in profiles])
        isempty(missing) || throw(ArgumentError("unknown models: $(join(missing, ", "))"))
        filter!(p -> p.name in models, profiles)
    end
    hits = SearchHit[]
    for p in profiles
        append!(hits, scan_fasta(build_plan7(p; target_length=target_length), fasta; kwargs...))
    end
    sort!(hits; by=h -> (h.target_name,h.target_start,h.strand,h.query_name))
end