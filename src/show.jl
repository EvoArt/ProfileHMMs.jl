Base.show(io::IO, h::Plan7HMM) = print(io, "Plan7HMM(\"", h.name, "\", ", h.M, " states)")
function Base.show(io::IO, ::MIME"text/plain", h::Union{Plan7HMM,HMMProfile})
    _title(io, string(nameof(typeof(h))), h.name)
    _kv(io, "states", h isa Plan7HMM ? h.M : h.length)
    _kv(io, "alphabet", h isa Plan7HMM ? "nucleotide" : h.alphabet)
    _kv(io, "scoring", "MSV / single-hit Viterbi")
end
function Base.show(io::IO, ::MIME"text/plain", h::SearchHit)
    _title(io, "SearchHit", h.query_name)
    _kv(io, "target", h.target_name)
    _kv(io, "coordinates", "$(h.target_start):$(h.target_end) ($(h.strand))")
    _kv(io, "score", "$(_g(h.bit_score)) bits")
    _kv(io, "approx. E", _g(h.evalue))
end
