# Search profile in HMMER3 local mode, as log-odds scores (nats, Float32,
# higher is better). As in HMMER, insert emissions score zero, so a profile's
# insert-state composition cannot reward biased target composition, and local
# entry into match state k is weighted by the occupancy of k.

"""
    Plan7HMM

A nucleotide profile configured for local search. Build one with
[`build_plan7`](@ref). Scores are natural-log odds against a uniform background.
"""
struct Plan7HMM
    name::String
    M::Int
    maxl::Int                   # MAXL: longest expected alignment (sets search windows)
    msc::Matrix{Float32}        # M × 5 match log-odds; column = base (A, C, G, T, ambiguous = 0)
    tmm::Vector{Float32}        # into M_k from M_{k-1}, I_{k-1}, D_{k-1} (index k; k=1 unused)
    tim::Vector{Float32}
    tdm::Vector{Float32}
    tmi::Vector{Float32}        # out of node k: M_k→I_k, I_k→I_k, M_k→D_{k+1}, D_k→D_{k+1}
    tii::Vector{Float32}
    tmd::Vector{Float32}
    tdd::Vector{Float32}
    tbm::Vector{Float32}        # local entry B→M_k
    stats_msv::Tuple{Float64,Float64}
    stats_viterbi::Tuple{Float64,Float64}
    stats_forward::Tuple{Float64,Float64}
end

"Map a nucleotide to 1-4 (A, C, G, T/U); 0 for an ambiguous base."
@inline function base_index(nt)
    nt == DNA_A && return 0x01
    nt == DNA_C && return 0x02
    nt == DNA_G && return 0x03
    nt == DNA_T && return 0x04
    return 0x00
end

"Encode a DNA sequence as base indices (`UInt8`, 1-4, 0 = ambiguous)."
seq_to_indices(seq) = UInt8[base_index(nt) for nt in seq]

_lp(x) = Float32(-x)   # file values are -ln(p); "*" was parsed as Inf

# Expected occupancy of each match state (HMMER's p7_hmm_CalculateOccupancy).
function _occupancy(t)
    M = size(t, 2) - 1
    p(r, k) = exp(-t[r, k+1])          # transition row r out of node k
    occ = zeros(M)
    occ[1] = p(2, 0) + p(1, 0)         # 1 - P(B→D_1)
    for k in 2:M
        occ[k] = occ[k-1] * (p(1, k-1) + p(2, k-1)) + (1 - occ[k-1]) * p(6, k-1)
    end
    occ
end

"""
    build_plan7(profile::HMMProfile) -> Plan7HMM

Configure a parsed profile for local, single-hit nucleotide search.
"""
function build_plan7(profile::HMMProfile)
    M = profile.length
    M > 0 || throw(ArgumentError("profile length must be positive"))
    t = profile.transitions
    msc = zeros(Float32, M, 5)
    for k in 1:M, b in 1:4
        msc[k, b] = Float32(log(4) - profile.match_emissions[b, k])
    end
    into(r) = Float32[k == 1 ? -Inf32 : _lp(t[r, k]) for k in 1:M]
    out(r) = Float32[_lp(t[r, k+1]) for k in 1:M]
    occ = _occupancy(t)
    Z = sum(occ[k] * (M - k + 1) for k in 1:M)
    tbm = Float32.(log.(occ ./ Z))
    maxl = profile.maxl > 0 ? profile.maxl : 2M
    Plan7HMM(profile.name, M, maxl, msc, into(1), into(4), into(6),
             out(2), out(5), out(3), out(7), tbm,
             profile.stats_msv, profile.stats_viterbi, profile.stats_forward)
end
