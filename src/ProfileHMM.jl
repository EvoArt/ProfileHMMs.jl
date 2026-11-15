# Search profile in HMMER3 local mode, as log-odds scores (nats, Float32,
# higher is better). As in HMMER, insert emissions score zero, so a profile's
# insert-state composition cannot reward biased target composition, and local
# entry into match state k is weighted by the occupancy of k.

# SIMD width (Float32 lanes) of the striped Viterbi kernel.
const _FW = 8

# Striped layout (Farrar 2007, as in HMMER): state k lives in vector
# q = (k-1) % Q, lane l = (k-1) ÷ Q. State k-1 is then in vector q-1 of the
# same lane, except for q = 0, which needs the last vector shifted by a lane.
struct _VitProfile
    Q::Int
    msc::Vector{Float32}     # 5 blocks of Q*_FW: emission by base (5 = ambiguous, 0)
    tbm::Vector{Float32}     # B → M_k
    tmm::Vector{Float32}     # M_{k-1} → M_k
    tim::Vector{Float32}     # I_{k-1} → M_k
    tdm::Vector{Float32}     # D_{k-1} → M_k
    tmi::Vector{Float32}     # M_k → I_k
    tii::Vector{Float32}     # I_k → I_k
    tmd::Vector{Float32}     # M_k → D_{k+1}
    tdd::Vector{Float32}     # D_k → D_{k+1}
end

@inline _stripe(k, Q) = ((k - 1) % Q) * _FW + (k - 1) ÷ Q + 1

function _vitprofile(M, msc, tbm, tmm, tim, tdm, tmi, tii, tmd, tdd)
    Q = cld(M, _FW); n = Q * _FW
    function s(v)        # padding states k > M are unreachable (-Inf transitions)
        a = fill(-Inf32, n)
        for k in 1:M
            a[_stripe(k, Q)] = v[k]
        end
        a
    end
    # Plan7 has no states beyond node M: in the striped layout a transition out
    # of node M would feed a padding lane, so remove it (the files may give one).
    tmd, tdd, tmi, tii = (setindex!(copy(v), -Inf32, M) for v in (tmd, tdd, tmi, tii))
    e = zeros(Float32, 5n)
    for b in 1:5, k in 1:M
        e[(b-1)*n+_stripe(k, Q)] = msc[k, b]
    end
    _VitProfile(Q, e, s(tbm), s(tmm), s(tim), s(tdm), s(tmi), s(tii), s(tmd), s(tdd))
end

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
    vp::_VitProfile             # the same scores, striped for the SIMD Viterbi
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
    tmm, tim, tdm = into(1), into(4), into(6)
    tmi, tii, tmd, tdd = out(2), out(5), out(3), out(7)
    Plan7HMM(profile.name, M, maxl, msc, tmm, tim, tdm, tmi, tii, tmd, tdd, tbm,
             profile.stats_msv, profile.stats_viterbi, profile.stats_forward,
             _vitprofile(M, msc, tbm, tmm, tim, tdm, tmi, tii, tmd, tdd))
end
