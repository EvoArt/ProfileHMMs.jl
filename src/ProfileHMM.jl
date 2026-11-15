"""
    ProfileHMM.jl

Profile HMM data structures and construction from parsed HMM files.
Implements the HMMER3 Plan7 architecture with SIMD-optimized layout:
  - Float32 throughout (8-wide AVX2 SIMD)
  - M×4 emission layout for contiguous column access
  - Pre-extracted transition vectors for cache-friendly sequential access
  - Shared emission storage for the ungapped filter and Viterbi
"""

export Plan7HMM, build_plan7, seq_to_indices, base_index

const DNA_BASES = 4

# ──────────────────────────────────────────────────────────────────────
#  Optimized Plan7 HMM structure
# ──────────────────────────────────────────────────────────────────────

"""
    Plan7HMM

SIMD-optimized profile HMM using Float32 with M×4 emission layout.
MSV and Viterbi use negative log probabilities (lower is better).
"""
struct Plan7HMM
    name::String
    M::Int

    # ── Viterbi: neglog emissions + pre-extracted transitions (Float32) ──
    vit_match::Matrix{Float32}    # M × 4: match emission costs
    vit_insert::Matrix{Float32}   # M × 4: insert emission costs
    # Transitions pre-extracted per-node for sequential access:
    tr_mm::Vector{Float32}        # M: m→m from prev node (into node k)
    tr_im::Vector{Float32}        # M: i→m from prev node
    tr_dm::Vector{Float32}        # M: d→m from prev node
    tr_mi::Vector{Float32}        # M: m→i at current node
    tr_ii::Vector{Float32}        # M: i→i at current node
    tr_md::Vector{Float32}        # M: m→d from current node
    tr_dd::Vector{Float32}        # M: d→d from current node

    # ── Special state parameters (Float32 for speed, Float64 for E-value) ──
    nj_move::Float32              # N/C/J self-loop
    nj_exit::Float32              # N/C → next state
    null1::Float32                # background emission per residue (-log(0.25))
    t_bm::Float32                 # B→M entry cost

    # ── Float64 for precise E-value calculation ──
    null1_f::Float64
    stats_msv::Tuple{Float64, Float64}
    stats_viterbi::Tuple{Float64, Float64}
    stats_forward::Tuple{Float64, Float64}
end

# ──────────────────────────────────────────────────────────────────────
#  Utility functions
# ──────────────────────────────────────────────────────────────────────

"""Map DNA nucleotide to index 1-4. Returns 0 for ambiguous."""
@inline function base_index(nt)
    nt == DNA_A && return 1
    nt == DNA_C && return 2
    nt == DNA_G && return 3
    nt == DNA_T && return 4
    return 0
end

"""Convert a BioSequences DNA sequence to a vector of UInt8 indices (1-4, 0=ambiguous)."""
function seq_to_indices(seq)
    n = length(seq)
    indices = Vector{UInt8}(undef, n)
    @inbounds for i in 1:n
        indices[i] = base_index(seq[i])
    end
    indices
end

# ──────────────────────────────────────────────────────────────────────
#  Model construction
# ──────────────────────────────────────────────────────────────────────

"""
    build_plan7(profile::HMMProfile; target_length=400) -> Plan7HMM

Construct an optimized Plan7HMM with Float32 SIMD-friendly layout.
"""
function build_plan7(profile::HMMProfile; target_length=400)
    target_length > 0 || throw(ArgumentError("target_length must be positive"))
    profile.length > 0 || throw(ArgumentError("profile length must be positive"))
    M = profile.length
    K = DNA_BASES

    null1 = Float32(-log(0.25))   # 1.3863f0
    t_bm = Float32(-log(2.0 / (M * (M + 1))))

    L = target_length
    nj_move = Float32(-log(L / (L + 1)))
    nj_exit = Float32(-log(1 / (L + 1)))

    # ── Viterbi match/insert emissions (neglog, Float32) ──
    vit_match = Matrix{Float32}(undef, M, K)
    vit_insert = Matrix{Float32}(undef, M, K)
    for b in 1:K
        @inbounds for k in 1:M
            vit_match[k, b] = Float32(profile.match_emissions[b, k])
            icol = min(k + 1, size(profile.insert_emissions, 2))
            vit_insert[k, b] = Float32(profile.insert_emissions[b, icol])
        end
    end

    # ── Viterbi transitions (pre-extract per-node, Float32) ──
    tr = profile.transitions  # 7 × (M+1), neglog nats
    tr_mm = Vector{Float32}(undef, M)
    tr_im = Vector{Float32}(undef, M)
    tr_dm = Vector{Float32}(undef, M)
    tr_mi = Vector{Float32}(undef, M)
    tr_ii = Vector{Float32}(undef, M)
    tr_md = Vector{Float32}(undef, M)
    tr_dd = Vector{Float32}(undef, M)
    for k in 1:M
        tcol = k + 1
        # Transitions INTO node k FROM node k-1
        tr_mm[k] = k > 1 ? Float32(tr[1, tcol-1]) : 0.0f0
        tr_im[k] = k > 1 ? Float32(tr[4, tcol-1]) : 0.0f0
        tr_dm[k] = k > 1 ? Float32(tr[6, tcol-1]) : 0.0f0
        # Transitions FROM node k
        tr_mi[k] = Float32(tr[2, tcol])
        tr_ii[k] = Float32(tr[5, tcol])
        tr_md[k] = Float32(tr[3, tcol])
        tr_dd[k] = Float32(tr[7, tcol])
    end

    Plan7HMM(
        profile.name, M,
        vit_match, vit_insert,
        tr_mm, tr_im, tr_dm, tr_mi, tr_ii, tr_md, tr_dd,
        nj_move, nj_exit, null1, t_bm,
        Float64(null1),
        profile.stats_msv, profile.stats_viterbi, profile.stats_forward
    )
end
