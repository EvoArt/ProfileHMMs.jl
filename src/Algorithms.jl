# Dynamic programming on a local Plan7 profile.
#
#   viterbi_score   single-hit local Viterbi, exact endpoints, O(M) memory
#   ssv_seeds       ungapped diagonal filter (HMMER's SSV), 8-bit SIMD
#
# Scores are log-odds (higher is better). Ambiguous bases (index 0) score zero.

"""
    ViterbiResult

Best local alignment of a profile to a sequence: raw log-odds `score` (nats),
`bit_score` (with the length-dependent null model of HMMER), and the 1-based
inclusive span `target_start:target_end` of the aligned residues (0:0 if none).
"""
struct ViterbiResult
    score::Float64
    bit_score::Float64
    target_start::Int
    target_end::Int
end

function _check_indices(xs)
    all(x -> x <= 0x04, xs) || throw(ArgumentError("base indices must be in 0:4"))
end

mutable struct _ViterbiBuffer
    m::Vector{Float32}; i::Vector{Float32}; d::Vector{Float32}
    m2::Vector{Float32}; i2::Vector{Float32}
    om::Vector{Int32}; oi::Vector{Int32}; od::Vector{Int32}
    om2::Vector{Int32}; oi2::Vector{Int32}
end
function _ViterbiBuffer(M)
    f() = Vector{Float32}(undef, M); o() = Vector{Int32}(undef, M)
    _ViterbiBuffer(f(), f(), f(), f(), f(), o(), o(), o(), o(), o())
end

# NaN-free maximum: compiles to a vector max instruction.
@inline _mx(a, b) = ifelse(a > b, a, b)

function _rowmax(v)
    acc = -Inf32
    @inbounds @simd for k in eachindex(v)
        acc = _mx(acc, v[k])
    end
    acc
end

"""
    viterbi_score(hmm, xs) -> ViterbiResult

Score the best single local alignment of `hmm` to base indices `xs`. The
flanking N/C states and the null model are configured for `length(xs)` as in
HMMER, so scores do not drift with the length of the searched window.
"""
viterbi_score(hmm::Plan7HMM, xs) = (_check_indices(xs); _viterbi(hmm, xs, _ViterbiBuffer(hmm.M)))

# Score only: the filter applied to every seeded window.
_viterbi_bits(hmm::Plan7HMM, xs, buf::_ViterbiBuffer) = _viterbi(hmm, xs, buf, Val(false)).bit_score

# With `origins`, each state also carries the residue at which its alignment
# entered the model, giving exact endpoints without a traceback matrix.
function _viterbi(hmm::Plan7HMM, xs, buf::_ViterbiBuffer, ::Val{origins}=Val(true)) where {origins}
    L = length(xs); M = hmm.M
    L == 0 && return ViterbiResult(-Inf, -Inf, 0, 0)
    loop = log(L / (L + 2)); move = Float32(log(2 / (L + 2)))
    nullsc = L * log(L / (L + 1)) + log(1 / (L + 1))
    m, ins, d, m2, i2 = buf.m, buf.i, buf.d, buf.m2, buf.i2
    om, oi, od, om2, oi2 = buf.om, buf.oi, buf.od, buf.om2, buf.oi2
    fill!(m, -Inf32); fill!(ins, -Inf32); fill!(d, -Inf32)
    msc, tmm, tim, tdm, tmi, tii = hmm.msc, hmm.tmm, hmm.tim, hmm.tdm, hmm.tmi, hmm.tii
    tmd, tdd, tbm = hmm.tmd, hmm.tdd, hmm.tbm
    # C state: best exit so far plus C→C loops since, in Float64 because the
    # per-residue loop cost is below Float32 resolution at large scores.
    cbest = -Inf; rexit = 0; first_res = 0; last_res = 0
    @inbounds for r in 1:L
        x = xs[r]; col = x == 0 ? 5 : Int(x)
        enter = Float32((r - 1) * loop) + move              # N^(r-1) → B
        # Match and insert states read only the previous row: vectorised.
        m2[1] = enter + tbm[1] + msc[1, col]
        i2[1] = _mx(m[1] + tmi[1], ins[1] + tii[1])
        if origins
            om2[1] = r; oi2[1] = ifelse(m[1] + tmi[1] >= ins[1] + tii[1], om[1], oi[1])
        end
        @simd for k in 2:M
            s0 = enter + tbm[k]
            s1 = m[k-1] + tmm[k]; s2 = ins[k-1] + tim[k]; s3 = d[k-1] + tdm[k]
            best = _mx(_mx(s0, s1), _mx(s2, s3))
            m2[k] = best + msc[k, col]
            a = m[k] + tmi[k]; b = ins[k] + tii[k]
            i2[k] = _mx(a, b)
            if origins
                om2[k] = ifelse(best == s1, om[k-1], ifelse(best == s2, oi[k-1],
                         ifelse(best == s3, od[k-1], Int32(r))))
                oi2[k] = ifelse(a >= b, om[k], oi[k])
            end
        end
        # Deletions run along the current row: a serial chain.
        d[1] = -Inf32
        for k in 2:M
            a = m2[k-1] + tmd[k-1]; b = d[k-1] + tdd[k-1]
            d[k] = _mx(a, b)
            origins && (od[k] = ifelse(a >= b, om2[k-1], od[k-1]))
        end
        m, m2 = m2, m; ins, i2 = i2, ins
        origins && ((om, om2) = (om2, om); (oi, oi2) = (oi2, oi))
        # Local exit M_k → E → C. Exits through D cannot score higher.
        best = _rowmax(m)
        if best > cbest + (r - rexit) * loop
            cbest = Float64(best); rexit = r
            if origins
                last_res = r; first_res = om[findfirst(==(best), m)]
            end
        end
    end
    buf.m, buf.m2, buf.i, buf.i2 = m, m2, ins, i2
    buf.om, buf.om2, buf.oi, buf.oi2 = om, om2, oi, oi2
    raw = cbest + (L - rexit) * loop + move - nullsc
    ViterbiResult(raw, raw / log(2), first_res, last_res)
end

# ── SSV: ungapped local diagonals, 32 lanes of Int8 ─────────────────────────
#
# For each residue r, lane k holds the best ungapped diagonal score ending at
# (r, k), floored at zero. A diagonal that reaches the threshold yields a seed
# (r, k) and restarts from zero, so long matches give a seed every few dozen
# residues rather than one per residue. Scores are scaled so that the
# threshold plus one emission fits in Int8, as in HMMER's 8-bit filters.

const _W = 32
const _V8 = Vec{_W,Int8}

struct _SSVProfile
    M::Int
    nchunk::Int
    cap::Int8             # threshold in scaled units
    emit::Vector{Int8}    # 5 columns of nchunk*W lanes; column 5 = ambiguous base
end

function _SSVProfile(hmm::Plan7HMM, threshold_bits::Real)
    threshold_bits > 0 || throw(ArgumentError("seed threshold must be positive"))
    scale = min(4.0, 100 / threshold_bits) / log(2)          # units per nat
    nchunk = cld(hmm.M, _W)
    n = nchunk * _W
    emit = fill(Int8(-100), 5n)                              # padding lanes stay at 0
    for k in 1:hmm.M
        for b in 1:4
            emit[(b-1)*n+k] = Int8(clamp(round(Int, scale * hmm.msc[k, b]), -100, 27))
        end
        emit[4n+k] = Int8(0)
    end
    cap = Int8(round(Int, scale * log(2) * threshold_bits))
    _SSVProfile(hmm.M, nchunk, cap, emit)
end

"""
    ssv_seeds(hmm, xs; seed_bits) -> Vector{Tuple{Int,Int}}

Positions `(residue, model state)` at which an ungapped diagonal of `hmm`
against base indices `xs` reaches `seed_bits`. Used to choose the windows that
are scored with [`viterbi_score`](@ref).
"""
ssv_seeds(hmm::Plan7HMM, xs; seed_bits=default_seed_bits(hmm)) =
    (_check_indices(xs); _ssv!(Tuple{Int,Int}[], _SSVProfile(hmm, seed_bits), xs, 0))

function _ssv!(seeds, p::_SSVProfile, xs, offset)
    n = p.nchunk * _W
    prev = zeros(Int8, n + _W)     # index j holds state k = j-1; index 1 is the k=0 sentinel
    curr = zeros(Int8, n + _W)
    emit = p.emit; cap = p.cap
    vcap = _V8(cap); vzero = _V8(0)
    @inbounds for (r, x) in enumerate(xs)
        col = (x == 0 ? 4 : Int(x) - 1) * n
        acc = vzero
        for c in 0:p.nchunk-1
            j = c * _W + 1
            v = vload(_V8, prev, j) + vload(_V8, emit, col + j)
            v = min(max(v, vzero), vcap)
            vstore(v, curr, j + 1)
            acc = max(acc, v)
        end
        if maximum(acc) >= cap
            for k in 1:p.M
                if curr[k+1] >= cap
                    push!(seeds, (offset + r, k))
                    curr[k+1] = 0
                end
            end
        end
        prev, curr = curr, prev
    end
    seeds
end

"""
    default_seed_bits(hmm)

Ungapped seed threshold in bits. Chance seeds arise at a rate proportional to
the model length `M`, and each costs a Viterbi window of about `M × 1.5M`
cells, so the threshold grows as `3 log2(M)`: about 19 bits for 5S rRNA, 30
for 16S and 32 for 23S. Real hits have ungapped stretches far above these.
"""
default_seed_bits(hmm::Plan7HMM) = clamp(3log2(hmm.M) - 2, 14.0, 40.0)
