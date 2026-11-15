# SSV: ungapped local diagonals in 8-bit SIMD, striped as HMMER's SSV filter.
#
# For each residue r, the lane for state k holds the best ungapped diagonal
# score ending at (r, k), floored at zero. A diagonal that reaches the
# threshold yields a seed (r, k) and restarts from zero, so a long match gives
# a seed every few dozen residues. Scores are scaled so that the threshold
# plus one emission fits in Int8. In the striped layout state k-1 sits in the
# previous vector of the same lane, so the update is one aligned add and max
# per vector, done in place; only the first vector of a row needs a shift.

const _W = 64          # 64 Int8 lanes: one AVX-512 register, two AVX2 registers
const _V8 = Vec{_W,Int8}

const _SHR8 = Val(ntuple(i -> i + _W - 2, _W))
@inline _shr8(v::_V8) = shufflevector(_V8(0), v, _SHR8)

struct _SSVProfile
    M::Int
    Q::Int
    cap::Int8             # threshold in scaled units
    emit::Vector{Int8}    # 5 striped blocks of Q*W lanes; block 5 = ambiguous base
end

@inline _sstripe(k, Q) = ((k - 1) % Q) * _W + (k - 1) ÷ Q + 1

function _SSVProfile(hmm::Plan7HMM, threshold_bits::Real)
    threshold_bits > 0 || throw(ArgumentError("seed threshold must be positive"))
    upb = min(4.0, 100 / threshold_bits)                      # units per bit
    Q = cld(hmm.M, _W); n = Q * _W
    emit = fill(Int8(-100), 5n)                               # padding states stay at 0
    for k in 1:hmm.M
        j = _sstripe(k, Q)
        for b in 1:4
            emit[(b-1)*n+j] = Int8(clamp(round(Int, upb / log(2) * hmm.msc[k, b]), -100, 27))
        end
        emit[4n+j] = Int8(0)
    end
    _SSVProfile(hmm.M, Q, Int8(round(Int, upb * threshold_bits)), emit)
end

"Scaled SSV emission of state `k` for base column `col` (1-5), for tests."
_ssv_emit(p::_SSVProfile, col, k) = p.emit[(col-1)*p.Q*_W+_sstripe(k, p.Q)]

"""
    ssv_seeds(hmm, xs; seed_bits) -> Vector{Tuple{Int,Int}}

Positions `(residue, model state)` at which an ungapped diagonal of `hmm`
against base indices `xs` reaches `seed_bits`. Used to choose the windows that
are scored with [`viterbi_score`](@ref).
"""
ssv_seeds(hmm::Plan7HMM, xs; seed_bits=default_seed_bits(hmm)) =
    (_check_indices(xs); sort!(_ssv!(Tuple{Int,Int}[], _SSVProfile(hmm, seed_bits), xs, 0)))

# Values never exceed cap + 27 <= 127: a lane at or above the threshold is
# reset to zero in the same row, so no upper clamp is needed.
function _ssv!(seeds, p::_SSVProfile, xs, offset)
    Q = p.Q; n = Q * _W
    dp = zeros(Int8, n)
    emit = p.emit; cap = p.cap; vzero = _V8(0); vcap = _V8(cap)
    last = (Q - 1) * _W + 1
    @inbounds for (r, x) in enumerate(xs)
        eo = (x == 0 ? 4 : Int(x) - 1) * n
        mpv = _shr8(vload(_V8, dp, last))
        acc = vzero
        for q in 0:Q-1
            j = q * _W + 1
            sv = max(mpv + vload(_V8, emit, eo + j), vzero)
            mpv = vload(_V8, dp, j)
            vstore(sv, dp, j)
            acc = max(acc, sv)
        end
        if any(acc >= vcap)
            for j in 1:n
                if dp[j] >= cap
                    q, l = divrem(j - 1, _W)
                    k = l * Q + q + 1
                    k <= p.M && push!(seeds, (offset + r, k))
                    dp[j] = 0
                end
            end
        end
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
