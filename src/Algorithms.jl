# Nucleotide Plan7 dynamic programming. Scores are negative log probabilities.
# Carry alignment origins with each DP state, avoiding a full traceback matrix.
const INF32 = Float32(Inf)

struct ViterbiResult
    score::Float64
    bit_score::Float64
    path::Vector{Tuple{Int,Symbol}} # Reserved; no traceback is materialised.
    target_start::Int
    target_end::Int
end

function _check_indices(xs)
    all(x -> 0 <= x <= 4, xs) || throw(ArgumentError("indices must be in 0:4"))
end

"Ungapped MSV filter score; ambiguous bases emit with background probability."
function msv_score(hmm::Plan7HMM, xs)
    _check_indices(xs)
    M = hmm.M
    prev = fill(INF32, M)
    curr = similar(prev)
    n = 0.0f0
    c = INF32
    @inbounds for x in xs
        entry = n + hmm.nj_exit + hmm.t_bm
        curr[1] = entry + (x == 0 ? hmm.null1 : hmm.vit_match[1,x])
        @simd for k in 2:M
            emission = x == 0 ? hmm.null1 : hmm.vit_match[k,x]
            curr[k] = min(prev[k-1], entry) + emission
        end
        c = min(minimum(curr), c + hmm.nj_move + hmm.null1)
        n += hmm.nj_move + hmm.null1
        prev, curr = curr, prev
    end
    Float64(c + hmm.nj_exit)
end

# The null model includes the geometric length probability as well as base
# emissions. Omitting this term penalises short RNAs in long search windows.
function _null_score(hmm::Plan7HMM, len)
    len >= 0 || throw(ArgumentError("sequence length must be nonnegative"))
    length_cost = len == 0 ? 0.0 : log1p(len) + len * log1p(1.0 / len)
    len * hmm.null1_f + length_cost
end

function _tail_probability(bits, stats)
    mu, lambda = stats
    isfinite(mu) && isfinite(lambda) && lambda > 0 ||
        throw(ArgumentError("profile lacks valid calibration statistics"))
    -expm1(-exp(-lambda * (bits - mu)))
end

"Approximate MSV E-value using bit-score calibration and an explicit search count."
function msv_evalue(hmm::Plan7HMM, raw, len; db_size=1)
    db_size >= 1 || throw(ArgumentError("db_size must be at least one"))
    bits = (_null_score(hmm,len) - raw) / log(2.0)
    db_size * _tail_probability(bits, hmm.stats_msv)
end

"""
    viterbi_score(hmm, indices)

Single-hit local Viterbi score and exact endpoints of its optimal path.
Insertions consume a residue from the previous row; deletions consume a model
state in the current row. Ambiguous bases have neutral emission odds.
The path field is empty: only endpoints are tracked, in O(model length) memory.
"""
function viterbi_score(hmm::Plan7HMM, xs)
    _check_indices(xs)
    M = hmm.M
    pm = fill(INF32, M); pi = copy(pm); pd = copy(pm)
    cm = copy(pm); ci = copy(pm); cd = copy(pm)
    sm = zeros(Int, M); si = copy(sm); sd = copy(sm)
    tm = copy(sm); ti = copy(sm); td = copy(sm)
    n = 0.0f0; c = INF32
    firstpos = 0; lastpos = 0
    @inbounds for (i, x) in enumerate(xs)
        entry = n + hmm.nj_exit + hmm.t_bm
        for k in 1:M
            best = entry; start = i
            if k > 1
                a = pm[k-1] + hmm.tr_mm[k]
                b = pi[k-1] + hmm.tr_im[k]
                d = pd[k-1] + hmm.tr_dm[k]
                if a < best; best = a; start = sm[k-1]; end
                if b < best; best = b; start = si[k-1]; end
                if d < best; best = d; start = sd[k-1]; end
            end
            cm[k] = best + (x == 0 ? hmm.null1 : hmm.vit_match[k,x])
            tm[k] = start
            a = pm[k] + hmm.tr_mi[k]
            b = pi[k] + hmm.tr_ii[k]
            ci[k] = min(a,b) + (x == 0 ? hmm.null1 : hmm.vit_insert[k,x])
            ti[k] = a <= b ? sm[k] : si[k]
            if k == 1
                cd[k] = INF32; td[k] = 0
            else
                a = cm[k-1] + hmm.tr_md[k-1]
                b = cd[k-1] + hmm.tr_dd[k-1]
                cd[k] = min(a,b)
                td[k] = a <= b ? tm[k-1] : td[k-1]
            end
        end
        # Local exits from match states. A terminal deletion cannot improve
        # a local exit because transition costs are nonnegative.
        best, k = findmin(cm)
        old = c + hmm.nj_move + hmm.null1
        if best < old
            c = best; firstpos = tm[k]; lastpos = i
        else
            c = old
        end
        n += hmm.nj_move + hmm.null1
        pm,cm = cm,pm; pi,ci = ci,pi; pd,cd = cd,pd
        sm,tm = tm,sm; si,ti = ti,si; sd,td = td,sd
    end
    raw = Float64(c + hmm.nj_exit)
    bits = (_null_score(hmm,length(xs)) - raw) / log(2.0)
    ViterbiResult(raw, bits, Tuple{Int,Symbol}[], firstpos, lastpos)
end
