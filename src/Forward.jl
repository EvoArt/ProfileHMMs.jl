# nhmmer's rescoring of an alignment (HMMER 3.4, p7_domaindef.c
# rescore_isolated_domain and p7_pipeline.c, long-target case):
#
# 1. The envelope is the alignment plus up to 20 residues either side.
# 2. Forward score of the envelope, summing over all alignments, with N/C
#    states configured for the envelope's length.
# 3. Composition bias: the background is re-estimated as 25% uniform + 75% the
#    envelope's own base composition, match emission odds are recomputed
#    against it, and Forward is run again. The bias is the drop in score
#    (never negative), and the lower score is kept.
# 4. The score is moved to the scale of a MAXL-length window (entry/exit and
#    flanking costs of MAXL instead of the envelope) and the MAXL null model
#    is subtracted. P-values use the profile's Forward calibration.
#
# Forward is striped like the Viterbi kernel (HMMER's forward_engine), in
# probability (odds) space with Float32 lanes. The main loop carries only
# M → D; D → D paths are added by serial passes over the row, each crossing one
# more lane boundary, until they no longer change anything. Rows are rescaled
# when the exit probability exceeds 1e4 (HMMER's sparse rescaling).

struct _OddsProfile
    Q::Int
    em::Vector{Float32}      # 5 striped blocks of emission odds (block 5: ambiguous, 1)
    tbm::Vector{Float32}; tmm::Vector{Float32}; tim::Vector{Float32}; tdm::Vector{Float32}
    tmi::Vector{Float32}; tii::Vector{Float32}; tmd::Vector{Float32}; tdd::Vector{Float32}
end
_OddsProfile(h::Plan7HMM) = (vp = h.vp; _OddsProfile(vp.Q, exp.(vp.msc),
    exp.(vp.tbm), exp.(vp.tmm), exp.(vp.tim), exp.(vp.tdm),
    exp.(vp.tmi), exp.(vp.tii), exp.(vp.tmd), exp.(vp.tdd)))

"""
    ForwardResult

nhmmer-style score of an alignment: `envsc`, the composition-adjusted Forward
log-odds score of the envelope (nats); `bias`, the composition correction it
contains (nats; the unadjusted score minus `envsc`); and `bit_score`, the final
score in bits on nhmmer's scale.
"""
struct ForwardResult
    envsc::Float64
    bias::Float64
    bit_score::Float64
end

const _ENV_EXTRA = 20       # nhmmer's max_env_extra
const _BG_SMOOTH = 0.25     # nhmmer's bg_smooth for windows of >= 100 residues

"""
    forward_score(hmm, xs, a=1, b=length(xs); bias=true) -> ForwardResult

Rescore the alignment of `hmm` to `xs[a:b]` as nhmmer does: Forward score of
the envelope `a-20:b+20`, composition-bias correction (unless `bias=false`),
and normalisation to a window of the profile's MAXL.
"""
function forward_score(hmm::Plan7HMM, xs, a::Integer=1, b::Integer=length(xs); bias::Bool=true)
    _check_indices(xs)
    1 <= a <= b <= length(xs) || throw(ArgumentError("alignment span $a:$b is outside 1:$(length(xs))"))
    _nhmmer_score(hmm, _OddsProfile(hmm), xs, a, b, bias)
end

function _nhmmer_score(hmm::Plan7HMM, op::_OddsProfile, xs, a, b, bias::Bool)
    i = max(1, a - _ENV_EXTRA); j = min(length(xs), b + _ENV_EXTRA)
    env = view(xs, i:j); Le = j - i + 1; La = b - a + 1
    envsc = _forward(op, op.em, env)
    corr = 0.0
    if bias
        cnt = zeros(4)
        for x in env
            x == 0 || (cnt[x] += 1)
        end
        f = _BG_SMOOTH * 0.25 .+ (1 - _BG_SMOOTH) .* (sum(cnt) > 0 ? cnt ./ sum(cnt) : fill(0.25, 4))
        n = op.Q * _FW
        em = copy(op.em)
        for c in 1:4
            em[(c-1)*n+1:c*n] .*= Float32(0.25 / f[c])    # odds against f, not uniform
        end
        corr = max(0.0, envsc - _forward(op, em, env))
        envsc -= corr
    end
    X = hmm.maxl; Lx = max(X, Le)
    s = envsc - 2log(2 / (Le + 2)) + 2log(2 / (X + 2)) -
        (Le - La) * log(Le / (Le + 2)) + (Lx - La) * log(X / (X + 2))
    nullsc = Lx * log(Lx / (Lx + 1)) + log(1 / (Lx + 1))
    ForwardResult(envsc, corr, (s - nullsc) / log(2))
end

const _ZERO = _VF(0.0f0)
@inline _shr0(v::_VF) = shufflevector(_ZERO, v, _SHRIDX)

# Unihit local Forward (log-odds, nats) of `xs` with N/C configured for
# length(xs), using emission odds `em` (striped, as `op.em`). D → D chains
# decay into subnormal floats, which are very slow on x86, so subnormals are
# flushed to zero while the kernel runs (as HMMER does); the kernel does not
# yield, so the thread's previous setting is restored before any other task runs.
function _forward(op::_OddsProfile, em, xs)
    old = get_zero_subnormals()
    set_zero_subnormals(true)
    try
        return _forward_kernel(op, em, xs)
    finally
        set_zero_subnormals(old)
    end
end

function _forward_kernel(op::_OddsProfile, em, xs)
    L = length(xs)
    L == 0 && return -Inf
    Q = op.Q; n = Q * _FW; last = (Q - 1) * _FW + 1
    loop = L / (L + 2); move = 2 / (L + 2)
    tbm, tmm, tim, tdm = op.tbm, op.tmm, op.tim, op.tdm
    tmi, tii, tmd, tdd = op.tmi, op.tii, op.tmd, op.tdd
    m = zeros(Float32, n); ins = zeros(Float32, n); d = zeros(Float32, n)
    N = 1.0; C = 0.0; totscale = 0.0
    @inbounds for r in 1:L
        x = xs[r]; eo = (x == 0 ? 4 : Int(x) - 1) * n
        xB = _VF(Float32(N * move))
        mpv = _shr0(vload(_VF, m, last)); ipv = _shr0(vload(_VF, ins, last)); dpv = _shr0(vload(_VF, d, last))
        dcv = _ZERO; xE = _ZERO
        for q in 0:Q-1
            j = q * _FW + 1
            sv = (xB * vload(_VF, tbm, j) + mpv * vload(_VF, tmm, j) +
                  ipv * vload(_VF, tim, j) + dpv * vload(_VF, tdm, j)) * vload(_VF, em, eo + j)
            mpv = vload(_VF, m, j); ipv = vload(_VF, ins, j); dpv = vload(_VF, d, j)
            vstore(sv, m, j)
            vstore(dcv, d, j)
            vstore(mpv * vload(_VF, tmi, j) + ipv * vload(_VF, tii, j), ins, j)
            dcv = sv * vload(_VF, tmd, j)
            xE += sv
        end
        # D → D: one full pass (within lanes, plus the wrap), then passes that
        # carry the chain across further lane boundaries until it is negligible.
        dcv = _shr0(dcv)
        for q in 0:Q-1
            j = q * _FW + 1
            dq = vload(_VF, d, j) + dcv
            vstore(dq, d, j)
            dcv = dq * vload(_VF, tdd, j)
        end
        for _ in 2:_FW
            dcv = _shr0(dcv)
            changed = false
            for q in 0:Q-1
                j = q * _FW + 1
                dq = vload(_VF, d, j)
                nq = dq + dcv
                changed |= any(nq > dq)
                vstore(nq, d, j)
                dcv = dcv * vload(_VF, tdd, j)
            end
            changed || break
        end
        for q in 0:Q-1
            xE += vload(_VF, d, q * _FW + 1)           # local exits D_k → E
        end
        E = Float64(sum(xE))
        N *= loop; C = C * loop + E
        if E > 1e4
            s = 1 / E; sv = _VF(Float32(s))
            for q in 0:Q-1
                j = q * _FW + 1
                vstore(vload(_VF, m, j) * sv, m, j)
                vstore(vload(_VF, ins, j) * sv, ins, j)
                vstore(vload(_VF, d, j) * sv, d, j)
            end
            N *= s; C *= s; totscale += log(E)
        end
    end
    log(C * move) + totscale
end
