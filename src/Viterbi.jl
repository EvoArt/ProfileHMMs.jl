# Single-hit local Viterbi on the striped profile (Farrar 2007; HMMER's
# p7_ViterbiFilter), in Float32 with 8 lanes. Match and insert states of a row
# depend only on the previous row; delete states chain along the row, which
# the main loop handles within each lane and the "lazy F" passes complete
# across lanes, usually after a single vector.
#
# With `trace`, every cell also stores 1 byte of traceback pointers, from which
# the exact start of the alignment, in the sequence and the model, is recovered.

"""
    ViterbiResult

Best local alignment of a profile to a sequence: raw log-odds `score` (nats),
`bit_score` (with HMMER's length-dependent null model), the 1-based inclusive
aligned span `target_start:target_end` of the sequence and `model_start:model_end`
of the profile (all 0 if there is no alignment).
"""
struct ViterbiResult
    score::Float64
    bit_score::Float64
    target_start::Int
    target_end::Int
    model_start::Int
    model_end::Int
end

function _check_indices(xs)
    all(x -> x <= 0x04, xs) || throw(ArgumentError("base indices must be in 0:4"))
end

const _VF = Vec{_FW,Float32}
const _VB = Vec{_FW,UInt8}
const _VI = Vec{_FW,Int32}
const _NEG = _VF(-Inf32)

# Shift lanes up by one (lane l takes lane l-1), filling lane 0.
const _SHRIDX = Val(ntuple(i -> i + _FW - 2, _FW))
@inline _shr(v::_VF) = shufflevector(_NEG, v, _SHRIDX)
@inline _shr(v::_VB) = shufflevector(_VB(0), v, _SHRIDX)

# Pointer byte: bits 0-1 the source of M (0 B, 1 M, 2 I, 3 D), bit 2 the source
# of I (0 M, 1 I), bit 3 the source of D (0 M, 1 D).
const _P1, _P2, _P3, _PI, _PD = _VB(1), _VB(2), _VB(3), _VB(4), _VB(8)

mutable struct _ViterbiBuffer
    m::Vector{Float32}; i::Vector{Float32}; d::Vector{Float32}
    ptr::Vector{UInt8}
end
_ViterbiBuffer(hmm::Plan7HMM) = (n = hmm.vp.Q * _FW;
    _ViterbiBuffer(Vector{Float32}(undef, n), Vector{Float32}(undef, n), Vector{Float32}(undef, n), UInt8[]))

"""
    viterbi_score(hmm, xs) -> ViterbiResult

Score the best single local alignment of `hmm` to base indices `xs`. The
flanking N/C states and the null model are configured for `length(xs)` as in
HMMER, so scores do not drift with the length of the searched window.
"""
viterbi_score(hmm::Plan7HMM, xs) = (_check_indices(xs); _viterbi(hmm, xs, _ViterbiBuffer(hmm)))

# Score only: the filter applied to every seeded window.
_viterbi_bits(hmm::Plan7HMM, xs, buf::_ViterbiBuffer) = _viterbi(hmm, xs, buf, Val(false)).bit_score

function _viterbi(hmm::Plan7HMM, xs, buf::_ViterbiBuffer, ::Val{trace}=Val(true)) where {trace}
    L = length(xs)
    L == 0 && return ViterbiResult(-Inf, -Inf, 0, 0, 0, 0)
    vp = hmm.vp; Q = vp.Q; n = Q * _FW
    loop = log(L / (L + 2)); move = log(2 / (L + 2))
    nullsc = L * log(L / (L + 1)) + log(1 / (L + 1))
    m, ins, d = buf.m, buf.i, buf.d
    fill!(m, -Inf32); fill!(ins, -Inf32); fill!(d, -Inf32)
    ptr = buf.ptr
    trace && length(ptr) < L * n && resize!(ptr, L * n)
    msc, tbm, tmm, tim, tdm = vp.msc, vp.tbm, vp.tmm, vp.tim, vp.tdm
    tmi, tii, tmd, tdd = vp.tmi, vp.tii, vp.tmd, vp.tdd
    # C state: best exit so far plus C→C loops since, in Float64 because the
    # per-residue loop cost is below Float32 resolution at large scores.
    cbest = -Inf; rexit = 0; kexit = 0
    @inbounds for r in 1:L
        x = xs[r]; eo = (x == 0 ? 4 : Int(x) - 1) * n
        xB = _VF(Float32((r - 1) * loop + move))          # N^(r-1) → B
        row = (r - 1) * n
        last = (Q - 1) * _FW + 1
        mpv = _shr(vload(_VF, m, last)); ipv = _shr(vload(_VF, ins, last)); dpv = _shr(vload(_VF, d, last))
        dcv = _NEG; acc = _NEG; dptr = _VB(0); idx = _VI(0)
        for q in 0:Q-1
            j = q * _FW + 1
            s0 = xB + vload(_VF, tbm, j)
            s1 = mpv + vload(_VF, tmm, j)
            s2 = ipv + vload(_VF, tim, j)
            s3 = dpv + vload(_VF, tdm, j)
            best = max(max(s0, s1), max(s2, s3))
            sv = best + vload(_VF, msc, eo + j)
            mpv = vload(_VF, m, j); ipv = vload(_VF, ins, j); dpv = vload(_VF, d, j)
            vstore(sv, m, j)
            vstore(dcv, d, j)
            a = mpv + vload(_VF, tmi, j); b = ipv + vload(_VF, tii, j)
            vstore(max(a, b), ins, j)
            dm = sv + vload(_VF, tmd, j); dd = dcv + vload(_VF, tdd, j)
            if trace
                pm = vifelse(best == s1, _P1, vifelse(best == s2, _P2, vifelse(best == s3, _P3, _VB(0))))
                vstore(pm | vifelse(a >= b, _VB(0), _PI) | dptr, ptr, row + j)
                dptr = vifelse(dm >= dd, _VB(0), _PD)
                idx = vifelse(sv > acc, _VI(q), idx)
            end
            acc = max(acc, sv)
            dcv = max(dm, dd)
        end
        # Lazy F: carry D→D paths across lanes until no lane improves.
        dcv = _shr(dcv); trace && (dptr = _shr(dptr))
        while true
            complete = true
            for q in 0:Q-1
                j = q * _FW + 1
                dold = vload(_VF, d, j)
                up = dcv > dold
                if !any(up)
                    complete = false
                    break
                end
                dnew = max(dcv, dold)
                vstore(dnew, d, j)
                if trace
                    p = vload(_VB, ptr, row + j)
                    vstore(vifelse(up, (p & _VB(7)) | dptr, p), ptr, row + j)
                    dptr = _PD
                end
                dcv = dnew + vload(_VF, tdd, j)
            end
            complete || break
            dcv = _shr(dcv); trace && (dptr = _shr(dptr))
        end
        # Local exit M_k → E → C. Exits through D cannot score higher.
        best = maximum(acc)
        if best > cbest + (r - rexit) * loop
            cbest = Float64(best); rexit = r
            if trace
                l = 1                       # a loop, not a closure: closures box `acc`
                while acc[l] != best
                    l += 1
                end
                kexit = (l - 1) * Q + Int(idx[l]) + 1
            end
        end
    end
    raw = cbest + (L - rexit) * loop + move - nullsc
    trace || return ViterbiResult(raw, raw / log(2), 0, rexit, 0, 0)
    rexit == 0 && return ViterbiResult(raw, raw / log(2), 0, 0, 0, 0)
    rstart, kstart = _traceback(ptr, Q, rexit, kexit)
    ViterbiResult(raw, raw / log(2), rstart, rexit, kstart, kexit)
end

# Follow pointers from M_kexit at row rexit back to the B → M entry.
function _traceback(ptr, Q, r, k)
    n = Q * _FW
    state = 1                                   # 1 M, 2 I, 3 D
    @inbounds while true
        p = ptr[(r-1)*n+_stripe(k, Q)]
        if state == 1
            src = p & 0x03
            src == 0 && return (r, k)
            r -= 1; k -= 1; state = Int(src)
        elseif state == 2
            r -= 1; state = (p & 0x04) == 0 ? 1 : 2
        else
            k -= 1; state = (p & 0x08) == 0 ? 1 : 3
        end
    end
end
