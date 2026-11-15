using Test, Random, BioSequences, ProfileHMMs

function toy_profile(consensus="ACGTTGCAACGT")
    M = length(consensus)
    emit = fill(-log(0.01/3), 4, M)
    for (k, c) in enumerate(consensus)
        emit[findfirst(==(c), "ACGT"), k] = -log(0.99)
    end
    trans = repeat(-log.([0.9, 0.05, 0.05, 0.8, 0.2, 0.8, 0.2]), 1, M+1)
    HMMProfile("toy", "", "synthetic test model", M, M+10, "DNA", 1, 1.0,
        (-10.0, 0.7), (-10.0, 0.7), (-10.0, 0.7), fill(log(4),4),
        emit, fill(log(4),4,M+1), trans)
end

# Independent exhaustive path enumeration on small inputs, no DP recurrence.
# Enumerates every local start, match/insert/delete path, and local end.
function exhaustive(h, xs)
    best = (Inf, 0, 0)
    L = length(xs)
    emit(k,i) = xs[i] == 0 ? h.null1 : h.vit_match[k,xs[i]]
    function visit(state, k, i, start, score)
        if state == :M
            total = score + (L-i)*(h.nj_move+h.null1) + h.nj_exit
            total < best[1] && (best = (total,start,i))
        end
        if i < L && k < h.M
            cost = state == :M ? h.tr_mm[k+1] : state == :I ? h.tr_im[k+1] : h.tr_dm[k+1]
            visit(:M,k+1,i+1,start,score+cost+emit(k+1,i+1))
        end
        if i < L && state != :D
            cost = state == :M ? h.tr_mi[k] : h.tr_ii[k]
            e = xs[i+1] == 0 ? h.null1 : h.vit_insert[k,xs[i+1]]
            visit(:I,k,i+1,start,score+cost+e)
        end
        if k < h.M && state != :I
            cost = state == :M ? h.tr_md[k] : h.tr_dd[k]
            visit(:D,k+1,i,start,score+cost)
        end
    end
    for i in 1:L, k in 1:h.M
        visit(:M,k,i,i,(i-1)*(h.nj_move+h.null1)+h.nj_exit+h.t_bm+emit(k,i))
    end
    best
end

@testset "Plan7 correctness" begin
    p = toy_profile(); h = build_plan7(p; target_length=20)
    @test seq_to_indices(dna"ACGTN") == [1,2,3,4,0]
    @test_throws ArgumentError build_plan7(p; target_length=0)
    @test_throws ArgumentError viterbi_score(h, [5])
    @test viterbi_score(h, UInt8[]).target_start == 0
    rng = MersenneTwister(41)
    small = build_plan7(toy_profile("ACG"); target_length=5)
    for _ in 1:100
        xs = rand(rng, UInt8(0):UInt8(4),5)
        ref = exhaustive(small,xs)
        got = viterbi_score(small,xs)
        @test got.score ≈ ref[1] atol=1e-5
        @test 1 <= got.target_start <= got.target_end <= length(xs)
    end
    for s in ("ACGTTGCAACGT", "ACGTTTGCAACGT", "ACGTGCAACGT", "ACGTTNCAACGT")
        xs = seq_to_indices(LongDNA{4}("NNN" * s * "NNN"))
        v = viterbi_score(h,xs)
        @test (v.target_start,v.target_end) == (4,3+length(s))
    end
    @test isfinite(msv_score(h, zeros(UInt8,10)))
    @test msv_evalue(h, 1000.0, 10) > 0.9
    @test msv_evalue(h, 0.0, 100; db_size=10) ≈ 10msv_evalue(h,0.0,100)
    # Independent probability calculation for the geometric null length model.
    L = 20; p_continue = L/(L+1)
    null = L*log(4) - log(p_continue^L * (1-p_continue))
    @test ProfileHMMs._null_score(h,L) ≈ null atol=1e-6
end

@testset "Parser rejects bad inputs" begin
    mktemp() do path, io
        write(io,"not an HMM"); close(io)
        @test_throws ArgumentError parse_hmm_file(path)
    end
end

@testset "Search, windows and progress" begin
    h = build_plan7(toy_profile(); target_length=20)
    seq = dna"NNNACGTTGCAACGTNNNNNNACGTTGCAACGTNNN"
    opts = (msv_threshold=Inf, evalue_threshold=1.0, min_length=10)
    hits = search_sequence(h, seq; opts...)
    @test count(x -> x.strand == '+', hits) == 2
    mktemp() do path, io
        println(io,">copies\n",seq); close(io)
        serial = scan_fasta(h,path; window=24,overlap=16,threaded=false,progress=false,opts...)
        parallel = scan_fasta(h,path; window=24,overlap=16,threaded=true,progress=false,opts...)
        sig(xs) = [(x.target_start,x.target_end,x.strand,x.bit_score) for x in xs]
        @test sig(serial) == sig(parallel)
        @test count(x -> x.strand == '+', serial) == 2
        @test_throws ArgumentError scan_fasta(h,path; window=10,overlap=10)
        buf = IOBuffer()
        scan_fasta(h,path; progress=true,io=buf,opts...)
        shown = String(take!(buf))
        @test occursin("Searching",shown)
        @test !occursin('\e',shown)
        @test count(==('\n'),shown) == 1
    end
    @test !occursin('\e',sprint(show,MIME"text/plain"(),h))
end
