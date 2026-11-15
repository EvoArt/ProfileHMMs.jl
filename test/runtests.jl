using Test, Random, BioSequences, ProfileHMMs
const P = ProfileHMMs

function toy_profile(consensus="ACGTTGCAAGGT"; insert=fill(log(4), 4))
    M = length(consensus)
    emit = fill(-log(0.01 / 3), 4, M)
    for (k, c) in enumerate(consensus)
        emit[findfirst(==(c), "ACGT"), k] = -log(0.99)
    end
    trans = repeat(-log.([0.9, 0.05, 0.05, 0.8, 0.2, 0.8, 0.2]), 1, M + 1)
    HMMProfile("toy", "", "synthetic test model", M, M + 10, "DNA", 1, 1.0,
        (-10.0, 0.7), (-10.0, 0.7), (-10.0, 0.7), fill(log(4), 4),
        emit, repeat(insert, 1, M + 1), trans)
end

# Every local path (start residue/state, M/I/D moves, end) enumerated
# explicitly, scored from the definitions rather than the DP recurrences.
function exhaustive(h, xs)
    L = length(xs)
    loop = log(L / (L + 2)); move = log(2 / (L + 2))
    null = L * log(L / (L + 1)) + log(1 / (L + 1))
    e(k, i) = xs[i] == 0 ? 0.0 : Float64(h.msc[k, xs[i]])
    best = (-Inf, 0, 0)
    function visit(state, k, i, start, s)
        if state == :M
            total = s + (L - i) * loop + move - null
            total > best[1] && (best = (total, start, i))
        end
        if i < L && k < h.M
            t = state == :M ? h.tmm[k+1] : state == :I ? h.tim[k+1] : h.tdm[k+1]
            visit(:M, k + 1, i + 1, start, s + t + e(k + 1, i + 1))
        end
        if i < L && state != :D
            visit(:I, k, i + 1, start, s + (state == :M ? h.tmi[k] : h.tii[k]))
        end
        if k < h.M && state != :I
            visit(:D, k + 1, i, start, s + (state == :M ? h.tmd[k] : h.tdd[k]))
        end
    end
    for i in 1:L, k in 1:h.M
        visit(:M, k, i, i, (i - 1) * loop + move + h.tbm[k] + e(k, i))
    end
    best
end

# Scalar statement of the seed rule, for the SIMD kernel.
function ssv_reference(p::P._SSVProfile, xs)
    n = p.nchunk * P._W
    prev = zeros(Int, p.M + 1)
    seeds = Tuple{Int,Int}[]
    for (r, x) in enumerate(xs)
        col = (x == 0 ? 4 : x - 1) * n
        curr = zeros(Int, p.M + 1)
        for k in 1:p.M
            curr[k+1] = clamp(prev[k] + p.emit[col+k], 0, p.cap)
            if curr[k+1] >= p.cap
                push!(seeds, (r, k)); curr[k+1] = 0
            end
        end
        prev = curr
    end
    seeds
end

@testset "Profile configuration" begin
    h = build_plan7(toy_profile())
    @test h.M == 12 && h.maxl == 22
    @test seq_to_indices(dna"ACGTN") == [1, 2, 3, 4, 0]
    # Occupancy-weighted local entry is a distribution over (start, length).
    @test sum(exp(Float64(h.tbm[k])) * (h.M - k + 1) for k in 1:h.M) ≈ 1 atol = 1e-5
    @test all(iszero, h.msc[:, 5])
    @test h.msc[1, 1] ≈ log(4 * 0.99) atol = 1e-6
    # Insert emissions never score, whatever composition the file gives them.
    gc = build_plan7(toy_profile(; insert=-log.([0.1, 0.4, 0.4, 0.1])))
    seq = seq_to_indices(LongDNA{4}("ACGTTG" * "GC"^30 * "CAAGGT"))
    @test viterbi_score(gc, seq).score ≈ viterbi_score(h, seq).score
    @test_throws ArgumentError viterbi_score(h, UInt8[5])
end

@testset "Viterbi against exhaustive enumeration" begin
    rng = MersenneTwister(41)
    small = build_plan7(toy_profile("ACG"))
    for _ in 1:150
        xs = rand(rng, UInt8(0):UInt8(4), rand(rng, 1:5))
        ref = exhaustive(small, xs)
        got = viterbi_score(small, xs)
        @test got.score ≈ ref[1] atol = 1e-4
        @test got.bit_score ≈ ref[1] / log(2) atol = 1e-4
        @test P._viterbi_bits(small, xs, P._ViterbiBuffer(small.M)) ≈ got.bit_score
        @test 1 <= got.target_start <= got.target_end <= length(xs)
    end
    h = build_plan7(toy_profile())
    @test viterbi_score(h, UInt8[]).target_start == 0
    for s in ("ACGTTGCAAGGT", "ACGTTTGCAAGGT", "ACGTGCAAGGT", "ACGTTNCAAGGT")
        v = viterbi_score(h, seq_to_indices(LongDNA{4}("NNN" * s * "NNN")))
        @test (v.target_start, v.target_end) == (4, 3 + length(s))
    end
    # Length-configured N/C states, as HMMER: the score falls by about log2 of the
    # ratio of window lengths, not in proportion to the flanking length.
    core = viterbi_score(h, seq_to_indices(dna"ACGTTGCAAGGT")).bit_score
    long = viterbi_score(h, seq_to_indices(LongDNA{4}("N"^5000 * "ACGTTGCAAGGT" * "N"^5000))).bit_score
    @test 0 < core - long < log2(10012 / 12) + 3
end

@testset "SIMD seed scan" begin
    rng = MersenneTwister(7)
    for M in (5, 31, 32, 33, 70)
        cons = String(rand(rng, collect("ACGT"), M))
        h = build_plan7(toy_profile(cons))
        p = P._SSVProfile(h, 8.0)
        planted = string(cons[1:min(M, 20)])
        xs = seq_to_indices(LongDNA{4}(String(randdnaseq(rng, 300)) * "NN" * planted * String(randdnaseq(rng, 200))))
        got = P._ssv!(Tuple{Int,Int}[], p, xs, 0)
        @test got == ssv_reference(p, xs)
        @test !isempty(got)
        @test ssv_seeds(h, xs; seed_bits=8.0) == got
    end
    h = build_plan7(toy_profile())
    @test isempty(ssv_seeds(h, zeros(UInt8, 500)))
    @test_throws ArgumentError ssv_seeds(h, UInt8[1, 2]; seed_bits=0)
    @test 14 <= default_seed_bits(h) <= 40
end

@testset "Statistics" begin
    h = build_plan7(toy_profile())
    @test viterbi_pvalue(h, 1000.0) < 1e-100
    @test viterbi_pvalue(h, -100.0) ≈ 1
    @test search_space(h, 1100) ≈ 100
    @test search_space(h, 1) == 1
end

@testset "Search, threads and progress" begin
    h = build_plan7(toy_profile())
    seq = dna"NNNACGTTGCAAGGTNNNNNNNNNNNNNNNNNNNNNNNNNNACGTTGCAAGGTNNN"
    opts = (evalue_threshold=1.0, min_length=10)
    hits = search_sequence(h, seq; opts...)
    @test [(x.target_start, x.target_end, x.strand) for x in hits] == [(4, 15, '+'), (42, 53, '+')]
    rc = search_sequence(h, reverse_complement(seq); opts...)
    @test [(x.target_start, x.target_end, x.strand) for x in rc] ==
          [(length(seq) - 52, length(seq) - 41, '-'), (length(seq) - 14, length(seq) - 3, '-')]
    @test all(x -> x.seq_length == length(seq), hits)
    @test isempty(search_sequence(h, seq; evalue_threshold=1e-300))
    @test_throws ArgumentError search_sequence(h, seq; evalue_threshold=0)
    mktemp() do path, io
        println(io, ">copies\n", seq, "\n>empty\n\n>other\n", reverse_complement(seq)); close(io)
        serial = scan_fasta(h, path; threaded=false, progress=false, opts...)
        parallel = scan_fasta(h, path; threaded=true, progress=false, opts...)
        sig(xs) = [(x.target_name, x.target_start, x.target_end, x.strand, x.bit_score) for x in xs]
        @test sig(serial) == sig(parallel)
        @test length(serial) == 4
        @test all(x -> x.evalue ≈ search_space(h, 2length(seq)) * viterbi_pvalue(h, x.bit_score), serial)
        buf = IOBuffer()
        scan_fasta(h, path; progress=true, io=buf, opts...)
        shown = String(take!(buf))
        @test occursin("Searching", shown)
        @test !occursin('\e', shown)
        @test count(==('\n'), shown) == 1
        @test occursin("toy", shown)
    end
    mktemp() do path, io
        println(io, ">a\nACGT\n>a\nACGT"); close(io)
        @test_throws ArgumentError scan_fasta(h, path; progress=false)
    end
    @test !occursin('\e', sprint(show, MIME"text/plain"(), h))
end

@testset "Parser rejects bad inputs" begin
    mktemp() do path, io
        write(io, "not an HMM"); close(io)
        @test_throws ArgumentError parse_hmm_file(path)
    end
end
