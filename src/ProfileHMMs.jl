"""
    ProfileHMMs

Search nucleotide sequences with HMMER3 profile HMMs in native Julia: an
8-bit SIMD ungapped seed filter, then local Viterbi alignment of the seeded
windows, after nhmmer. No external executables; runs on Windows, macOS and Linux.
"""
module ProfileHMMs

using BioSequences, FASTX
using SIMD: Vec, vload, vstore
using Printf: @sprintf
using PrecompileTools: @setup_workload, @compile_workload
include("style.jl")
include("progress.jl")
include("HMMFile.jl")
include("ProfileHMM.jl")
include("Algorithms.jl")
include("Search.jl")
include("show.jl")
include("precompile.jl")

export HMMProfile, parse_hmm_file, Plan7HMM, build_plan7, seq_to_indices, base_index
export ViterbiResult, viterbi_score, ssv_seeds, default_seed_bits, viterbi_pvalue, search_space
export SearchHit, search_sequence, scan_fasta, nhmmscan, setprogress!
end
