module ProfileHMMs

using BioSequences, FASTX
using Printf: @sprintf
include("style.jl")
include("progress.jl")
include("HMMFile.jl")
include("ProfileHMM.jl")
include("Algorithms.jl")
include("Search.jl")
include("show.jl")

export HMMProfile, parse_hmm_file, Plan7HMM, build_plan7, seq_to_indices, base_index
export msv_score, msv_evalue, viterbi_score, ViterbiResult
export SearchHit, search_sequence, scan_fasta, nhmmscan, setprogress!
end
