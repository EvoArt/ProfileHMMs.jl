# ProfileHMMs.jl

Search DNA sequences with HMMER3 nucleotide profile HMMs, in native Julia.
No HMMER executables are needed, so it runs the same on Windows, macOS and
Linux. [FindRNA.jl](https://github.com/EvoArt/FindRNA.jl) builds an rRNA
finder on it.

```julia
using ProfileHMMs
hits = nhmmscan("models.hmm", "genome.fna")              # every profile in the file
hmm = build_plan7(first(parse_hmm_file("models.hmm")))
hits = scan_fasta(hmm, "genome.fna"; evalue_threshold=1e-5)
hits = scan_fasta(hmm, "genome.fna"; scoring=:forward)   # nhmmer's scores and E-values
v = viterbi_score(hmm, seq_to_indices(dna"ACGT..."))     # score and endpoints
```

Each `SearchHit` gives the sequence, profile, 1-based inclusive coordinates on
the input sequence and on the profile, strand, bit score, bias and E-value.
Both strands are searched.

## How it searches

The pipeline follows nhmmer, with HMMER's striped SIMD layout (Farrar 2007) in
every dynamic-programming kernel:

1. **SSV seeds.** One ungapped pass over each strand finds diagonals scoring
   above `seed_bits` (default `3·log2(M) − 2` for a model of length `M`), with
   8-bit scores, 64 per vector.
2. **Viterbi.** Windows around the seeds (the model length plus room for
   insertions) are aligned by single-hit local Viterbi on a HMMER3 local-mode
   profile: occupancy-weighted entry, insert emissions scoring zero, and N/C
   states configured for each window's length. A 1-byte-per-cell traceback
   gives the exact start and end in the sequence and the profile. The flanks of
   each alignment are searched again, so that nearby copies are found.
3. **Scores and E-values.** With `scoring=:viterbi` (the default), the
   alignment's Viterbi score, with P-values from the profile's
   `STATS LOCAL VITERBI` calibration. With `scoring=:forward`, the alignment
   is rescored exactly as nhmmer 3.4 does for DNA: Forward score of the
   alignment ±20 bases, composition-bias correction (the score again with the
   background re-estimated from the region's own composition), normalisation
   to a MAXL-length window, and `STATS LOCAL FORWARD` P-values. This costs
   about 10% more time.

E-values are P-values times the search space, `2 × residues / MAXL`, as nhmmer.
Searches use all Julia threads (`threaded=false` to disable). In interactive
sessions a progress bar is shown on `stderr` (`progress=false`, or
`ProfileHMMs.setprogress!(false)` for every call).

## How it compares with HMMER

On five bacterial genomes searched with the bacterial rRNA profiles, the hits
are the same as nhmmer 3.4's, with endpoints within 6 bases. With
`scoring=:forward`, bit scores and bias corrections are within about one bit
of nhmmer's and E-values within a factor of two. The remaining differences:
alignment ends come from Viterbi, not nhmmer's optimal-accuracy alignment, and
the envelope is the alignment ±20 bases rather than one found by posterior
decoding.

Only DNA/RNA profiles are supported; profile building and protein alphabets
are not implemented.

## Licence

MIT. Display and progress code is shared with DenoiseDNA and AssignTaxonomy.
