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
v = viterbi_score(hmm, seq_to_indices(dna"ACGT..."))     # score and endpoints
```

Each `SearchHit` gives the sequence, profile, 1-based inclusive coordinates on
the input sequence, strand, bit score and E-value. Both strands are searched.

## How it searches

The pipeline follows nhmmer:

1. **SSV seeds.** One ungapped pass over each strand finds diagonals scoring
   above `seed_bits` (default `3·log2(M) − 2` for a model of length `M`). The
   kernel uses 8-bit scores, 32 per SIMD instruction.
2. **Viterbi.** Windows around the seeds (the model length plus room for
   insertions) are scored by single-hit local Viterbi on a HMMER3 local-mode
   profile: occupancy-weighted entry, insert emissions scoring zero, and N/C
   states configured for each window's length. The flanks of each accepted
   alignment are searched again, so that nearby copies are found.
3. **E-values.** Viterbi P-values from the profile's `STATS LOCAL VITERBI`
   calibration are multiplied by the search space, `2 × residues / MAXL`.

Searches use all Julia threads (`threaded=false` to disable). In interactive
sessions a progress bar is shown on `stderr` (`progress=false`, or
`ProfileHMMs.setprogress!(false)` for every call).

## How it differs from HMMER

Scores are Viterbi scores, not nhmmer's Forward scores, and there is no
composition-bias (null2) correction. Bit scores therefore differ from nhmmer's
by a few percent, and E-values only approximate nhmmer's. Only DNA/RNA profiles
are supported; profile building and protein alphabets are not implemented.
On the five bacterial genomes FindRNA uses as a benchmark, the hits are the
same as nhmmer 3.4's, with endpoints within a few bases.

## Licence

MIT. Display and progress code is shared with DenoiseDNA and AssignTaxonomy.
