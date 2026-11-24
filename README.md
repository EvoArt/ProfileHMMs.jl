# ProfileHMMs.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/ProfileHMMs.jl/dev/)
[![CI](https://github.com/EvoArt/ProfileHMMs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/EvoArt/ProfileHMMs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/EvoArt/ProfileHMMs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EvoArt/ProfileHMMs.jl)

Search DNA sequences with HMMER3 nucleotide profile HMMs, in native Julia.
No HMMER executables are needed, so it runs the same on Windows, macOS and
Linux. [RibosomalRNAFinder.jl](https://github.com/EvoArt/RibosomalRNAFinder.jl) builds an rRNA
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

## Speed compared with nhmmer

Searching a whole bacterial genome (both strands) with three rRNA profiles
(16S, 23S, 5S; barrnap 0.9's `bac.hmm`), E-value ≤ 1e-5, on one Intel
Skylake-AVX512 node. Seconds, median of three runs; ProfileHMMs after a
warm-up call.

| genome (size) | threads | ProfileHMMs | ProfileHMMs `scoring=:forward` | nhmmer 3.4 | nhmmer 3.5-dev |
|---|---|---|---|---|---|
| CP148752.1 (6.8 Mb) | 1 | 1.46 | 1.71 | 10.4 | 10.5 |
| CP148752.1 (6.8 Mb) | 4 | 0.86 | 0.89 | 4.48 | 4.66 |
| NZ_JBBHKP010000001.1 (3.7 Mb) | 1 | 0.82 | 0.92 | 4.61 | 4.73 |
| NZ_JBBHKP010000001.1 (3.7 Mb) | 4 | 0.44 | 0.46 | 2.13 | 2.44 |

All five genomes, the exact commands and versions are in
[benchmark/results.md](benchmark/results.md). nhmmer 3.4 is the latest
release; 3.5-dev is HMMER's development branch. On every genome both tools
report the same hits (69 in all), with endpoints within 6 bases.

Why it is faster, and what that costs:

- **Wider vectors.** HMMER's x86 nhmmer uses 128-bit SSE (its experimental
  AVX branch's nhmmer crashes), while ProfileHMMs uses 512-bit vectors for the
  seed scan where available. The instructions are the same kind; Julia
  compiles through LLVM, as clang does.
- **A stricter first filter.** nhmmer's first stage passes on ~2–3% of the
  sequence to its later stages, to catch weak, divergent hits; ProfileHMMs'
  seed threshold passes on ~0.2% (16S, 23S) to 0.7% (5S). Most of the gap comes
  from this. On these genomes nothing is lost, but for highly divergent or
  fragmentary matches nhmmer may be more sensitive; lower `seed_bits` to move
  ProfileHMMs towards nhmmer's sensitivity, at some cost in speed.
- **Less work per hit.** nhmmer also runs posterior decoding and an
  optimal-accuracy alignment for every hit.

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
