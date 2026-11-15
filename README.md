# ProfileHMMs.jl

A native Julia engine for reading HMMER3 nucleotide profiles, ungapped filtering,
local Viterbi scoring and FASTA searches. Works without external executables on
Windows, Linux and macOS. FindRNA builds its RNA-specific workflow on this package.

```julia
using ProfileHMMs, BioSequences
profiles = parse_hmm_file("models.hmm")
hmm = build_plan7(first(profiles); target_length=1500)
score = viterbi_score(hmm, seq_to_indices(dna"ACGTACGT"))
hits = scan_fasta(hmm, "sequences.fna"; progress=true)
ProfileHMMs.setprogress!(false)
```

Install this subdirectory with `Pkg.develop(path="packages/ProfileHMMs")` while
it is unregistered. It is an independent Julia package: it does not import
FindRNA or read files from its parent directory. Run `Pkg.test("ProfileHMMs")`.

Supports DNA/RNA profiles only. Inputs use DNA symbols; index 0 denotes an
ambiguous base with neutral emission odds. Match/insert/delete recurrences
track alignment origins in O(model length) memory; `ViterbiResult.path` is
empty because full traceback is not retained. Search returns 1-based inclusive
coordinates on the original sequence and a `+` or `-` strand.

The algorithm is HMMER-inspired, not a complete HMMER port. MSV and Viterbi
statistics are approximate; Forward scoring, posterior domain decoding,
composition-bias correction, profile training and protein alphabets are not
implemented. These differences can affect significance, sensitivity and speed.

Displays and thread-safe progress bars follow DenoiseDNA/AssignTaxonomy. Logs
receive a single completion line; colours respect the IO context. CI disables
progress by default. Code is MIT-licensed; copied display helpers retain the
same author's original license.
