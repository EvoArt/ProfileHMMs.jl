# Benchmark: ProfileHMMs and nhmmer

Athena (University of Exeter) compute node, job 4052420, 2026-11-15. Julia 1.12.7.
nhmmer 3.4 is the latest HMMER release; 3.5-dev is HMMER's `develop` branch
(commit 3221e22, banner "HMMER 3.5 (July 2026)"), the newest HMMER code. Both use
HMMER's SSE kernels: the `h3-avx` branch (AVX2/AVX-512) builds, but its nhmmer
segfaults on every input.

Input: barrnap 0.9's `bac.hmm` (Rfam 16S, 23S, 5S rRNA profiles) against five
bacterial genomes from NCBI (3.7–6.8 Mb), both strands, E-value ≤ 1e-5.
nhmmer: `--cpu 0` (serial) for 1 thread, `--cpu 4` for 4; `--noali --tblout`.
ProfileHMMs: `scan_fasta` with `julia -t 1` or `-t 4`, after one warm-up call.
Reproduce with `sbatch benchmark/athena.sh` (needs the two nhmmer builds).

CPU: `skylake-avx512`. Median of repeated runs, whole-genome search with all three profiles.

**1 thread** (seconds)

| genome | ProfileHMMs viterbi | ProfileHMMs forward | nhmmer 3.4 | nhmmer 3.5-dev |
|---|---|---|---|---|
| CP148752.1 | 1.46 | 1.71 | 10.4 | 10.5 |
| CP148753.1 | 1.29 | 1.42 | 8.02 | 8.13 |
| JBBHKO010000001.1 | 1.13 | 1.23 | 6.76 | 6.87 |
| NZ_JBBHKP010000001.1 | 0.82 | 0.92 | 4.61 | 4.73 |
| NZ_JBBHKQ010000001.1 | 0.86 | 0.97 | 5.3 | 5.42 |

**4 threads** (seconds)

| genome | ProfileHMMs viterbi | ProfileHMMs forward | nhmmer 3.4 | nhmmer 3.5-dev |
|---|---|---|---|---|
| CP148752.1 | 0.86 | 0.89 | 4.48 | 4.66 |
| CP148753.1 | 0.72 | 0.78 | 3.58 | 3.75 |
| JBBHKO010000001.1 | 0.61 | 0.64 | 2.87 | 3.16 |
| NZ_JBBHKP010000001.1 | 0.44 | 0.46 | 2.13 | 2.44 |
| NZ_JBBHKQ010000001.1 | 0.45 | 0.49 | 2.55 | 2.69 |

**Agreement** with ProfileHMMs' hits (matched / missing / extra, largest endpoint difference)

| genome | nhmmer 3.4 | nhmmer 3.5-dev |
|---|---|---|
| CP148752.1 | 25 / 0 / 0 (2 bp) | 25 / 0 / 0 (2 bp) |
| CP148753.1 | 13 / 0 / 0 (3 bp) | 13 / 0 / 0 (3 bp) |
| JBBHKO010000001.1 | 9 / 0 / 0 (3 bp) | 9 / 0 / 0 (3 bp) |
| NZ_JBBHKP010000001.1 | 10 / 0 / 0 (6 bp) | 10 / 0 / 0 (6 bp) |
| NZ_JBBHKQ010000001.1 | 12 / 0 / 0 (4 bp) | 12 / 0 / 0 (4 bp) |

**Work**: fraction of residues (both strands) passed on by the first filter, mean over genomes

| tool | 16S_rRNA | 23S_rRNA | 5S_rRNA |
|---|---|---|---|
| ProfileHMMs | 0.19% | 0.2% | 0.67% |
| nhmmer 3.4 | 2.61% | 2.37% | 1.63% |
| nhmmer 3.5-dev | 2.61% | 2.35% | 1.63% |
