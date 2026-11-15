#!/bin/bash
#SBATCH --job-name=profilehmms-bench
#SBATCH --partition=hmq
#SBATCH --exclude=athena-h03
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
# ProfileHMMs against nhmmer 3.4 (latest release) and nhmmer from HMMER's
# h3-avx development branch (AVX2/AVX-512 kernels), at 1 and 4 threads.
# Submit from a checkout: sbatch benchmark/athena.sh
# Expects $TOOLS/bin/nhmmer (3.4) and $TOOLS/h3avx/bin/nhmmer.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
TOOLS="${TOOLS:-$HOME/findrna-tools}"
export PATH="$HOME/julia-1.12.7/bin:$PATH"
export NHMMERS="3.4=$TOOLS/bin/nhmmer,h3-avx=$TOOLS/h3avx/bin/nhmmer"
export BENCH_REPEATS=3
data=benchmark/data; mkdir -p "$data"
# barrnap 0.9's bacterial rRNA profiles (Rfam 16S, 23S, 5S)
if [ ! -s "$data/bac.hmm" ]; then
    curl -fsSL https://raw.githubusercontent.com/tseemann/barrnap/0.9/db/bac.hmm -o "$data/bac.hmm"
fi
echo "63e4093feeb24c6fb96f768df91ca9360210ed523af6151baaab4f2095e07bc3  $data/bac.hmm" | sha256sum -c -
for acc in CP148752.1 CP148753.1 JBBHKO010000001.1 NZ_JBBHKP010000001.1 NZ_JBBHKQ010000001.1; do
    [ -s "$data/$acc.fasta" ] && continue
    curl -fsSL --retry 5 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" -o "$data/$acc.fasta"
    sleep 1
done
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
out="benchmark/results/athena-$SLURM_JOB_ID"
for t in 1 4; do
    julia -t "$t" --startup-file=no --project=. benchmark/run.jl "$out" "$data/bac.hmm" "$data"/*.fasta
done
julia --startup-file=no benchmark/summary.jl "$out" | tee "$out/summary.md"
