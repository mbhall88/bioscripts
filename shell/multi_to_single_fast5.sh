#!/usr/bin/env bash
# https://github.com/nanoporetech/ont_fast5_api#multi_to_single_fast5

if [[ $# -lt 2 || $# -gt 3 ]]; then
	echo "Error: Illegal number of parameters: $#"
	echo "Usage: $0 <indir> <outdir> [<threads>]"
	exit 2
fi

set -eux

# make sure singularity v3 is being used
module load singularity/3.5.0

container="docker://quay.io/biocontainers/ont-fast5-api:3.0.1--py_0"

indir="$1"
if [ ! -d "$indir" ]; then
	echo "Error: fast5 input directory does not exist: $indir"
	exit 1
fi

outdir="$2"
if [ ! -d "$outdir" ]; then
	echo "Output directory $outdir does not exist. Creating..."
	mkdir "$outdir"
fi

threads="${3:-1}"

singularity exec "$container" multi_to_single_fast5 --recursive \
	--input_path "$indir" \
	--save_path "$outdir" \
	--threads "$threads"
