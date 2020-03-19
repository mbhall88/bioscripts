[TOC]: #

# Table of Contents
- [Basecall with guppy-gpu on an LSF cluster](#basecall-with-guppy-gpu-on-an-lsf-cluster)
- [Demultiplex with guppy on an LSF cluster](#demultiplex-with-guppy-on-an-lsf-cluster)
- [Convert multi-fast5 files to single-fast5s](#convert-multi-fast5-files-to-single-fast5s)
- [Convert single-fast5 files to multi-fast5s](#convert-single-fast5-files-to-multi-fast5s)


## Basecall with guppy-gpu on an LSF cluster
[`guppy_gpu.sh`][1]

```
$ ./guppy_gpu.sh --help
A script to assist with submitting a guppy basecalling job on an LSF cluster. The script
submits the job to run in a singularity (v3) container and makes some assumptions about the
location of CUDA libraries and binaries. Update if necessary.

Usage: guppy_gpu.sh -i <fast5_dir> -o <outdir> -j <jobname>

     -h|--help                  Displays this help
     -v|--verbose               Displays verbose output
    -nc|--no-colour             Disables colour output
    -cr|--cron                  Run silently unless we encounter an error
     -i|--input                 Directory containing fast5 files [Required]
     -o|--outdir                Directory for guppy output [Default: <current directory>]
     -j|--jobname               Name for the LSF job [Default: basecall]
     -g|--gpus                  Number of GPUs to use [Default: 1]
     -m|--memory                Memory (GB) to allocate for the job [Default: 4]
     -c|--config                Guppy config to use for basecalling [Default: dna_r9.4.1_450bps_hac_prom.cfg]
    -gv|--guppy-version         Version of guppy to use. Check valid versions at
                                https://cloud.sylabs.io/library/mbhall88/default/guppy-gpu
                                [Default: 3.4.5]
        --host                  GPU host to submit the job to [Default: gpu74-v100-hosts]
```

## Demultiplex with guppy on an LSF cluster
[`guppy_demux.sh`][4]

```
$ ./guppy_demux.sh --help
A script to assist with submitting a guppy demultiplexing job on an LSF cluster. The script
submits the job to run in a singularity (v3) container.

Usage: guppy_demux.sh -i <fastq_dir> -o <outdir> -j <jobname>

     -h|--help                  Displays this help
     -v|--verbose               Displays verbose output
    -nc|--no-colour             Disables colour output
    -cr|--cron                  Run silently unless we encounter an error
     -i|--input                 Directory containing fastq files [required]
     -o|--outdir                Directory for guppy output [<current directory>]
     -j|--jobname               Name for the LSF job [demux]
     -t|--threads               Number of threads to use [1]
     -m|--memory                Memory (GB) to allocate for the job [4]
     -k|--barcode-kits          Space separated list of barcoding kit(s) or expansion
                                kit(s) to detect against. Must be in double quotes ["EXP-NBD103"]
     -z|--compress              Compress the output files
        --trim                  Trim the barcodes from the output sequences
    -nr|--no-recursive          Don't recursively search for fastq in input directory
    -gv|--guppy-version         Version of guppy to use. Check valid versions at
                                https://cloud.sylabs.io/library/mbhall88/default/guppy-cpu
                                [3.4.5]
```

## Convert multi-fast5 files to single-fast5s
Because sometimes life just sucks!

Runs the job in a container (downloaded from biocontainers) using singularity.
The script will try and load singularity from a module with `module load singularity/3.5.0`. If
this is not relevant to your environment, then comment that line out.

[`multi_to_single_fast5.sh`][2]

```
$ ./multi_to_single_fast5.sh
Usage: ./multi_to_single_fast5.sh <indir> <outdir> [<threads>]
```

## Convert single-fast5 files to multi-fast5s
Runs the job in a container (downloaded from biocontainers) using singularity.
The script will try and load singularity from a module with `module load singularity/3.5.0`. If
this is not relevant to your environment, then comment that line out.

[`single_to_multi_fast5.sh`][3]

```
$ ./single_to_multi_fast5.sh
Usage: ./single_to_multi_fast5.sh <indir> <outdir> [<threads> <batch_size> <filename_base>]
```

[1]: https://github.com/mbhall88/bioscripts/blob/master/shell/guppy_gpu.sh
[2]: https://github.com/mbhall88/bioscripts/blob/master/shell/multi_to_single_fast5.sh
[3]: https://github.com/mbhall88/bioscripts/blob/master/shell/single_to_multi_fast5.sh
[4]: https://github.com/mbhall88/bioscripts/blob/master/shell/guppy_demux.sh