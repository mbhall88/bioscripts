[TOC]: #

# Table of Contents
- [Basecall with guppy-gpu on an LSF cluster](#basecall-with-guppy-gpu-on-an-lsf-cluster)
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