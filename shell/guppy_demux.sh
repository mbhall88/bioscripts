#!/usr/bin/env bash

# template taken from https://github.com/ralish/bash-script-template/blob/stable/template.sh
set -o errexit  # Exit on most errors (see the manual)
set -o errtrace # Make sure any error trap is inherited
set -o nounset  # Disallow expansion of unset variables
set -o pipefail # Use last non-zero exit code in a pipeline
#set -o xtrace          # Trace the execution of the script (debug)

# DESC: Handler for unexpected errors
# ARGS: $1 (optional): Exit code (defaults to 1)
# OUTS: None
function script_trap_err() {
	local exit_code=1

	# Disable the error trap handler to prevent potential recursion
	trap - ERR

	# Consider any further errors non-fatal to ensure we run to completion
	set +o errexit
	set +o pipefail

	# Validate any provided exit code
	if [[ ${1-} =~ ^[0-9]+$ ]]; then
		exit_code="$1"
	fi

	# Output debug data if in Cron mode
	if [[ -n ${cron-} ]]; then
		# Restore original file output descriptors
		if [[ -n ${script_output-} ]]; then
			exec 1>&3 2>&4
		fi

		# Print basic debugging information
		printf '%b\n' "$ta_none"
		printf '***** Abnormal termination of script *****\n'
		printf 'Script Path:            %s\n' "$script_path"
		printf 'Script Parameters:      %s\n' "$script_params"
		printf 'Script Exit Code:       %s\n' "$exit_code"

		# Print the script log if we have it. It's possible we may not if we
		# failed before we even called cron_init(). This can happen if bad
		# parameters were passed to the script so we bailed out very early.
		if [[ -n ${script_output-} ]]; then
			printf 'Script Output:\n\n%s' "$(cat "$script_output")"
		else
			printf 'Script Output:          None (failed before log init)\n'
		fi
	fi

	# Exit with failure status
	exit "$exit_code"
}

# DESC: Handler for exiting the script
# ARGS: None
# OUTS: None
function script_trap_exit() {
	cd "$orig_cwd"

	# Remove Cron mode script log
	if [[ -n ${cron-} && -f ${script_output-} ]]; then
		rm "$script_output"
	fi

	# Remove script execution lock
	if [[ -d ${script_lock-} ]]; then
		rmdir "$script_lock"
	fi

	# Restore terminal colours
	printf '%b' "$ta_none"
}

# DESC: Exit script with the given message
# ARGS: $1 (required): Message to print on exit
#       $2 (optional): Exit code (defaults to 0)
# OUTS: None
# NOTE: The convention used in this script for exit codes is:
#       0: Normal exit
#       1: Abnormal exit due to external error
#       2: Abnormal exit due to script error
function script_exit() {
	if [[ $# -eq 1 ]]; then
		printf '%s\n' "$1"
		exit 0
	fi

	if [[ ${2-} =~ ^[0-9]+$ ]]; then
		printf '%b\n' "$1"
		# If we've been provided a non-zero exit code run the error trap
		if [[ $2 -ne 0 ]]; then
			script_trap_err "$2"
		else
			exit 0
		fi
	fi

	script_exit 'Missing required argument to script_exit()!' 2
}

# DESC: Generic script initialisation
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: $orig_cwd: The current working directory when the script was run
#       $script_path: The full path to the script
#       $script_dir: The directory path of the script
#       $script_name: The file name of the script
#       $script_params: The original parameters provided to the script
#       $ta_none: The ANSI control code to reset all text attributes
# NOTE: $script_path only contains the path that was used to call the script
#       and will not resolve any symlinks which may be present in the path.
#       You can use a tool like realpath to obtain the "true" path. The same
#       caveat applies to both the $script_dir and $script_name variables.
# shellcheck disable=SC2034
function script_init() {
	# Useful paths
	readonly orig_cwd="$PWD"
	readonly script_path="${BASH_SOURCE[0]}"
	readonly script_dir="$(dirname "$script_path")"
	readonly script_name="$(basename "$script_path")"
	readonly script_params="$*"

	# Important to always set as we use it in the exit handler
	readonly ta_none="$(tput sgr0 2>/dev/null || true)"
}

# DESC: Initialise colour variables
# ARGS: None
# OUTS: Read-only variables with ANSI control codes
# NOTE: If --no-colour was set the variables will be empty
# shellcheck disable=SC2034
function colour_init() {
	if [[ -z ${no_colour-} ]]; then
		# Text attributes
		readonly ta_bold="$(tput bold 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly ta_uscore="$(tput smul 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly ta_blink="$(tput blink 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly ta_reverse="$(tput rev 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly ta_conceal="$(tput invis 2>/dev/null || true)"
		printf '%b' "$ta_none"

		# Foreground codes
		readonly fg_black="$(tput setaf 0 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_blue="$(tput setaf 4 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_cyan="$(tput setaf 6 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_green="$(tput setaf 2 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_magenta="$(tput setaf 5 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_red="$(tput setaf 1 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_white="$(tput setaf 7 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly fg_yellow="$(tput setaf 3 2>/dev/null || true)"
		printf '%b' "$ta_none"

		# Background codes
		readonly bg_black="$(tput setab 0 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_blue="$(tput setab 4 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_cyan="$(tput setab 6 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_green="$(tput setab 2 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_magenta="$(tput setab 5 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_red="$(tput setab 1 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_white="$(tput setab 7 2>/dev/null || true)"
		printf '%b' "$ta_none"
		readonly bg_yellow="$(tput setab 3 2>/dev/null || true)"
		printf '%b' "$ta_none"
	else
		# Text attributes
		readonly ta_bold=''
		readonly ta_uscore=''
		readonly ta_blink=''
		readonly ta_reverse=''
		readonly ta_conceal=''

		# Foreground codes
		readonly fg_black=''
		readonly fg_blue=''
		readonly fg_cyan=''
		readonly fg_green=''
		readonly fg_magenta=''
		readonly fg_red=''
		readonly fg_white=''
		readonly fg_yellow=''

		# Background codes
		readonly bg_black=''
		readonly bg_blue=''
		readonly bg_cyan=''
		readonly bg_green=''
		readonly bg_magenta=''
		readonly bg_red=''
		readonly bg_white=''
		readonly bg_yellow=''
	fi
}

# DESC: Initialise Cron mode
# ARGS: None
# OUTS: $script_output: Path to the file stdout & stderr was redirected to
function cron_init() {
	if [[ -n ${cron-} ]]; then
		# Redirect all output to a temporary file
		readonly script_output="$(mktemp --tmpdir "$script_name".XXXXX)"
		exec 3>&1 4>&2 1>"$script_output" 2>&1
	fi
}

# DESC: Acquire script lock
# ARGS: $1 (optional): Scope of script execution lock (system or user)
# OUTS: $script_lock: Path to the directory indicating we have the script lock
# NOTE: This lock implementation is extremely simple but should be reliable
#       across all platforms. It does *not* support locking a script with
#       symlinks or multiple hardlinks as there's no portable way of doing so.
#       If the lock was acquired it's automatically released on script exit.
function lock_init() {
	local lock_dir
	if [[ $1 = 'system' ]]; then
		lock_dir="/tmp/$script_name.lock"
	elif [[ $1 = 'user' ]]; then
		lock_dir="/tmp/$script_name.$UID.lock"
	else
		script_exit 'Missing or invalid argument to lock_init()!' 2
	fi

	if mkdir "$lock_dir" 2>/dev/null; then
		readonly script_lock="$lock_dir"
		verbose_print "Acquired script lock: $script_lock"
	else
		script_exit "Unable to acquire script lock: $lock_dir" 1
	fi
}

# DESC: Pretty print the provided string
# ARGS: $1 (required): Message to print (defaults to a green foreground)
#       $2 (optional): Colour to print the message with. This can be an ANSI
#                      escape code or one of the prepopulated colour variables.
#       $3 (optional): Set to any value to not append a new line to the message
# OUTS: None
function pretty_print() {
	if [[ $# -lt 1 ]]; then
		script_exit 'Missing required argument to pretty_print()!' 2
	fi

	if [[ -z ${no_colour-} ]]; then
		if [[ -n ${2-} ]]; then
			printf '%b' "$2"
		else
			printf '%b' "$fg_green"
		fi
	fi

	# Print message & reset text attributes
	if [[ -n ${3-} ]]; then
		printf '%s%b' "$1" "$ta_none"
	else
		printf '%s%b\n' "$1" "$ta_none"
	fi
}

# DESC: Only pretty_print() the provided string if verbose mode is enabled
# ARGS: $@ (required): Passed through to pretty_print() function
# OUTS: None
function verbose_print() {
	if [[ -n ${verbose-} ]]; then
		pretty_print "$@"
	fi
}

# DESC: Combines two path variables and removes any duplicates
# ARGS: $1 (required): Path(s) to join with the second argument
#       $2 (optional): Path(s) to join with the first argument
# OUTS: $build_path: The constructed path
# NOTE: Heavily inspired by: https://unix.stackexchange.com/a/40973
function build_path() {
	if [[ $# -lt 1 ]]; then
		script_exit 'Missing required argument to build_path()!' 2
	fi

	local new_path path_entry temp_path

	temp_path="$1:"
	if [[ -n ${2-} ]]; then
		temp_path="$temp_path$2:"
	fi

	new_path=
	while [[ -n $temp_path ]]; do
		path_entry="${temp_path%%:*}"
		case "$new_path:" in
		*:"$path_entry":*) ;;
		*)
			new_path="$new_path:$path_entry"
			;;
		esac
		temp_path="${temp_path#*:}"
	done

	# shellcheck disable=SC2034
	build_path="${new_path#:}"
}

# DESC: Check a binary exists in the search path
# ARGS: $1 (required): Name of the binary to test for existence
#       $2 (optional): Set to any value to treat failure as a fatal error
# OUTS: None
function check_binary() {
	if [[ $# -lt 1 ]]; then
		script_exit 'Missing required argument to check_binary()!' 2
	fi

	if ! command -v "$1" >/dev/null 2>&1; then
		if [[ -n ${2-} ]]; then
			script_exit "Missing dependency: Couldn't locate $1." 1
		else
			verbose_print "Missing dependency: $1" "${fg_red-}"
			return 1
		fi
	fi

	verbose_print "Found dependency: $1"
	return 0
}

# DESC: Validate we have superuser access as root (via sudo if requested)
# ARGS: $1 (optional): Set to any value to not attempt root access via sudo
# OUTS: None
function check_superuser() {
	local superuser
	if [[ $EUID -eq 0 ]]; then
		superuser=true
	elif [[ -z ${1-} ]]; then
		if check_binary sudo; then
			verbose_print 'Sudo: Updating cached credentials ...'
			if ! sudo -v; then
				verbose_print "Sudo: Couldn't acquire credentials ..." \
					"${fg_red-}"
			else
				local test_euid
				test_euid="$(sudo -H -- "$BASH" -c 'printf "%s" "$EUID"')"
				if [[ $test_euid -eq 0 ]]; then
					superuser=true
				fi
			fi
		fi
	fi

	if [[ -z ${superuser-} ]]; then
		verbose_print 'Unable to acquire superuser credentials.' "${fg_red-}"
		return 1
	fi

	verbose_print 'Successfully acquired superuser credentials.'
	return 0
}

# DESC: Run the requested command as root (via sudo if requested)
# ARGS: $1 (optional): Set to zero to not attempt execution via sudo
#       $@ (required): Passed through for execution as root user
# OUTS: None
function run_as_root() {
	if [[ $# -eq 0 ]]; then
		script_exit 'Missing required argument to run_as_root()!' 2
	fi

	if [[ ${1-} =~ ^0$ ]]; then
		local skip_sudo=true
		shift
	fi

	if [[ $EUID -eq 0 ]]; then
		"$@"
	elif [[ -z ${skip_sudo-} ]]; then
		sudo -H -- "$@"
	else
		script_exit "Unable to run requested command as root: $*" 1
	fi
}

# CLI defaults
guppy_version="3.4.5"
memory=4000
kits="EXP-NBD103"
compress=false
recurse=true
trim=false
threads=1
jobname="demux"
outdir=$(pwd)

# DESC: Usage help
# ARGS: None
# OUTS: None
function script_usage() {
	cat <<EOF
A script to assist with submitting a guppy demultiplexing job on an LSF cluster. The script
submits the job to run in a singularity (v3) container.

Usage: $(basename "$0") -i <fastq_dir> -o <outdir> -j <jobname>

     -h|--help                  Displays this help
     -v|--verbose               Displays verbose output
    -nc|--no-colour             Disables colour output
    -cr|--cron                  Run silently unless we encounter an error
     -i|--input                 Directory containing fastq files [required]
     -o|--outdir                Directory for guppy output [${outdir}]
     -j|--jobname               Name for the LSF job [${jobname}]
     -t|--threads               Number of threads to use [${threads}]
     -m|--memory                Memory (GB) to allocate for the job [$((memory / 1000))]
     -k|--barcode-kits          Space separated list of barcoding kit(s) or expansion
                                kit(s) to detect against. Must be in double quotes ["${kits}"]
     -z|--compress              Compress the output files
        --trim                  Trim the barcodes from the output sequences
    -nr|--no-recursive          Don't recursively search for fastq in input directory
    -gv|--guppy-version         Version of guppy to use. Check valid versions at
                                https://cloud.sylabs.io/library/mbhall88/default/guppy-cpu
                                [${guppy_version}]
EOF
}

# DESC: Check correct Singularity version is available or try and load module
# ARGS: None
# OUTS: None
function singularity_init() {
	regex="([0-9]+)\.[0-9]+\.[0-9]+"
	command -v singularity >/dev/null 2>&1 || script_exit "ERROR: Singularity is not available on your PATH" 1
	sing_version=$(singularity --version)
	if [[ ${sing_version} =~ $regex ]]; then
		major_version="${BASH_REMATCH[1]}"
		if ((major_version < 3)); then
			module load singularity/3.5.0 || script_exit "ERROR: Your singularity version is less than v3 and could not load v3 module" 1
		fi
	else
		script_exit "ERROR: Couldn't parse singularity major version from $sing_version" 1
	fi
}

# DESC: Parameter parser
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: Variables indicating command-line parameters and options
function parse_params() {
	local param
	while [[ $# -gt 0 ]]; do
		param="$1"
		shift  # past argument
		case ${param} in
		-h | --help)
			script_usage
			exit 0
			;;
		-v | --verbose)
			verbose=true
			;;
		-nc | --no-colour)
			no_colour=true
			;;
		-cr | --cron)
			cron=true
			;;
		-i | --input)
			input="$1"
			if [[ ! -d "$input" ]]; then
				script_exit "ERROR: Input directory $input does not exist" 1
			fi
            shift # past value
			;;
		-o | --outdir)
			outdir="$1"
			if [[ ! -d "$outdir" ]]; then
				pretty_print "WARNING: Output directory $outdir does not exist. Creating..."
				mkdir "$outdir"
			fi
            shift # past value
			;;
		-j | --jobname)
			jobname="$1"
            shift # past value
			;;
		-t | --threads)
			threads="$1"
            shift # past value
			;;
		-m | --memory)
            val="$1"
			memory=$((val * 1000))
            shift # past value
			;;
		-k | --barcode-kits)
			kits="$1"
            shift # past value
			;;
		-z | --compress)
			compress=true
			;;
		--trim)
		    trim=true
		    ;;
		-nr | --no-recursive)
		    recurse=false
		    ;;
		-gv | --guppy-version)
			guppy_version="$1"
            shift # past value
			;;
		*)
			script_exit "ERROR: Invalid parameter was provided: $param" 1
			;;
		esac
	done
	if [[ -z "${input+x}" ]]; then
		script_usage
		script_exit "ERROR: Input [-i|--input] is a required parameter and was not provided" 1
	fi
}

# DESC: Main control flow
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: None
function main() {
	trap script_trap_err ERR
	trap script_trap_exit EXIT

	script_init "$@"
	parse_params "$@"
	cron_init
	colour_init
	#lock_init system
	singularity_init

	container="library://mbhall88/default/guppy-cpu:${guppy_version}"

	flags=""
	if [[ "$recurse" = true ]]; then
	    flags="${flags} --recursive"
    fi
    if [[ "$compress" = true ]]; then
        flags="${flags} --compress_fastq"
    fi
    if [[ "$trim" = true ]]; then
        flags="${flags} --trim_barcodes"
    fi

	bsub -R "select[mem>${memory}] rusage[mem=${memory}]" \
    -M "$memory" \
    -o "$jobname".o \
    -e "$jobname".e \
    -J "$jobname" \
    singularity exec "$container" \
    guppy_barcoder \
        --input_path "$input" \
        --save_path "$outdir" \
        --barcode_kits "$kits" \
        --worker_threads "$threads" \
        "${flags}"
}

# Make it rain
main "$@"

# vim: syntax=sh cc=80 tw=79 ts=4 sw=4 sts=4 et sr