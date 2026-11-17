#!/usr/bin/env bash

usage() {
    cat << 'EOF'
  Usage:
search_homologes.bash -i QUERY.fa 

Runs a homologe search using RNAHub (and nhmmer implicitly).
The alignment of homologes and additional files created by RNAHub are stored in the current directory.

Options:
  -i FILE   query fasta file
  -h        show this help
EOF
}

while getopts ":i:h" opt; do
    case "$opt" in
        i) input=$OPTARG ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Missing argument for -$OPTARG" >&2
            usage >&2
            exit 2
            ;;
        \?)
            echo "Unknown option: -$OPTARG" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$input" ]]; then
    usage >&2
    exit 2
fi
if [[ ! -s "$HOME/shared/Bioinf26/Data/Transcripts/sequences/$input" ]]; then
    echo "Input does not exist or is empty: $input" >&2
    exit 1
fi

rnahub --db "$HOME/shared/Bioinf26/Data/Database/genomes" --verbose --cpus 11 --job-folder . --dev-skip-rscape --dev-skip-infernal --input "$HOME/shared/Bioinf26/Data/Transcripts/sequences/$input" --minidb --minidb-jobs 11