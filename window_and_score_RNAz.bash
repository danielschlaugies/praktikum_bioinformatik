#!/usr/bin/env bash

usage() {
    cat << 'EOF'
  Usage:
window_and_score-RNAz.bash -i ALIGNMENT.sto 

The input must be a curated Stockholm alignment. For each input it:
  1. reports alignment statistics with esl-alistat;
  2. creates RNAz windows with rnazWindow.pl;
  3. scores those windows with RNAz;

Options:
  -i FILE   input Stockholm alignment
  -h        show this help
EOF
}

input=

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
if [[ ! -s "$input" ]]; then
    echo "Input does not exist or is empty: $input" >&2
    exit 1
fi
if ! grep -q '^# STOCKHOLM 1\.0' "$input"; then
    echo "Input is not a Stockholm 1.0 alignment: $input" >&2
    exit 1
fi

# remove file extension
name="$(basename -s .sto $input)"

# print alignment statistics 
esl-alistat "$input" > "$name.alistat.txt"
echo "input alignment statistics:"
cat "${name}.alistat.txt"

# convert to clustal
clustal=$name.aln
esl-reformat clustal "$input" > "$clustal"

# window the alignment
windows="$name.windows.aln"
rnazWindow.pl \
    --window=120 \
    --slide=40 \
    --no-reference \
    "$clustal" > "$windows"

# run RNAz
RNAz -b -o "$name.rnaz.txt" "$windows"
