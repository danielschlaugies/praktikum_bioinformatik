#!/usr/bin/env python3

import argparse
import re
import os
import subprocess
from pathlib import Path

PROBABILITY_RE = re.compile(r"SVM RNA-class probability:\s*([0-9.eE+-]+)")
DIRECTION_RE = re.compile(r"Reading direction:\s*(forward|reverse)")
COORDINATE_RE = re.compile(r"(?:_|/)(\d+)-(\d+)(?:\s|$)")

def read_window_coordinates(path: Path) -> list[tuple[int, int]]:
    coordinates: list[tuple[int, int]] = []
    in_block = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("CLUSTAL"):
                in_block = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            if not in_block:
                continue
            name = line.split()[0]
            match = COORDINATE_RE.search(name)
            if match:
                coordinates.append((int(match.group(1)), int(match.group(2))))
                in_block = False
    return coordinates

def read_rnaz_results(path: Path) -> list[tuple[str, float]]:
    results: list[tuple[str, float]] = []
    direction = "forward"
    with path.open() as handle:
        for line in handle:
            direction_match = DIRECTION_RE.search(line)
            if direction_match:
                direction = direction_match.group(1)
            match = PROBABILITY_RE.search(line)
            if match:
                results.append((direction, float(match.group(1))))
    return results

def preprocess(path: Path, start: int, end: int) -> str:
    # select columns in the alignment
    specified_columns = subprocess.run(
        [
            "esl-alimask",
            "-t",
            str(path),
            f"{start+1}..{end}"
        ],
        check=True,
        stdout=subprocess.PIPE,
        text=True
    )

    # remove all gap sequences
    empty_removed = subprocess.run(
        [
            "esl-alimanip",
            "--lmin",
            "1",
            "-"
        ],
        check=True,
        input=specified_columns.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )

    # now filter out columns that contain more than 50% gaps
    gap_column_filtered = subprocess.run(
        [
            "esl-alimask",
            # this defaults to 50% gaps
            "-g",
            "-"
        ],
        check=True,
        input=empty_removed.stdout,
        stdout=subprocess.PIPE,
        text=True
    )

    # TODO: might need to filter for gap only sequences at this point again

    return gap_column_filtered.stdout



def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alignment", type=Path, required=True)
    parser.add_argument("-w", "--windows", type=Path, required=True)
    parser.add_argument("-r", "--rnaz-report", type=Path, required=True)
    parser.add_argument("-t", "--threshold", type=float, default=0.5)
    args = parser.parse_args()

    coordinates = read_window_coordinates(args.windows)
    scores = read_rnaz_results(args.rnaz_report)

    for (i, (direction, value)) in enumerate(scores):
        # RNAz processes each window first in forward direction then in
        # reverse direction
        # hence we can get the start, respective endpoint,
        # from the position in in the coordinates list
        score_index = i // 2
        start, end = coordinates[score_index]

        if value > args.threshold:

            # skip reverse orientation if this has a lower value than forward
            if direction == "reverse" and value < scores[i-1][1]:
                continue

            base_name = f"{str(args.alignment)[:-4]}_{start + 1}-{end}"
            newdir = os.path.join(os.getcwd(), base_name)
            if not os.path.exists(newdir):
                os.mkdir(newdir)

            preprocessed_alignment = preprocess(args.alignment, start, end)
            preprocessed_filepath = os.path.join(newdir, f"{base_name}.sto")

            with open(preprocessed_filepath, 'w') as handle:
                handle.write(preprocessed_alignment)

            cacofold_out_dir = os.path.join(newdir, 'rscape_cacofold')
            if not os.path.exists(cacofold_out_dir):
                os.mkdir(cacofold_out_dir)

            # run R-scape and predict consensus structure with cacofold
            rcape = subprocess.run(
                [
                    "R-scape",
                    "--cacofold",
                    "--outdir",
                    cacofold_out_dir,
                    preprocessed_filepath
                ]
            )


    return 0

if __name__ == '__main__':
    raise SystemExit(main())

