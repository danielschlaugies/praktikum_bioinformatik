#!/usr/bin/env python3

import argparse
import csv
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


WINDOW_SIZE = 120
WINDOW_SLIDE = 40
MAX_SEQUENCES = 6
OPT_IDENTITY = 60


@dataclass(frozen=True)
class RnazCandidate:
    transcript_id: str
    start: int
    end: int
    direction: str
    n_sequences: int
    columns: int
    average_identity: float
    mean_z_score: float
    sci: float
    svm_value: float
    rna_probability: float
    prediction: str

    @property
    def window(self) -> str:
        return f"{self.start}-{self.end}"


class ToolError(RuntimeError):
    pass


def tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise ToolError(
            f"Required tool '{name}' is not available. "
            "Run this script inside the praktikum_bioinformatik Conda environment."
        )
    return executable


def run_tool(
    args: list[str],
    *,
    input_text: str | None = None,
    cwd: Path | None = None,
) -> str:
    try:
        result = subprocess.run(
            args,
            input=input_text,
            text=True,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip() or "no diagnostics"
        raise ToolError(f"{Path(args[0]).name} failed: {detail}") from error
    return result.stdout


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def parse_rnaz_report(
    report_path: Path,
    transcript_id: str,
    allowed_windows: set[tuple[int, int]],
) -> list[RnazCandidate]:
    report = report_path.read_text()
    chunks = re.split(r"#+\s+RNAz 2\.1\.1\s+#+", report)
    candidates: list[RnazCandidate] = []

    for chunk in chunks:
        direction_match = re.search(r"Reading direction:\s+(forward|reverse)", chunk)
        if not direction_match:
            continue
        direction = direction_match.group(1)
        coordinate_match = re.search(
            r"^>(?!consensus\b).*?/([0-9]+)-([0-9]+)(?=\s|$)",
            chunk,
            re.MULTILINE,
        )
        if not coordinate_match:
            continue
        start = int(coordinate_match.group(1))
        end = int(coordinate_match.group(2))
        if allowed_windows and (start, end) not in allowed_windows:
            continue

        def value(pattern: str) -> float:
            match = re.search(pattern, chunk, re.MULTILINE)
            return float(match.group(1)) if match else 0.0

        def integer(pattern: str) -> int:
            match = re.search(pattern, chunk, re.MULTILINE)
            return int(match.group(1)) if match else 0

        prediction_match = re.search(
            r"^\s*Prediction:\s*(RNA|OTHER)\s*$",
            chunk,
            re.MULTILINE,
        )
        prediction = prediction_match.group(1) if prediction_match else "OTHER"

        candidates.append(
            RnazCandidate(
                transcript_id=transcript_id,
                start=start,
                end=end,
                direction=direction,
                n_sequences=integer(r"^\s*Sequences:\s+(\d+)"),
                columns=integer(r"^\s*Columns:\s+(\d+)"),
                average_identity=value(r"^\s*Mean pairwise identity:\s+([0-9.eE+-]+)"),
                mean_z_score=value(r"^\s*Mean z-score:\s+([0-9.eE+-]+)"),
                sci=value(r"^\s*Structure conservation index:\s+([0-9.eE+-]+)"),
                svm_value=value(r"^\s*SVM decision value:\s+([0-9.eE+-]+)"),
                rna_probability=value(r"^\s*SVM RNA-class probability:\s+([0-9.eE+-]+)"),
                prediction=prediction,
            )
        )

    unique: dict[tuple[str, int, int, str], RnazCandidate] = {}
    for candidate in candidates:
        unique[(candidate.transcript_id, candidate.start, candidate.end, candidate.direction)] = candidate
    return sorted(
        unique.values(),
        key=lambda candidate: (
            candidate.transcript_id,
            candidate.start,
            candidate.end,
            candidate.direction,
        ),
    )


def select_positive_window_winners(
    candidates: list[RnazCandidate],
    rna_cutoff: float,
) -> list[RnazCandidate]:
    winners: dict[tuple[str, int, int], RnazCandidate] = {}
    for candidate in candidates:
        if candidate.rna_probability < rna_cutoff:
            continue
        key = (candidate.transcript_id, candidate.start, candidate.end)
        current = winners.get(key)
        if (
            current is None
            or candidate.rna_probability > current.rna_probability
            or (
                candidate.rna_probability == current.rna_probability
                and candidate.direction == "forward"
                and current.direction == "reverse"
            )
        ):
            winners[key] = candidate

    return sorted(
        winners.values(),
        key=lambda candidate: (
            candidate.transcript_id,
            candidate.start,
            candidate.end,
        ),
    )


def preprocess_candidate(
    source_alignment: Path,
    start: int,
    end: int,
    candidate_alignment: Path,
    gap_threshold: float,
) -> str:
    extracted = run_tool(
        [
            tool("esl-alimask"),
            "-t",
            "--rna",
            "--outformat",
            "stockholm",
            str(source_alignment),
            f"{start}..{end}",
        ]
    )

    sequences_removed = run_tool(
        [tool("esl-alimanip"), "--lmin", "1", "-"],
        input_text=extracted,
    )

    gap_filtered = run_tool(
        [
            tool("esl-alimask"),
            "-g",
            "--gapthresh",
            str(gap_threshold),
            "--rna",
            "--outformat",
            "stockholm",
            "-",
        ],
        input_text=sequences_removed,
    )

    sequences_removed_again = run_tool(
        [tool("esl-alimanip"), "--lmin", "1", "-"],
        input_text=gap_filtered,
    )

    write_text(candidate_alignment, sequences_removed_again)
    return run_tool(
        [tool("esl-alistat"), "--rna", str(candidate_alignment)]
    )


def parse_rscape_record(
    rscape_dir: Path,
    candidate_name: str,
) -> dict[str, str]:
    record = {
        "rscape_mode": "cacofold",
        "tested_pairs": "0",
        "significant_pairs": "0",
        "observed_covarying_pairs": "0",
        "expected_covarying_pairs": "",
        "mean_power": "",
        "target_evalue": "",
        "coding_flags": "rna,cacofold",
    }

    sorted_cov = rscape_dir / f"{candidate_name}.sorted.cacofold.cov"
    if sorted_cov.is_file():
        text = sorted_cov.read_text()
        lines = [
            line for line in text.splitlines()
            if line.startswith("*") or line.startswith("~")
        ]
        record["tested_pairs"] = str(len(lines))
        record["significant_pairs"] = str(
            sum(line.startswith("*") for line in lines)
        )
        observed_match = re.search(
            r"# BPAIRS observed to covary\s+(\d+)",
            text,
        )
        if observed_match:
            record["observed_covarying_pairs"] = observed_match.group(1)
        target_match = re.search(
            r"Target_E-val\s+\[?[^\]]*\]?\s+\[([^\]]+)\]",
            text,
        )
        if target_match:
            record["target_evalue"] = target_match.group(1)

        powers = []
        for line in lines:
            fields = line.split()
            if len(fields) >= 9:
                try:
                    powers.append(float(fields[-1]))
                except ValueError:
                    pass
        if powers:
            record["mean_power"] = f"{sum(powers) / len(powers):.6f}"

    power_file = rscape_dir / f"{candidate_name}.cacofold.power"
    if power_file.is_file():
        text = power_file.read_text()
        expected_match = re.search(
            r"# BPAIRS expected to covary\s+([0-9.]+)",
            text,
        )
        observed_match = re.search(
            r"# BPAIRS observed to covary\s+(\d+)",
            text,
        )
        if expected_match:
            record["expected_covarying_pairs"] = expected_match.group(1)
        if observed_match:
            record["observed_covarying_pairs"] = observed_match.group(1)

    return record


def classify_rscape(
    status: str,
    pair_count: str,
    record: dict[str, str],
) -> tuple[str, str]:
    if status != "completed":
        return "uninformative", "R-scape assessment missing or failed"

    significant = int(record["significant_pairs"] or 0)
    tested = int(record["tested_pairs"] or 0)
    try:
        mean_power = float(record["mean_power"]) if record["mean_power"] else 0.0
    except ValueError:
        mean_power = 0.0

    if significant > 0:
        return "supported", f"{significant} significant covarying pair(s) reported"
    if tested > 0 and mean_power >= 0.5:
        return (
            "not_supported_with_power",
            f"No significant pairs; mean power across {tested} tested pair(s) is {mean_power:.2f}",
        )
    return "uninformative", "No significant pairs and power/coverage is insufficient or unavailable"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_window_coordinates(windows_alignment: Path) -> set[tuple[int, int]]:
    coordinates: set[tuple[int, int]] = set()
    text = windows_alignment.read_text()
    for match in re.finditer(
        r"^(?:>)?[^\r\n]*[/](\d+)-(\d+)(?:\s|$)",
        text,
        re.MULTILINE,
    ):
        coordinates.add((int(match.group(1)), int(match.group(2))))
    return coordinates


def discover_targets(root: Path, names: list[str] | None) -> list[Path]:
    if names:
        targets = [root / name for name in names]
    else:
        targets = sorted(
            path
            for path in root.glob("AFUB*")
            if path.is_dir()
        )
    for target in targets:
        if not (target / "v3_rm.sto").is_file():
            raise ToolError(f"Missing RepeatMasker alignment: {target / 'v3_rm.sto'}")
    return targets


def analyze_target(
    target_dir: Path,
    rna_cutoff: float,
    rscape_evalue: float,
    gap_threshold: float,
    jobs: int = 1,
) -> tuple[list[RnazCandidate], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    alignment = target_dir / "v3_rm.sto"
    if not alignment.is_file():
        raise ToolError(f"Missing alignment: {alignment}")

    output_dir = target_dir / "analysis_rm"
    output_dir.mkdir(parents=True, exist_ok=True)

    source_stats = run_tool([tool("esl-alistat"), "--rna", str(alignment)])
    write_text(output_dir / "v3_rm.alistat.txt", source_stats)

    clustal_alignment = output_dir / "v3_rm.clustal.aln"
    write_text(
        clustal_alignment,
        run_tool([tool("esl-reformat"), "clustal", str(alignment)]),
    )

    windows_alignment = output_dir / "v3_rm.windows.aln"
    window_result = subprocess.run(
        [
            tool("rnazWindow.pl"),
            "--window",
            str(WINDOW_SIZE),
            "--slide",
            str(WINDOW_SLIDE),
            "--no-reference",
            "--no-rangecheck",
            "--max-seqs",
            str(MAX_SEQUENCES),
            "--opt-id",
            str(OPT_IDENTITY),
            str(clustal_alignment),
        ],
        text=True,
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    write_text(windows_alignment, window_result.stdout)
    write_text(output_dir / "rnazWindow.log", window_result.stderr)

    if window_result.returncode != 0:
        raise ToolError(
            f"rnazWindow.pl failed: {window_result.stderr.strip() or 'no diagnostics'}"
        )

    window_count = window_result.stdout.count("CLUSTAL")
    rnaz_report = output_dir / "v3_rm.rnaz.txt"
    if window_count == 0:
        write_text(
            rnaz_report,
            "No RNAz windows were produced; RNAz was not run.\n",
        )
        candidates: list[RnazCandidate] = []
        reason = "rnazWindow.pl produced no eligible windows"
    else:
        run_tool(
            [
                tool("RNAz"),
                "-b",
                "-o",
                str(rnaz_report),
                str(windows_alignment),
            ]
        )
        allowed_windows = parse_window_coordinates(windows_alignment)
        candidates = parse_rnaz_report(
            rnaz_report,
            target_dir.name,
            allowed_windows,
        )
        reason = ""

    rscape_candidates = select_positive_window_winners(candidates, rna_cutoff)
    valid_candidate_dirs = {
        f"window_{candidate.start}-{candidate.end}_{candidate.direction}"
        for candidate in rscape_candidates
    }
    for stale_dir in output_dir.glob("window_*"):
        if stale_dir.is_dir() and stale_dir.name not in valid_candidate_dirs:
            shutil.rmtree(stale_dir)

    classification_rows: list[dict[str, str]] = []
    rscape_rows: list[dict[str, str]] = []
    machine_rows: list[dict[str, str]] = []

    def analyze_candidate(candidate: RnazCandidate) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
        screen_class = (
            "rnaz_positive"
            if candidate.rna_probability >= rna_cutoff
            else "rnaz_negative"
        )
        candidate_name = (
            f"window_{candidate.start}-{candidate.end}_{candidate.direction}"
        )
        candidate_dir: Path | None = None
        candidate_alignment: Path | None = None
        rscape_dir: Path | None = None

        rscape_status = "not_run"
        rscape_error = ""
        rscape_pair_count = ""
        rscape_record = {
            "rscape_mode": "cacofold",
            "tested_pairs": "0",
            "significant_pairs": "0",
            "observed_covarying_pairs": "0",
            "expected_covarying_pairs": "",
            "mean_power": "",
            "target_evalue": "",
            "coding_flags": "rna,cacofold",
        }

        if screen_class == "rnaz_positive":
            candidate_dir = output_dir / candidate_name
            candidate_dir.mkdir(parents=True, exist_ok=True)
            candidate_alignment = candidate_dir / f"{candidate_name}.sto"
            rscape_dir = candidate_dir / "rscape_cacofold"
            rscape_dir.mkdir(parents=True, exist_ok=True)
            try:
                candidate_stats = preprocess_candidate(
                    alignment,
                    candidate.start,
                    candidate.end,
                    candidate_alignment,
                    gap_threshold,
                )
                write_text(candidate_dir / "esl-alistat.txt", candidate_stats)
                run_tool(
                    [
                        tool("R-scape"),
                        "--cacofold",
                        "--rna",
                        "-E",
                        str(rscape_evalue),
                        "--outdir",
                        str(rscape_dir),
                        "--outname",
                        candidate_name,
                        str(candidate_alignment),
                    ]
                )
                rscape_status = "completed"
                rscape_record = parse_rscape_record(rscape_dir, candidate_name)
                rscape_pair_count = rscape_record["observed_covarying_pairs"]
            except ToolError as error:
                rscape_status = "failed"
                rscape_error = str(error)
                write_text(candidate_dir / "error.txt", rscape_error)

            final_outcome, final_reason = classify_rscape(
                rscape_status,
                rscape_pair_count,
                rscape_record,
            )
        else:
            final_outcome = "uninformative"
            final_reason = (
                "R-scape was not run because this window is RNAz-negative "
                f"at cutoff {rna_cutoff:.6f}"
            )

        classification_row = {
                "transcript_id": candidate.transcript_id,
                "window": candidate.window,
                "direction": candidate.direction,
                "rna_screen_class": screen_class,
                "final_outcome": final_outcome,
                "reason": final_reason,
            }

        rscape_row = {
                "transcript_id": candidate.transcript_id,
                "window": candidate.window,
                "direction": candidate.direction,
                "rna_screen_class": screen_class,
                "final_outcome": final_outcome,
                "reason": final_reason,
            }
        machine_row = {
                "transcript_id": candidate.transcript_id,
                "window": candidate.window,
                "direction": candidate.direction,
                "rna_cutoff": f"{rna_cutoff:.6f}",
                "rna_class_probability": f"{candidate.rna_probability:.6f}",
                "rna_prediction": candidate.prediction,
                "rscape_status": rscape_status,
                "rscape_error": rscape_error,
                "rscape_mode": rscape_record["rscape_mode"],
                "tested_pairs": rscape_record["tested_pairs"],
                "significant_pairs": rscape_record["significant_pairs"],
                "observed_covarying_pairs": rscape_record["observed_covarying_pairs"],
                "expected_covarying_pairs": rscape_record["expected_covarying_pairs"],
                "mean_power": rscape_record["mean_power"],
                "target_evalue": rscape_record["target_evalue"],
                "coding_flags": rscape_record["coding_flags"],
            }
        return classification_row, rscape_row, machine_row

    with ThreadPoolExecutor(max_workers=jobs) as candidate_executor:
        futures = {
            candidate_executor.submit(analyze_candidate, candidate): candidate
            for candidate in rscape_candidates
        }
        completed = []
        for future in as_completed(futures):
            completed.append((futures[future], future.result()))
        completed.sort(
            key=lambda item: (
                item[0].transcript_id,
                item[0].start,
                item[0].end,
                item[0].direction,
            )
        )
        for _, (classification_row, rscape_row, machine_row) in completed:
            classification_rows.append(classification_row)
            rscape_rows.append(rscape_row)
            machine_rows.append(machine_row)

    return candidates, classification_rows, rscape_rows, machine_rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analyze v3_rm Stockholm alignments with RNAz and R-scape CaCoFold."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Directory containing AFUB* target directories and combined output tables",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=None,
        help="AFUB target name; repeat to analyze multiple targets",
    )
    parser.add_argument(
        "--rna-cutoff",
        type=float,
        default=0.01,
        help="Minimum RNAz RNA-class probability for a positive screen",
    )
    parser.add_argument(
        "--rscape-evalue",
        type=float,
        default=0.05,
        help="R-scape E-value target",
    )
    parser.add_argument(
        "--gap-threshold",
        type=float,
        default=0.5,
        help="Maximum fraction of gaps allowed in a column during RNAz-hit preprocessing (default: 0.5)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Maximum parallel target and candidate workers",
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if not 0.0 <= args.gap_threshold <= 1.0:
        parser.error("--gap-threshold must be between 0 and 1")

    try:
        targets = discover_targets(args.root.resolve(), args.target)
    except ToolError as error:
        print(error, file=sys.stderr)
        return 1

    all_candidates: list[RnazCandidate] = []
    all_classification: list[dict[str, str]] = []
    all_rscape_machine: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    failures: list[str] = []

    def run_target(target: Path) -> tuple[list[RnazCandidate], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
        return analyze_target(
            target,
            args.rna_cutoff,
            args.rscape_evalue,
            args.gap_threshold,
            args.jobs,
        )

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_target = {
            executor.submit(run_target, target): target
            for target in targets
        }
        for future in as_completed(future_to_target):
            target = future_to_target[future]
            try:
                candidates, classification, rscape_rows, machine_rows = future.result()
                all_candidates.extend(candidates)
                all_classification.extend(classification)
                all_rscape_machine.extend(machine_rows)
                for row in rscape_rows:
                    summary_rows.append(row)
            except (ToolError, OSError) as error:
                failures.append(f"{target.name}: {error}")
                all_rscape_machine.append(
                    {
                        "transcript_id": target.name,
                        "window": "N/A",
                        "direction": "N/A",
                        "rna_cutoff": f"{args.rna_cutoff:.6f}",
                        "rna_class_probability": "",
                        "rna_prediction": "",
                        "rscape_status": "failed",
                        "rscape_error": str(error),
                        "rscape_mode": "cacofold",
                        "tested_pairs": "",
                        "significant_pairs": "",
                        "observed_covarying_pairs": "",
                        "expected_covarying_pairs": "",
                        "mean_power": "",
                        "target_evalue": "",
                        "coding_flags": "rna,cacofold",
                    }
                )
                print(f"{target.name}: {error}", file=sys.stderr)

    candidate_fields = [
        "Transcript ID",
        "Length (nt)",
        "Window",
        "Direction",
        "NSeq.",
        "Columns",
        "Avg. ID (%)",
        "Mean z-score",
        "SCI",
        "SVM value",
        "RNA-class probability",
        "Prediction",
    ]
    candidate_rows = [
        {
            "Transcript ID": candidate.transcript_id,
            "Length (nt)": str(candidate.end - candidate.start + 1),
            "Window": candidate.window,
            "Direction": candidate.direction,
            "NSeq.": str(candidate.n_sequences),
            "Columns": str(candidate.columns),
            "Avg. ID (%)": f"{candidate.average_identity:.2f}",
            "Mean z-score": f"{candidate.mean_z_score:.2f}",
            "SCI": f"{candidate.sci:.2f}",
            "SVM value": f"{candidate.svm_value:.6f}",
            "RNA-class probability": f"{candidate.rna_probability:.6f}",
            "Prediction": candidate.prediction,
        }
        for candidate in all_candidates
    ]

    output_root = args.root.resolve()
    candidate_table = output_root / "rnaz_candidates.tsv"
    write_tsv(candidate_table, candidate_fields, candidate_rows)

    classification_fields = [
        "Transcript ID",
        "Window",
        "Direction",
        "RNAz screen class",
        "Final outcome",
        "Reason / interpretation",
    ]
    classification_rows = [
        {
            "Transcript ID": row["transcript_id"],
            "Window": row["window"],
            "Direction": row["direction"],
            "RNAz screen class": row["rna_screen_class"],
            "Final outcome": row["final_outcome"],
            "Reason / interpretation": row["reason"],
        }
        for row in all_classification
    ]
    classification_table = output_root / "rscape_classification.tsv"
    write_tsv(classification_table, classification_fields, classification_rows)

    rscape_machine_fields = [
        "transcript_id",
        "window",
        "direction",
        "rna_cutoff",
        "rna_class_probability",
        "rna_prediction",
        "rscape_status",
        "rscape_error",
        "rscape_mode",
        "tested_pairs",
        "significant_pairs",
        "observed_covarying_pairs",
        "expected_covarying_pairs",
        "mean_power",
        "target_evalue",
        "coding_flags",
    ]
    rscape_machine = output_root / "rscape_records.tsv"
    write_tsv(rscape_machine, rscape_machine_fields, all_rscape_machine)

    print(f"Wrote {candidate_table}")
    print(f"Wrote {classification_table}")
    print(f"Wrote {rscape_machine}")

    if failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
