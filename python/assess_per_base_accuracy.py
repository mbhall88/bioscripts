import json
import logging
from collections import Counter, defaultdict
from enum import Enum
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import List, Tuple

import click
import pandas as pd
import pysam

# Pileup is 1-based positions, whereas BED is 0-based positions
PILEUP_TO_BED_OFFSET = 1


class NoBamIndex(Exception):
    pass


class PileupColumn:
    def __init__(
        self,
        chromosome: str,
        ref_pos: int,
        ref_base: str,
        depth: int,
        read_bases: str,
        base_qualities: str,
        mapping_qualities: str = "",
    ):
        self.chromosome = chromosome
        self.ref_pos = ref_pos
        self.ref_base = ref_base
        self.depth = depth
        self.read_bases = read_bases
        self.base_qualities = base_qualities
        self.mapping_qualities = mapping_qualities

    @staticmethod
    def from_string(line: str) -> "PileupColumn":
        fields = line.rstrip().split()
        mapping_qualities = "" if len(fields) < 7 else fields[6]
        return PileupColumn(
            chromosome=fields[0],
            ref_pos=int(fields[1]),
            ref_base=fields[2],
            depth=int(fields[3]),
            read_bases=fields[4],
            base_qualities=fields[5],
            mapping_qualities=mapping_qualities,
        )

    def matches(self) -> int:
        base_counts = Counter(self.read_bases)
        matches_forward = base_counts["."]
        matches_reverse = base_counts[","]
        return matches_forward + matches_reverse

    def match_ratio(self) -> float:
        try:
            return self.matches() / self.depth
        except ZeroDivisionError:
            return 0


def validate_file(ctx, param, value):
    if value.exists():
        return value
    else:
        raise click.BadParameter(f"{value} does not exist!")


def pileup_column_agrees_with_reference(column: PileupColumn, quorum: int) -> bool:
    """quorum says that this percent of bases in the column must agree
    with the reference base for it to be considered consensus."""
    match_percent = column.match_ratio() * 100
    logging.debug(
        f"Position {column.ref_pos} on {column.chromosome} has {column.matches()} "
        f"matches out of {column.depth} reads covering the site ("
        f"{round(match_percent, 2)}% match)."
    )
    return match_percent >= quorum


class Interval(Enum):
    OPEN = (1, 1)
    CLOSED = (0, 0)
    RIGHT_OPEN = (0, 1)
    LEFT_OPEN = (1, 0)


def collapse_positions_into_intervals(
    data: List[int], interval_type: Interval = Interval.RIGHT_OPEN
) -> List[Tuple[int, int]]:
    left_offset, right_offset = interval_type.value
    ranges = []
    for k, g in groupby(enumerate(data), lambda x: x[0] - x[1]):
        group = map(itemgetter(1), g)
        group = list(map(int, group))
        ranges.append((group[0] - left_offset, group[-1] + right_offset))

    return ranges


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "--bam",
    help=(
        "Bam file to assess. It is up to the user whether to include "
        "secondary/supplementary alignments. These can be removed with "
        "`samtools view -bh -f 0 -F 256`"
    ),
    type=Path,
    required=True,
    callback=validate_file,
)
@click.option(
    "--pileup",
    help=(
        "Pileup file generated by samtools mpileup. The pileup MUST be generated with "
        "the --fasta-ref option - with the assembly being provided as the fasta ref. "
        "It is also recommended that the pileup is generated using the mpileup flag "
        "-aa to report absolutely all positions."
    ),
    type=Path,
    required=True,
    callback=validate_file,
)
@click.option(
    "-p",
    "--prefix",
    help="Path prefix to write output JSON and BED files to.",
    type=str,
    required=True,
)
@click.option(
    "--quorum",
    help="Percentage of reads that must agree with the assembly at each position.",
    type=click.IntRange(0, 100),
    default=95,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(bam: Path, pileup: Path, quorum: int, prefix: str, verbose: bool):
    """A script to assess the per-base quality of an assembly. There are two key
    metrics analysed here:\n
    1. Per-base consensus of reads mapped to the assembly. Consensus is defined as
    whether the percentage of reads matching the assembly at a site is
    greater than or equal to the value passed to --quorum\n
    2. Mapping quality summary statistics. i.e. mean, median, quantiles etc.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(format="%(levelname)s: %(message)s", level=log_level)

    pileup_stats = dict()

    logging.info("Assessing pileup...")
    disagreement_positions = defaultdict(list)
    num_disagreements = 0
    with pileup.open() as pileup_handle:
        for num_pileup_positions, line in enumerate(pileup_handle, start=1):
            column = PileupColumn.from_string(line)
            if not pileup_column_agrees_with_reference(column, quorum):
                num_disagreements += 1
                disagreement_positions[column.chromosome].append(
                    column.ref_pos - PILEUP_TO_BED_OFFSET
                )

    logging.info(f"There are {num_disagreements} disagreements in total.")
    pileup_stats["total_disagreements"] = num_disagreements
    logging.info(f"There are {num_pileup_positions} total positions in the pileup.")
    pileup_stats["total_pileup_positions"] = num_pileup_positions
    percent_disagree = num_disagreements / num_pileup_positions * 100
    logging.info(
        f"Therefore, {round(percent_disagree, 2)}% of positions did not reach quorum."
    )
    pileup_stats["percent_pileup_disagree"] = percent_disagree

    logging.info("Assessing mapping quality...")
    mapping_qualities = []
    with pysam.AlignmentFile(bam) as alignment:
        for num_bam_entries, record in enumerate(alignment, start=1):
            mapping_qualities.append(record.mapping_quality)

    mapq_summary = pd.Series(mapping_qualities).describe()
    logging.info(f"Mapping quality summary statistics:\n{mapq_summary.to_string()}")

    logging.info("Writing output files...")
    json_path = Path(prefix + f".json")

    summary_stats = dict()
    summary_stats["pileup_stats"] = {**pileup_stats}
    summary_stats["mapping_quality_stats"] = {**mapq_summary.to_dict()}
    with json_path.open("w") as json_out_handle:
        print(json.dumps(summary_stats, sort_keys=True, indent=4), file=json_out_handle)
    logging.info(f"Summary statistics written to {json_path}")

    bed_path = Path(prefix + ".bed")
    with bed_path.open("w") as bed_out_handle:
        for chromosome, positions in disagreement_positions:
            for chrom_start, chrom_end in collapse_positions_into_intervals(positions):
                print(f"{chromosome}\t{chrom_start}\t{chrom_end}", file=bed_out_handle)
    logging.info(f"Positions that disagree with the assembly are written to {bed_path}")

    logging.info("Finished!")


if __name__ == "__main__":
    main()
