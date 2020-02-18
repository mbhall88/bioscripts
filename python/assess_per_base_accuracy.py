import click
from collections import Counter
import logging
from pathlib import Path
import pysam
import json


class NoBamIndex(Exception):
    pass


def validate_file(ctx, param, value):
    if value.exists():
        return value
    else:
        raise click.BadParameter(f"{value} does not exist!")


def pileup_column_agrees_with_reference(
    column: pysam.PileupColumn, ref_base: str, quorum: float
) -> bool:
    """quorum says that this percent of bases in the column must agree
    with the reference base for it to be considered consensus."""
    base_counts = Counter(map(str.upper, column.get_query_sequences()))
    consensus = base_counts.most_common(1)

    if not consensus:
        return False

    consensus_base, consensus_count = consensus[0]
    consensus_ratio = consensus_count / column.n * 100
    if consensus_base != ref_base.upper() or consensus_ratio < quorum:
        return False

    return True


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "--bam",
    help="Bam file to assess.",
    type=Path,
    required=True,
    callback=validate_file,
)
@click.option(
    "--assembly",
    help="Fasta of assembly bed file is related to.",
    type=Path,
    required=True,
    callback=validate_file,
)
@click.option(
    "--outfile", help="Path to write output json file to.", type=Path, required=True
)
@click.option(
    "--quorum",
    help="Percentage of bases that must agree with reference at each position.",
    type=float,
    default=95.0,
    show_default=True,
)
def main(bam: Path, assembly: Path, quorum: float, outfile: Path):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

    with pysam.FastxFile(assembly) as fastx:
        assembly_index = {entry.name: entry.sequence for entry in fastx}

    num_bases = sum(len(seq) for seq in assembly_index.values())
    logging.info(f"{num_bases} bases found in assembly.")

    base_disagreements = 0
    positions_covered_by_pileup = 0
    mean_position_mapping_qual = []
    all_mapping_qualities = []
    with pysam.AlignmentFile(bam) as alignment:
        alignment.check_index()

        for column in alignment.pileup():
            ref_base = assembly_index[column.reference_name][column.reference_pos]
            if not pileup_column_agrees_with_reference(column, ref_base, quorum):
                base_disagreements += 1

            mapping_qualities = column.get_mapping_qualities()
            mean_mapq = sum(mapping_qualities) / column.n
            mean_position_mapping_qual.append(mean_mapq)
            positions_covered_by_pileup += 1

        for record in alignment:
            all_mapping_qualities.append(record.mapping_quality)

    base_disagreements += (num_bases - positions_covered_by_pileup)

    logging.info(
        f"There are {base_disagreements} positions that disagree with the assembly."
    )

    perc_disagree = base_disagreements / num_bases * 100

    logging.info(
        f"{perc_disagree}% of the assembly positions disagree with the reads in the BAM."
    )

    mean_per_base_mapping_quality = sum(mean_position_mapping_qual) / num_bases
    logging.info(f"Mean per-base mapping quality is {mean_per_base_mapping_quality}")
    mean_total_mapping_quality = sum(all_mapping_qualities) / len(all_mapping_qualities)
    logging.info(f"Overall mean mapping quality is {mean_total_mapping_quality}")

    with outfile.open("w") as handle:
        print(
            json.dumps(
                {
                    "base_disagreements": base_disagreements,
                    "total_bases": num_bases,
                    "percent_bases_disagree": perc_disagree,
                    "mean_per_base_mapping_quality": mean_per_base_mapping_quality,
                    "mean_per_read_mapping_quality": mean_total_mapping_quality,
                },
                sort_keys=True,
                indent=4,
            ),
            file=handle,
        )


if __name__ == "__main__":
    main()
