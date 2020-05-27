#!/usr/bin/env python3
import logging
from itertools import repeat
from typing import TextIO, Tuple, Dict, List
from dataclasses import dataclass
from enum import Enum

import click

Contig = str
Seq = str
Index = Dict[Contig, Seq]


class Strand(Enum):
    Forward = "+"
    Reverse = "-"
    NotRelevant = "."
    Unknown = "?"

    def __str__(self) -> str:
        return str(self.value)


@dataclass
class GffFeature:
    seqid: Contig
    source: str
    method: str  # correct term is type, but that is a python reserved variable name
    start: int
    end: int
    score: float
    strand: Strand
    phase: int
    attributes: Dict[str, str]

    @staticmethod
    def from_str(s: str) -> "GffFeature":
        fields = s.split("\t")
        score = 0 if fields[5] == "." else float(fields[5])
        phase = -1 if fields[7] == "." else int(fields[7])
        attr_fields = fields[-1].split(";")
        attributes = {k: v for k, v in map(str.split, attr_fields, repeat("="))}
        return GffFeature(
            seqid=fields[0],
            source=fields[1],
            method=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=score,
            strand=Strand(fields[6]),
            phase=phase,
            attributes=attributes,
        )


class DuplicateContigsError(Exception):
    pass


def is_header(s: str) -> bool:
    if not s:
        return False
    return s[0] == ">"


def index_fasta(stream: TextIO) -> Index:
    fasta_index: Index = dict()
    sequence: List[Seq] = []
    name: Contig = ""
    for line in map(str.rstrip, stream):
        if not line:
            continue
        if is_header(line):
            if sequence and name:
                fasta_index[name] = "".join(sequence)
                sequence = []
            name = line.split()[0][1:]
            if name in fasta_index:
                raise DuplicateContigsError(
                    f"Contig {name} occurs multiple times in the fasta file."
                )
            continue
        else:
            sequence.append(line)
    if name and sequence:
        fasta_index[name] = "".join(sequence)

    return fasta_index


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-f",
    "--fasta",
    help="FASTA file to split.",
    type=click.File(mode="r"),
    default="-",
    show_default=True,
    required=True,
)
@click.option(
    "-g",
    "--gff",
    help="GFF3 file to base split coordinates on.",
    type=click.File(mode="r"),
    required=True,
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(file_okay=False, resolve_path=True, writable=True),
    default=".",
    show_default=True,
    help="The directory to write the output files to.",
)
@click.option(
    "--types",
    help=(
        "The feature types to split on. Separate types by a space or pass option "
        "mutiple times."
    ),
    multiple=True,
    default=["gene"],
    show_default=True,
)
@click.option(
    "--min-igr-len",
    help="The minimum length of the intergenic regions to output.",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--max-igr-len",
    help="The maximum length of the intergenic regions to output. Set to 0 to disable "
    "IGR output.",
    type=float,
    default=float("inf"),
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    fasta: TextIO,
    gff: TextIO,
    outdir: str,
    types: Tuple[str],
    min_igr_len: int,
    max_igr_len: float,
    verbose: bool,
):
    """Splits a FASTA file into chunks based on a GFF3 file.
    The splits produced are based on the --types given and everything inbetween. For
    example, the default --types is 'gene'. In this case, the coordinates for each gene
    are cut out of the FASTA file, as well as the bits inbetween - intergenic regions
    (IGRs).
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    logging.info("Indexing fasta file...")
    index: Index = index_fasta(fasta)
    logging.info(f"{len(index)} contigs indexed in the input file.")

    for feature in map(GffFeature.from_str, gff):
        # todo


if __name__ == "__main__":
    main()
