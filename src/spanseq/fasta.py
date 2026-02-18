"""FASTA parser and writer — no BioPython dependency."""
from __future__ import annotations
from collections import namedtuple
from pathlib import Path
from typing import List

FastaRecord = namedtuple("FastaRecord", ["id", "description", "seq"])


def parse_fasta(path: Path | str) -> List[FastaRecord]:
    """Parse a FASTA file into a list of records."""
    records: List[FastaRecord] = []
    with open(path) as fh:
        header = None
        seq_parts: List[str] = []
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    desc = header[1:]
                    seq_id = desc.split()[0] if desc else ""
                    records.append(FastaRecord(seq_id, desc, "".join(seq_parts)))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            desc = header[1:]
            seq_id = desc.split()[0] if desc else ""
            records.append(FastaRecord(seq_id, desc, "".join(seq_parts)))
    return records


def write_fasta(handle, record: FastaRecord) -> None:
    """Write a single FASTA record to an open file handle."""
    handle.write(f">{record.description}\n{record.seq}\n")
