"""Shared fixtures for SpanSeq tests."""
import os
import pytest
from pathlib import Path


@pytest.fixture
def fake_exec(tmp_path):
    """Create a fake executable that does nothing."""
    exe = tmp_path / "fake_tool"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    return exe


@pytest.fixture
def sample_fasta(tmp_path):
    """Write a small 3-sequence FASTA and return its path."""
    fasta = tmp_path / "sample.fsa"
    fasta.write_text(
        ">seq1 gene product\n"
        "ATCGATCGATCG\n"
        ">seq2\n"
        "MLLKPPAA\n"
        "VVGGCC\n"
        ">seq3/variant\n"
        "ACGTACGT\n"
    )
    return fasta


@pytest.fixture
def sample_fasta_path(tmp_path):
    """Create a FASTA file and return the tmp_path (for config tests needing a real file)."""
    fasta = tmp_path / "input.fsa"
    fasta.write_text(">s1\nATCG\n")
    return fasta
