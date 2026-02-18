"""Tests for the FASTA parser and writer."""
from spanseq.fasta import parse_fasta, write_fasta, FastaRecord


def test_parse_multi_sequence(sample_fasta):
    records = parse_fasta(sample_fasta)
    assert len(records) == 3
    assert records[0].id == "seq1"
    assert records[0].description == "seq1 gene product"
    assert records[0].seq == "ATCGATCGATCG"


def test_parse_multiline_sequence(sample_fasta):
    records = parse_fasta(sample_fasta)
    assert records[1].id == "seq2"
    assert records[1].seq == "MLLKPPAAVVGGCC"


def test_parse_id_with_slash(sample_fasta):
    records = parse_fasta(sample_fasta)
    assert records[2].id == "seq3/variant"
    assert records[2].seq == "ACGTACGT"


def test_parse_empty_file(tmp_path):
    empty = tmp_path / "empty.fsa"
    empty.write_text("")
    assert parse_fasta(empty) == []


def test_write_parse_roundtrip(tmp_path, sample_fasta):
    original = parse_fasta(sample_fasta)

    out = tmp_path / "roundtrip.fsa"
    with open(out, "w") as fh:
        for r in original:
            write_fasta(fh, r)

    roundtripped = parse_fasta(out)
    assert len(roundtripped) == len(original)
    for a, b in zip(original, roundtripped):
        assert a.id == b.id
        assert a.description == b.description
        assert a.seq == b.seq


def test_parse_description_with_spaces(tmp_path):
    fasta = tmp_path / "desc.fsa"
    fasta.write_text(">gene_1 length=100 organism=E.coli\nATCG\n")
    records = parse_fasta(fasta)
    assert records[0].id == "gene_1"
    assert records[0].description == "gene_1 length=100 organism=E.coli"


def test_parse_single_sequence(tmp_path):
    fasta = tmp_path / "single.fsa"
    fasta.write_text(">only\nGGCC\n")
    records = parse_fasta(fasta)
    assert len(records) == 1
    assert records[0].id == "only"
    assert records[0].seq == "GGCC"
