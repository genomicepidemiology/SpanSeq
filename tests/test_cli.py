"""Tests for the SpanSeq CLI argument parser and main entry point."""
import pytest
import sys
from spanseq.cli import build_parser, main


@pytest.fixture
def parser():
    return build_parser()


class TestSplitParsing:
    def test_minimal_split(self, parser):
        args = parser.parse_args([
            "split", "-i", "input.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "5",
        ])
        assert args.action == "split"
        assert args.input_fasta == "input.fsa"
        assert args.seqtype == "nucleotides"
        assert args.min_dist == 0.3
        assert args.bins == "5"

    def test_split_defaults(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.5", "-b", "3",
        ])
        assert args.distanceMethod == "cosine"
        assert args.approach == "all"
        assert args.output_format == "minimal"
        assert args.threads == 1

    def test_split_all_distance_methods(self, parser):
        for method in ("jaccard", "szymkiewicz_simpson", "cosine",
                        "kmer_inv", "mash", "identity", "mmseqs2",
                        "mmseqs-fast"):
            args = parser.parse_args([
                "split", "-i", "x.fsa", "-s", "nucleotides",
                "-o", "out/", "-c", "0.3", "-b", "5",
                "-d", method,
            ])
            assert args.distanceMethod == method

    def test_split_tool_paths(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "5",
            "-KP", "/usr/bin/kma",
            "-SP", "/usr/bin/mmseqs",
            "-GP", "/usr/bin/ggsearch36",
        ])
        assert args.kmaPath == "/usr/bin/kma"
        assert args.mmseqsPath == "/usr/bin/mmseqs"
        assert args.GGSearchPath == "/usr/bin/ggsearch36"


class TestTreeParsing:
    def test_tree_flag(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "3", "--tree",
        ])
        assert args.tree is True

    def test_tree_default_off(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "3",
        ])
        assert args.tree is False

    def test_tree_method(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "3",
            "--tree", "--tree_method", "nj",
        ])
        assert args.tree_method == "nj"

    def test_tree_method_default(self, parser):
        args = parser.parse_args([
            "split", "-i", "x.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "3",
        ])
        assert args.tree_method == "dnj"


class TestReduceParsing:
    def test_minimal_reduce(self, parser):
        args = parser.parse_args([
            "reduce", "-i", "input.fsa", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "5",
        ])
        assert args.action == "reduce"


class TestInputGroup:
    def test_mutually_exclusive(self, parser):
        with pytest.raises(SystemExit):
            parser.parse_args([
                "split", "-i", "a.fsa", "-if", "folder/",
                "-s", "nucleotides", "-o", "out/", "-c", "0.3", "-b", "5",
            ])

    def test_folder_input(self, parser):
        args = parser.parse_args([
            "split", "-if", "seqs/", "-s", "aminoacids",
            "-o", "out/", "-c", "0.3", "-b", "5",
        ])
        assert args.input_folder == "seqs/"
        assert args.input_fasta is None

    def test_batch_input(self, parser):
        args = parser.parse_args([
            "split", "-ib", "batch.txt", "-s", "nucleotides",
            "-o", "out/", "-c", "0.3", "-b", "5",
        ])
        assert args.input_batch == "batch.txt"


class TestNoAction:
    def test_no_action_returns_none(self, parser):
        args = parser.parse_args([])
        assert args.action is None


class TestMain:
    def test_no_action_exits(self, mocker):
        mocker.patch("sys.argv", ["spanseq"])
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1

    def test_successful_split(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        fasta.write_text(">s1\nATCG\n")
        out = tmp_path / "output"

        mocker.patch("sys.argv", [
            "spanseq", "split",
            "-i", str(fasta), "-s", "nucleotides",
            "-o", str(out), "-c", "0.3", "-b", "3",
        ])

        mock_pipeline = mocker.MagicMock()
        mock_pipeline.run.return_value = {"clusters": "c.tsv", "makespan": "m.tsv"}
        mocker.patch("spanseq.cli.SpanSeqPipeline", return_value=mock_pipeline)

        main()
        mock_pipeline.run.assert_called_once()

    def test_successful_with_list_output(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        fasta.write_text(">s1\nATCG\n")
        out = tmp_path / "output"

        mocker.patch("sys.argv", [
            "spanseq", "split",
            "-i", str(fasta), "-s", "nucleotides",
            "-o", str(out), "-c", "0.3", "-b", "3",
        ])

        mock_pipeline = mocker.MagicMock()
        mock_pipeline.run.return_value = {
            "fasta_files": ["f1.fsa", "f2.fsa", "f3.fsa"],
        }
        mocker.patch("spanseq.cli.SpanSeqPipeline", return_value=mock_pipeline)

        main()
        mock_pipeline.run.assert_called_once()

    def test_file_not_found_exits(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        fasta.write_text(">s1\nATCG\n")

        mocker.patch("sys.argv", [
            "spanseq", "split",
            "-i", str(fasta), "-s", "nucleotides",
            "-o", str(tmp_path / "out"), "-c", "0.3", "-b", "3",
        ])
        mocker.patch(
            "spanseq.cli.SpanSeqPipeline",
            side_effect=FileNotFoundError("ccphylo not found"),
        )

        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1

    def test_value_error_exits(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        fasta.write_text(">s1\nATCG\n")

        mocker.patch("sys.argv", [
            "spanseq", "split",
            "-i", str(fasta), "-s", "nucleotides",
            "-o", str(tmp_path / "out"), "-c", "0.3", "-b", "3",
        ])
        mocker.patch(
            "spanseq.cli.SpanSeqPipeline",
            side_effect=ValueError("bad config"),
        )

        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1

    def test_unexpected_error_exits(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        fasta.write_text(">s1\nATCG\n")

        mocker.patch("sys.argv", [
            "spanseq", "split",
            "-i", str(fasta), "-s", "nucleotides",
            "-o", str(tmp_path / "out"), "-c", "0.3", "-b", "3",
        ])
        mocker.patch(
            "spanseq.cli.SpanSeqPipeline",
            side_effect=RuntimeError("unexpected"),
        )

        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
