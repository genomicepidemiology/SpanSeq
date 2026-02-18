"""Tests for CdHitApp and select_word_size."""
import pytest
from spanseq.applications.cdhit import CdHitApp, select_word_size


class TestSelectWordSize:
    def test_nucleotide_high_threshold(self):
        assert select_word_size(0.95, "nucleotides") == 8

    def test_nucleotide_mid_threshold(self):
        assert select_word_size(0.85, "nucleotides") == 6

    def test_nucleotide_low_threshold(self):
        assert select_word_size(0.75, "nucleotides") == 4

    def test_protein_high_threshold(self):
        assert select_word_size(0.80, "aminoacids") == 5

    def test_protein_low_threshold(self):
        assert select_word_size(0.40, "aminoacids") == 2

    def test_too_low_raises(self):
        with pytest.raises(ValueError, match="too low"):
            select_word_size(0.30, "aminoacids")


class TestCdHitApp:
    def test_nucleotide_uses_est(self, fake_exec_pair):
        est, aa = fake_exec_pair
        app = CdHitApp(exec_path_est=est, exec_path_aa=aa)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_file="/tmp/out.fsa",
            seq_type="nucleotides",
            threshold=0.9,
            threads=4,
        )
        assert cmd[0] == str(est)
        assert "-c" in cmd
        idx = cmd.index("-c")
        assert cmd[idx + 1] == "0.9"

    def test_aminoacid_uses_cdhit(self, fake_exec_pair):
        est, aa = fake_exec_pair
        app = CdHitApp(exec_path_est=est, exec_path_aa=aa)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_file="/tmp/out.fsa",
            seq_type="aminoacids",
            threshold=0.7,
        )
        assert cmd[0] == str(aa)

    def test_auto_word_size(self, fake_exec_pair):
        est, aa = fake_exec_pair
        app = CdHitApp(exec_path_est=est, exec_path_aa=aa)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_file="/tmp/out.fsa",
            seq_type="nucleotides",
            threshold=0.9,
        )
        idx = cmd.index("-n")
        assert cmd[idx + 1] == "8"  # word_size for >= 0.90

    def test_explicit_word_size(self, fake_exec_pair):
        est, aa = fake_exec_pair
        app = CdHitApp(exec_path_est=est, exec_path_aa=aa)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_file="/tmp/out.fsa",
            seq_type="nucleotides",
            threshold=0.9,
            word_size=5,
        )
        idx = cmd.index("-n")
        assert cmd[idx + 1] == "5"
