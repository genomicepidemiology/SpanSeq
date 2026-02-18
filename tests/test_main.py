"""Tests for spanseq.__main__ module."""
import pytest


class TestMainModule:
    def test_main_module_runs(self, mocker):
        """Invoking __main__ calls cli.main()."""
        mock_main = mocker.patch("spanseq.cli.main")
        import runpy
        runpy.run_module("spanseq", run_name="__main__")
        mock_main.assert_called_once()
