"""Shared fixtures for application runner tests."""
import pytest
from pathlib import Path


@pytest.fixture
def fake_exec(tmp_path):
    """Create a fake executable for ApplicationRunner validation."""
    exe = tmp_path / "fake_tool"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    return exe


@pytest.fixture
def fake_exec_pair(tmp_path):
    """Create two fake executables (for CdHitApp: cd-hit-est + cd-hit)."""
    est = tmp_path / "cd-hit-est"
    est.write_text("#!/bin/sh\nexit 0\n")
    est.chmod(0o755)
    aa = tmp_path / "cd-hit"
    aa.write_text("#!/bin/sh\nexit 0\n")
    aa.chmod(0o755)
    return est, aa
