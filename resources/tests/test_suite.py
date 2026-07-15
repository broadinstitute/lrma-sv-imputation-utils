"""Optional pytest entry point.

Thin wrapper so `pytest` discovers and runs the same shell suites CI runs. It
does not reimplement anything -- each test shells out to a per-tool run.sh and
asserts a zero exit status. Binary discovery honours $BIN_DIR (default
target/release) exactly as the shell harness does.

    pip install pytest
    BIN_DIR=target/debug pytest resources/tests/test_suite.py -v
"""

import os
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))

SUITES = {
    "extract_bubble_pls": "extract-bubble-PLs/run.sh",
    "pop_glimpse2": "pop-glimpse2/run.sh",
    "paste_vcfs": "paste-vcfs/run.sh",
}


@pytest.mark.parametrize("name,script", list(SUITES.items()), ids=list(SUITES))
def test_suite(name, script):
    proc = subprocess.run(
        ["bash", os.path.join(HERE, script)],
        cwd=HERE, capture_output=True, text=True,
        env={**os.environ},
    )
    if proc.returncode != 0:
        raise AssertionError(
            "suite %s failed (exit %d)\n--- stdout ---\n%s\n--- stderr ---\n%s"
            % (name, proc.returncode, proc.stdout, proc.stderr)
        )
