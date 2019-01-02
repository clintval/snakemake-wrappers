from snakescale.testing import run_tool_test
from snakescale import scale

from pathlib import Path


def test_bedtools_subtract():
    process_return = run_tool_test('bedtools', 'subtract')
    assert process_return.exit_code.is_ok(), process_return.stderr

def test_dwgsim():
    process_return = run_tool_test('dwgsim')
    assert process_return.exit_code.is_ok(), process_return.stderr
