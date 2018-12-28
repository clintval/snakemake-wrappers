from snakescale.testing import run_tool_test


def test_bedtools_subtract():
    process_return = run_tool_test('bedtools', 'subtract')
    assert process_return.exit_code.is_ok(), process_return.stderr
