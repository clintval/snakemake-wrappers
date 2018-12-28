from snakescale.testing import run_tool_test


def test_bedtools_subtract():
    run_tool_test('bedtools', 'subtract')
