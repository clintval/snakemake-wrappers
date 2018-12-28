from datetime import datetime, timedelta

from hypothesis import example, given
from hypothesis.strategies import integers

from snakescale.utils import ExitCode


@given(integers())
@example(0)
def test_exit_code(integer):
    exit_code = ExitCode(integer)

    if integer == 0:
        assert exit_code.is_ok()
    if integer != 0:
        assert not exit_code.is_ok()

    assert exit_code == ExitCode(integer)
    assert repr(exit_code) == f'ExitCode({integer})'
