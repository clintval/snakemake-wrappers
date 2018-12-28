import io
import os
import sys
from typing import NamedTuple, Type

from snakemake import main as cli

__all__ = ['ExitCode', 'ProcessReturn', 'call_snakemake']


class ExitCode(int):
    """The code returned to a parent process by an executable.

    Examples:
        >>> from subprocess import call
        >>> exit_code = ExitCode(call('ls'))
        >>> exit_code.is_ok()
        True

    """

    def __new__(cls, exit_code: int) -> Type['ExitCode']:
        """Make a new :class:`ExitCode`."""
        cls._code = exit_code
        return super().__new__(cls, exit_code)  # type: ignore

    def is_ok(self) -> bool:
        """Is this code zero."""
        return self._code == 0

    def __repr__(self) -> str:
        """Represent this object."""
        return f'{self.__class__.__qualname__}({self._code})'


class ProcessReturn(NamedTuple):
    """A collection of STDOUT, STDERR, and exit code from a terminated process."""

    stdout: str
    stderr: str
    exit_code: ExitCode


def call_snakemake(arguments: str) -> ProcessReturn:
    """Call Snakemake in this thread as if we are using the CLI.

    Args:
        arguments: A string of arguments to pass to Snakemake.

    Examples:
        >>> process_return = call_snakemake('')
        >>> process_return.exit_code.is_ok()
        False
        >>> process_return.stderr.strip()
        'Error: Snakefile "Snakefile" not present.'

    """
    current_environ = dict(os.environ)
    current_stdout = sys.stdout
    current_stderr = sys.stderr

    exit_code = ExitCode(0)
    stdout = ''
    stderr = ''

    try:
        env = os.environ.copy()
        env['LC_CTYPE'] = u'en_US.UTF'
        os.environ.update(env)

        sys.stdout = io.StringIO()
        sys.stderr = io.StringIO()

        try:
            cli(arguments)
        except SystemExit as e:
            exit_code = ExitCode(e.code)

        stdout = sys.stdout.getvalue()
        stderr = sys.stderr.getvalue()
    finally:
        sys.stdout = current_stdout
        sys.stderr = current_stderr
        os.environ.clear()
        os.environ.update(current_environ)

    process_return = ProcessReturn(stdout=stdout, stderr=stderr, exit_code=exit_code)
    return process_return
