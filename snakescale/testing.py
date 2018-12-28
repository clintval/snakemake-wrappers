from pathlib import Path
from typing import Optional

from snakemake.utils import available_cpu_count

from snakescale.utils import ProcessReturn, call_snakemake

TEST_ARGUMENTS_BASE = ''
TEST_ARGUMENTS_BASE += ' -F'
TEST_ARGUMENTS_BASE += ' --use-conda'
TEST_ARGUMENTS_BASE += ' --printshellcmds'
TEST_ARGUMENTS_BASE += f' --cores {available_cpu_count()}'


def run_tool_test(
    tool: str, subtool: Optional[str] = None, version: str = 'latest'
) -> ProcessReturn:
    """Execute a tool and optional subtool test against a specific wrapper version."""
    from snakescale import scale

    wrapper_dir = Path(scale(tool, subtool, version, as_uri=False)) / 'test'
    arguments = TEST_ARGUMENTS_BASE
    arguments += f' --directory {wrapper_dir}'
    arguments += f' --snakefile {wrapper_dir}/Snakefile'
    process_return = call_snakemake(arguments)
    return process_return
