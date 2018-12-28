from pathlib import Path
from typing import List, Optional

from snakescale import formatters
from snakescale import testing
from snakescale import utils
from snakescale import resources

__all__ = [
    'DEFAULT_SCALE_FILENAME',
    'DEFAULT_SCALE_VERSION',
    'AVAILABLE_WRAPPER_VERSIONS',
    'WRAPPER_ROOT_PATH',
    'available_tools',
    'available_subtools',
    'scale',
]

DEFAULT_SCALE_FILENAME = 'wrapper.py'
DEFAULT_SCALE_VERSION = 'latest'
WRAPPER_ROOT_PATH = Path(__file__).parent / 'wrappers'

AVAILABLE_WRAPPER_VERSIONS = set(map(lambda _: _.name, WRAPPER_ROOT_PATH.glob('*')))


def available_tools(version: str = DEFAULT_SCALE_VERSION) -> List[str]:
    if version not in AVAILABLE_WRAPPER_VERSIONS:
        raise ValueError(f'Version does not exist or is not supported: {version}')
    tools = sorted([path.name for path in (WRAPPER_ROOT_PATH / version).glob('*')])
    return tools


def available_subtools(tool: str, version: str = DEFAULT_SCALE_VERSION) -> Optional[List[str]]:
    if tool not in available_tools(version=version):
        raise ValueError(f'Tool {tool} does not exist for wrappers version {version}')
    subtools = sorted([path.name for path in (WRAPPER_ROOT_PATH / version / tool).glob('*')])
    if DEFAULT_SCALE_FILENAME in subtools:
        return None
    else:
        return subtools


def scale(
    tool: str,
    subtool: Optional[str] = None,
    version: str = DEFAULT_SCALE_VERSION,
    as_uri: bool = True,
) -> str:
    if version not in AVAILABLE_WRAPPER_VERSIONS:
        raise ValueError(f'Wrapper version is not valid or supported: {version}')

    wrapper_dir = WRAPPER_ROOT_PATH / version / tool
    wrapper_dir = wrapper_dir / subtool if subtool is not None else wrapper_dir
    wrapper_file = wrapper_dir / DEFAULT_SCALE_FILENAME

    assert wrapper_file.exists(), f'Wrapper does not exist: {wrapper_file}'
    return wrapper_dir.as_uri() if as_uri is True else str(wrapper_dir)
