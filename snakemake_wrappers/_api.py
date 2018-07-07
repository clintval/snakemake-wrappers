from pathlib import Path
from typing import List

import snakemake_wrappers as sw

__all__ = [
    'list_wrappers',
    'root_path',
    'wrappers'
]


root_path: Path = Path(sw.__file__).resolve().parent / '_wrappers'


class _Wrappers(object):
    def __init__(self):
        for wrapper in list_wrappers():
            setattr(self, wrapper, root_path / wrapper)


def list_wrappers() -> List:
    return sorted(map(lambda _: _.name, root_path.glob('*')))


wrappers = _Wrappers()
