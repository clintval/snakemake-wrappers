"""Snakemake wrapper for dwgsim."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell


def make_params(params):
    import types
    formatted_params = ''

    def clean_value(value):
        if value is True:
            return 'true'
        elif value is False:
            return 'false'
        elif value is None:
            return 'null'
        elif isinstance(value, (list, tuple, types.GeneratorType)):
            return list(map(clean_value, value))
        else:
            return value

    def make_key(key):
        return f'-{key}'

    for key, value in params.items():
        if key == 'extra':
            continue
        value = clean_value(value)
        if isinstance(value, list):
            formatted_params += ''.join(f' {make_key(key)} {v}' for v in value)
        else:
            formatted_params += f' {make_key(key)} {value}'
    return formatted_params

extra = snakemake.params.get('extra', '')
params = make_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'dwgsim'
    ' {extra}'
    ' {params}'
    ' {snakemake.input}'
    ' {snakemake.output}'
    ' {log}')

