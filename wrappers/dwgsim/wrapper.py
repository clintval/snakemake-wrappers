"""Snakemake wrapper for dwgsim."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

def make_params(params):
    formatted_params = ''

    for key, value in params.items():
        if key in ('extra', 'output_prefix'):
            continue

        key = '1' if key == 'r1' else key
        key = '2' if key == 'r2' else key

        if value is True:
            formatted_params += f' -{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' -{key} {value}'
    return formatted_params

output_prefix = snakemake.params["output_prefix"]

extra = snakemake.params.get('extra', '')
params = make_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'dwgsim'
    ' {extra}'
    ' {params}'
    ' {snakemake.input}'
    ' {output_prefix}'
    ' {log}')
