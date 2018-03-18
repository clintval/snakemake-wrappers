"""Snakemake wrapper for dwgsim."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

def make_params(params):
    formatted_params = ''

    for key, value in params.items():
        if key == 'extra':
            continue

        if key == 'r1':
            key = '1'
        elif key == 'r2':
            key = '2'

        if value is True:
            formatted_params += f' -{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' -{key} {value}'
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

