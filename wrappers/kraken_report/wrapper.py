"""Snakemake wrapper for Kraken report."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

def make_kraken_params(params):
    formatted_params = ''

    for key, value in params.items():
        if key == 'extra':
            continue

        if value is True:
            formatted_params += f' --{key.replace("_", "-")}'
        elif value is False:
            continue
        else:
            formatted_params += f' --{key.replace("_", "-")} {value}'
    return formatted_params

params = make_kraken_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get('extra', '')

shell(
    'kraken-report'
    ' {extra}'
    ' {params}'
    ' {snakemake.input}'
    ' > {snakemake.output}'
    ' {log}'
)
