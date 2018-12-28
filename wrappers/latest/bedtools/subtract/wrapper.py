"""Snakemake wrapper for bedtools subtract."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.formatters import format_bedtools_params

extra = snakemake.params.get('extra', '')
params = format_bedtools_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

command = (
    'bedtools subtract'
    ' {extra}'
    ' {params}'
    ' -a {snakemake.input.a}'
    ' -b {snakemake.input.b}'
    ' > {snakemake.output} {log}'
)

shell(command)
