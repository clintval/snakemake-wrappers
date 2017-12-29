"""Snakemake wrapper for bedtools subtract."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    'bedtools subtract'
    ' -a {snakemake.input.a}'
    ' -b {snakemake.input.b}'
    ' > {snakemake.output} {log}')
