"""Snakemake wrapper for fgbio ClipBam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell
from snakemake_wrappers import _make_fgbio_params

extra = snakemake.params.get('extra', '')
params = _make_fgbio_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'fgbio ClipBam'
    ' {extra}'
    ' -i {snakemake.input.bam}'
    ' -o {snakemake.output[0]}'
    ' -r {snakemake.input.reference}'
    ' {params}'
    ' {log}'
)
