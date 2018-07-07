"""Snakemake wrapper for picard CollectAlignmentSummaryMetrics."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'picard CollectAlignmentSummaryMetrics'
    ' {extra}'
    ' INPUT={snakemake.input.bam}'
    ' OUTPUT={snakemake.output}'
    ' REFERENCE_SEQUENCE={snakemake.input.reference}'
    ' {params}'
    ' {log}')
