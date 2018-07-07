"""Snakemake wrapper for picard CollectBaseDistributionByCycle."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'picard CollectBaseDistributionByCycle'
    ' {extra}'
    ' INPUT={snakemake.input}'
    ' OUTPUT={snakemake.output.metrics_file}'
    ' CHART_OUTPUT={snakemake.output.chart_output}'
    ' {params}'
    ' {log}')
