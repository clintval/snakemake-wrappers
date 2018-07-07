"""Snakemake wrapper for picard CollectQualityYieldMetrics."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

input_files = ''.join(f' INPUT={bam}' for bam in snakemake.input)

shell(
    'picard GatherBamFiles'
    ' {input_files}'
    ' OUTPUT=/dev/stdout'
    ' {log}'
    ' | picard CollectQualityYieldMetrics'
    ' {extra}'
    ' {params}'
    ' INPUT=/dev/stdin'
    ' OUTPUT={snakemake.output}'
    ' {log}')
