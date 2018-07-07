"""Snakemake wrapper for picard CollectIlluminaBasecallingMetrics."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'picard CollectIlluminaBasecallingMetrics'
    ' {extra}'
    ' INPUT={snakemake.input.barcode_file}'
    ' BARCODES_DIR={snakemake.input.barcodes_dir}'
    ' BASECALLS_DIR={snakemake.input.basecalls_dir}'
    ' OUTPUT={snakemake.output}'
    ' {params}'
    ' {log}')
