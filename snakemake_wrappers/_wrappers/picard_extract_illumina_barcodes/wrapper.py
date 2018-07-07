"""Snakemake wrapper for picard ExtractIlluminaBarcodes."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'picard ExtractIlluminaBarcodes'
    ' {extra}'
    ' BARCODE_FILE={snakemake.input.barcode_file}'
    ' BASECALLS_DIR={snakemake.input.basecalls_dir}'
    ' METRICS_FILE={snakemake.output.metrics_file}'
    ' OUTPUT_DIR={snakemake.output.barcodes_dir}'
    ' NUM_PROCESSORS={snakemake.threads}'
    ' {params}'
    ' {log}')

