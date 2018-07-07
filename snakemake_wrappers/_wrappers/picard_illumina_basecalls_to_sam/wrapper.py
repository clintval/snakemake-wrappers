"""Snakemake wrapper for picard IlluminaBasecallsToSam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'picard IlluminaBasecallsToSam'
    ' {extra}'
    ' BASECALLS_DIR={snakemake.input.basecalls_dir}'
    ' BARCODES_DIR={snakemake.input.barcodes_dir}'
    ' LIBRARY_PARAMS={snakemake.input.library_params}'
    ' NUM_PROCESSORS={snakemake.threads}'
    ' {params}'
    ' {log}')
