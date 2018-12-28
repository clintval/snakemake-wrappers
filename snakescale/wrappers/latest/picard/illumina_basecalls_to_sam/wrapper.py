"""Snakemake wrapper for picard IlluminaBasecallsToSam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.utils import collect_jvm_resources
from snakescale.utils import collect_picard_style_jvm_resources

extra = snakemake.params.get('extra', '')
extra += collect_jvm_resources()
extra += collect_picard_style_jvm_resources()
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
    ' {log}'
)
