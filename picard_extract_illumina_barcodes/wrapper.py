"""Snakemake wrapper for picard ExtractIlluminaBarcodes."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

lane = snakemake.params["lane"]
read_structure = snakemake.params['read_structure']

extra = snakemake.params.get("extra", "")
output_dir = snakemake.params.get("output_dir", "null")
compress_outputs = snakemake.params.get('compress_outputs', 'true')

shell(
    "picard ExtractIlluminaBarcodes"
    " {extra}"
    " BASECALLS_DIR={snakemake.input.basecalls_dir}"
    " OUTPUT_DIR={output_dir}"
    " LANE={lane}"
    " READ_STRUCTURE={read_structure}"
    " BARCODE_FILE={snakemake.input.barcode_file}"
    " METRICS_FILE={snakemake.output.metrics_file}"
    " COMPRESS_OUTPUTS={compress_outputs}"
    " NUM_PROCESSORS={snakemake.threads}"
    " {log}")

