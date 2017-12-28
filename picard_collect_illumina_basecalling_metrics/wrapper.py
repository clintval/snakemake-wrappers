"""Snakemake wrapper for picard CollectIlluminaBasecallingMetrics."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

lane = snakemake.params["lane"]
read_structure = snakemake.params["read_structure"]

extra = snakemake.params.get("extra", "")
barcodes_dir = snakemake.input.get("barcodes_dir", "null")
barcode_file = snakemake.input.get("barcode_file", "null")

shell(
    "picard CollectIlluminaBasecallingMetrics"
    " {extra}"
    " BASECALLS_DIR={snakemake.input.basecalls_dir}"
    " BARCODES_DIR={barcodes_dir}"
    " LANE={lane}"
    " INPUT={barcode_file}"
    " READ_STRUCTURE={read_structure}"
    " OUTPUT={snakemake.output}"
    " {log}")
