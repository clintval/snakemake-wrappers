"""Snakemake wrapper for picard SortVcf."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "picard"
    " SortVcf"
    " {extra}"
    " I={snakemake.input.vcf}"
    " O={snakemake.output}"
    " SD={snakemake.input.sequence_dictionary}"
    " {log}")
