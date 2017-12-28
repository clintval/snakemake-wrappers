"""Snakemake wrapper for picard CollectQualityYieldMetrics."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

input_files = "".join(f" INPUT={bam}" for bam in snakemake.input)

extra = snakemake.params.get("extra", "")
use_original_qualities = snakemake.params.get("use_original_qualities", True)
include_secondary_alignments = snakemake.params.get("include_secondary_alignments", False)
include_supplemental_alignments = snakemake.params.get("include_supplemental_alignments", False)
assume_sorted = snakemake.params.get("assume_sorted", True)
stop_after = snakemake.params.get("stop_after", 0)

if use_original_qualities != "null":
    use_original_qualities = "true" if use_original_qualities else "false"
if include_secondary_alignments != "null":
    include_secondary_alignments = "true" if include_secondary_alignments else "false"
if include_supplemental_alignments != "null":
    include_supplemental_alignments = "true" if include_supplemental_alignments else "false"
if assume_sorted != "null":
    assume_sorted = "true" if assume_sorted else "false"

shell(
    "picard GatherBamFiles"
    " {input_files}"
    " OUTPUT=/dev/stdout"
    " {log} | "
    "picard CollectQualityYieldMetrics"
    " {extra}"
    " USE_ORIGINAL_QUALITIES={use_original_qualities}"
    " INCLUDE_SECONDARY_ALIGNMENTS={include_secondary_alignments}"
    " INCLUDE_SUPPLEMENTAL_ALIGNMENTS={include_supplemental_alignments}"
    " INPUT=/dev/stdin"
    " OUTPUT={snakemake.output}"
    " ASSUME_SORTED={assume_sorted}"
    " STOP_AFTER={stop_after}"
    " {log}")
