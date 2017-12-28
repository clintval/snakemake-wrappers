"""Snakemake wrapper for picard MergeSamFiles."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


if isinstance(snakemake.input, dict):
    bams = snakemake.input["bams"]
    intervals = snakemake.input.pop("intervals", "null")
else:
    bams = snakemake.input
    intervals = "null"

if isinstance(bams, (list, tuple)):
    INPUT = "".join(f" INPUT={bam}" for bam in bams)
else:
    INPUT = f" INPUT={bams}"

extra = snakemake.params.get("extra", "")
sort_order = snakemake.params.get("sort_order", "coordinate")
assume_sorted = snakemake.params.get("assume_sorted", False)
merge_sequence_dictionaries = snakemake.params.get("merge_sequence_dictionaries", False)
use_threading = snakemake.params.get("use_threading", False)
comment = snakemake.params.get("comment", "null")

if assume_sorted != "null":
    assume_sorted = True if assume_sorted else False
if merge_sequence_dictionaries != "null":
    merge_sequence_dictionaries = True if merge_sequence_dictionaries else False
if use_threading != "null":
    use_threading = True if use_threading else False

shell(
    "picard MergeSamFiles"
    " {extra}"
    " {INPUT}"
    " OUTPUT={snakemake.output}"
    " SORT_ORDER={sort_order}"
    " ASSUME_SORTED={assume_sorted}"
    " MERGE_SEQUENCE_DICTIONARIES={merge_sequence_dictionaries}"
    " USE_THREADING={use_threading}"
    " COMMENT={comment}"
    " INTERVALS={intervals}"
    " {log}")
