"""Snakemake wrapper for fgbio ClipBam."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
metric_prefix = snakemake.params.get('metric_prefix', None)
soft_clip = snakemake.params.get('soft_clip', 'false')
clipping_mode = snakemake.params.get('clipping_mode', None)
auto_clip_attributes = snakemake.params.get('auto_clip_attributes', 'false')
upgrade_clipping = snakemake.params.get('upgrade_clipping', 'false')
read_one_five_prime = snakemake.params.get('read_one_five_prime', 0)
read_one_three_prime = snakemake.params.get('read_one_three_prime', 0)
read_two_five_prime = snakemake.params.get('read_two_five_prime', 0)
read_two_three_prime = snakemake.params.get('read_two_three_prime', 0)
clip_overlapping_reads = snakemake.params.get('clip_overlapping_reads', 'false')

command = (
    "fgbio ClipBam"
    " {extra}"
    " -i {snakemake.input.bam}"
    " -o {snakemake.output[0]}"
    " -r {snakemake.input.reference}"
    " --auto-clip-attributes={auto_clip_attributes}"
    " --upgrade-clipping={upgrade_clipping}"
    " --read-one-five-prime={read_one_five_prime}"
    " --read-one-three-prime={read_one_three_prime}"
    " --read-two-five-prime={read_two_five_prime}"
    " --read-two-three-prime={read_two_three_prime}"
    " --clip-overlapping-reads={clip_overlapping_reads}"
    " {log}")

if metric_prefix is not None:
    command += " --metrics={metric_prefix}"

if clipping_mode is not None:
    command += " --clipping-mode={clipping_mode}"
else:
    command += " --soft-clip={soft_clip}"

shell(command)
