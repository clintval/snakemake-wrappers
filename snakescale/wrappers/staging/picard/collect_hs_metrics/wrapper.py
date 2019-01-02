"""Snakemake wrapper for picard CollectHsMetrics."""

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

bait_intervals = snakemake.input['bait_intervals']
target_intervals = snakemake.input.get('target_intervals', 'null')

if isinstance(bait_intervals, (list, tuple)):
    bait_intervals = ''.join(f' BAIT_INTERVALS={interval}' for interval in bait_intervals)
else:
    bait_intervals = f' BAIT_INTERVALS={bait_intervals}'

if isinstance(target_intervals, (list, tuple)):
    target_intervals = ''.join(f' TARGET_INTERVALS={interval}' for interval in target_intervals)
elif target_intervals == 'null':
    target_intervals = bait_intervals.replace('BAIT_INTERVALS', 'TARGET_INTERVALS')
else:
    target_intervals = f' TARGET_INTERVALS={target_intervals}'

summary_output = snakemake.output.get('summary_output', snakemake.output[0])
per_target_coverage = snakemake.output.get('per_target_coverage', 'null')
per_base_coverage = snakemake.output.get('per_base_coverage', 'null')

shell(
    'picard CollectHsMetrics'
    ' {extra}'
    ' INPUT={snakemake.input.bam}'
    ' OUTPUT={summary_output}'
    ' REFERENCE_SEQUENCE={snakemake.input.reference}'
    ' PER_TARGET_COVERAGE={per_target_coverage}'
    ' PER_BASE_COVERAGE={per_base_coverage}'
    ' {bait_intervals}'
    ' {target_intervals}'
    ' {params}'
    ' {log}'
)
