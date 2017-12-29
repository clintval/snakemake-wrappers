"""Snakemake wrapper for picard CollectHsMetrics."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

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

def make_picard_params(params):
    import types
    formatted_params = ''

    def clean_value(value):
        if value is True:
            return 'true'
        elif value is False:
            return 'false'
        elif value is None:
            return 'null'
        elif isinstance(value, (list, tuple, types.GeneratorType)):
            return list(map(clean_value, value))
        else:
            return value

    for key, value in params.items():
        if key == 'extra':
            continue
        value = clean_value(value)
        if isinstance(value, list):
            formatted_params += ''.join(f' {key.upper()}={v}' for v in value)
        else:
            formatted_params += f' {key.upper()}={value}'
    return formatted_params

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.resources.get('gc_heap_free_limit'):
    extra += f' -XX:GCHeapFreeLimit={snakemake.resources.gc_heap_free_limit}'

if snakemake.resources.get('gc_time_limit'):
    extra += f' -XX:GCTimeLimit={snakemake.resources.gc_time_limit}'

if snakemake.resources.get('malloc'):
    extra += f' -Xmx{snakemake.resources.malloc}m'

if snakemake.resources.get('samjdk_buffer_size'):
    extra += f' -Dsamjdk.buffer_size={snakemake.resources.samjdk_buffer_size}'

if snakemake.resources.get('use_async_io_read_samtools') == 1:
    extra += ' -Dsamjdk.use_async_io_read_samtools=true'

if snakemake.resources.get('use_async_io_write_samtools') == 1:
    extra += ' -Dsamjdk.use_async_io_write_samtools=true'

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
    ' {log}')
