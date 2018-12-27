__all__ = ['collect_jvm_resources', 'collect_picard_style_jvm_resources']


def collect_jvm_resources() -> str:
    args = ''
    if snakemake.resources.get('gc_heap_free_limit'):
        args += f' -XX:GCHeapFreeLimit={snakemake.resources.gc_heap_free_limit}'

    if snakemake.resources.get('gc_time_limit'):
        args += f' -XX:GCTimeLimit={snakemake.resources.gc_time_limit}'

    if snakemake.resources.get('heap_size'):
        args += f' -Xmx{snakemake.resources.heap_size}m'

    return args


def collect_picard_style_jvm_resources() -> str:
    args = ''
    if snakemake.resources.get('samjdk_buffer_size'):
        args += f' -Dsamjdk.buffer_size={snakemake.resources.samjdk_buffer_size}'

    if snakemake.resources.get('use_async_io_read_samtools') == 1:
        args += ' -Dsamjdk.use_async_io_read_samtools=true'

    if snakemake.resources.get('use_async_io_write_samtools') == 1:
        args += ' -Dsamjdk.use_async_io_write_samtools=true'

    return args
