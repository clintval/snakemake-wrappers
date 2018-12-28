__all__ = ['collect_jvm_resources', 'collect_picard_style_jvm_resources']  # type: ignore


def collect_jvm_resources() -> str:
    """Collect JVM resources allocations the resource dictionary."""
    args = ''
    if snakemake.resources.get('gc_heap_free_limit'):  # type: ignore
        args += f' -XX:GCHeapFreeLimit={snakemake.resources.gc_heap_free_limit}'  # type: ignore

    if snakemake.resources.get('gc_time_limit'):  # type: ignore
        args += f' -XX:GCTimeLimit={snakemake.resources.gc_time_limit}'  # type: ignore

    if snakemake.resources.get('heap_size'):  # type: ignore
        args += f' -Xmx{snakemake.resources.heap_size}m'  # type: ignore

    return args


def collect_picard_style_jvm_resources() -> str:
    """Collect Picard specific JVM resources allocations the resource dictionary."""
    args = ''
    if snakemake.resources.get('samjdk_buffer_size'):  # type: ignore
        args += f' -Dsamjdk.buffer_size={snakemake.resources.samjdk_buffer_size}'  # type: ignore

    if snakemake.resources.get('use_async_io_read_samtools') == 1:  # type: ignore
        args += ' -Dsamjdk.use_async_io_read_samtools=true'

    if snakemake.resources.get('use_async_io_write_samtools') == 1:  # type: ignore
        args += ' -Dsamjdk.use_async_io_write_samtools=true'

    return args
