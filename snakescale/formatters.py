from types import GeneratorType as _generator_type
from typing import List, Mapping, Union


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


def clean_picard_style_value(value) -> Union[List[str], str]:
    if isinstance(value, (list, tuple, _generator_type)):
        return list(map(clean_picard_style_value, value))
    elif value is None:
        return 'null'
    elif value is True:
        return 'true'
    elif value is False:
        return 'false'
    else:
        return value


def snakecase_to_kebab_case(key: str) -> str:
    return f'--{key.lower().replace("_", "-")}'


def clean_picard_style_key(key: str) -> str:
    return key.upper()


def make_bwa_params(params: Mapping) -> str:
    formatted_params = ''

    for key, value in params.items():

        if key == 'extra':
            continue
        elif value is True:
            formatted_params += f' -{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' -{key} {value}'
    return formatted_params


def make_dwgsim_params(params: Mapping) -> str:
    formatted_params = ''

    for key, value in params.items():
        if key in ('extra', 'output_prefix'):
            continue

        key = '1' if key == 'r1' else key
        key = '2' if key == 'r2' else key

        if value is True:
            formatted_params += f' -{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' -{key} {value}'
    return formatted_params


def make_fgbio_params(params: Mapping) -> str:
    formatted_params = ''

    for key, value in params.items():
        key = snakecase_to_kebab_case(key)
        value = clean_picard_style_value(value)

        if key == 'extra':
            continue
        elif isinstance(value, list):
            formatted_params += ''.join(f' --{key}={v}' for v in value)
        else:
            formatted_params += f' --{key}={value}'
    return formatted_params


def make_kraken_params(params: Mapping) -> str:
    formatted_params = ''

    for key, value in params.items():
        key = snakecase_to_kebab_case(key)

        if key == 'extra':
            continue
        elif value is True:
            formatted_params += f' --{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' --{key} {value}'
    return formatted_params


def make_picard_params(params: Mapping) -> str:
    formatted_params = ''

    for key, value in params.items():
        key = clean_picard_style_key(key)
        value = clean_picard_style_value(value)

        if key == 'extra':
            continue
        elif isinstance(value, list):
            formatted_params += ''.join(f' {key}={v}' for v in value)
        else:
            formatted_params += f' {key}={value}'
    return formatted_params
