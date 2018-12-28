from types import GeneratorType
from typing import List, Mapping, Union

__all__ = [
    'clean_picard_style_value',
    'snakecase_to_kebab_case',
    'clean_picard_style_key',
    'format_bedtools_params',
    'format_bwa_params',
    'format_dwgsim_params',
    'format_fgbio_params',
    'format_kraken_params',
    'format_picard_params',
]


def clean_picard_style_value(value: Union[List[str], str]) -> Union[List[str], str]:
    """Clean a dictionary of Picard key-value pairs."""
    if isinstance(value, (list, tuple, GeneratorType)):
        return list(map(clean_picard_style_value, value))  # type: ignore
    elif value is None:
        return 'null'
    elif value is True:
        return 'true'
    elif value is False:
        return 'false'
    else:
        return value


def format_bed_key(key: str) -> str:
    """Clean a bedtools parameter key."""
    return '-' + key.replace('_', '')


def snakecase_to_kebab_case(key: str) -> str:
    """Convert snake_case to kebab-case."""
    return f'--{key.lower().replace("_", "-")}'


def clean_picard_style_key(key: str) -> str:
    """Clean a Picard parameter key."""
    return key.upper()


def format_bedtools_params(params: Mapping) -> str:
    """Clean a dictionary of bedtools key-value pairs."""
    formatted_params = ''

    for key, value in params.items():

        if key == 'extra':
            continue

        key = format_bed_key(key)
        if value is True:
            formatted_params += f' {key}'
        elif value is False:
            continue
        else:
            formatted_params += f' {key} {value}'
    return formatted_params


def format_bwa_params(params: Mapping) -> str:
    """Clean a dictionary of bwa key-value pairs."""
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


def format_dwgsim_params(params: Mapping) -> str:
    """Clean a dictionary of dwgsim key-value pairs."""
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


def format_fgbio_params(params: Mapping) -> str:
    """Clean a dictionary of fgbio key-value pairs."""
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


def format_kraken_params(params: Mapping) -> str:
    """Clean a dictionary of kraken key-value pairs."""
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


def format_picard_params(params: Mapping) -> str:
    """Clean a dictionary of picard key-value pairs."""
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
