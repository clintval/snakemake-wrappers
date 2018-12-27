from types import GeneratorType as _generator_type
from typing import List, Mapping, Union

__all__ = [
    'clean_picard_style_value',
    'snakecase_to_kebab_case',
    'clean_picard_style_key',
    'make_bwa_params',
    'make_dwgsim_params',
    'make_fgbio_params',
    'make_kraken_params',
    'make_picard_params',
]


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
