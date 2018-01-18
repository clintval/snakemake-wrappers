"""Snakemake wrapper for fgbio ClipBam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell


def make_fgbio_params(params):
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

    def clean_key(key):
        return f'--{key.lower().replace("_", "-")}'

    for key, value in params.items():
        if key == 'extra':
            continue
        value = clean_value(value)
        if isinstance(value, list):
            formatted_params += ''.join(f' {clean_key(key)}={v}' for v in value)
        else:
            formatted_params += f' {clean_key(key)}={value}'
    return formatted_params

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    'fgbio ClipBam'
    ' {extra}'
    ' -i {snakemake.input.bam}'
    ' -o {snakemake.output[0]}'
    ' -r {snakemake.input.reference}'
    ' {params}'
    ' {log}')

