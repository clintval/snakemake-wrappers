"""Snakemake wrapper for picard IlluminaBasecallsToSam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

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

if snakemake.resources.get('malloc'):
    extra += f' -Xmx{snakemake.resources.malloc}m'

shell(
    'picard IlluminaBasecallsToSam'
    ' {extra}'
    ' BASECALLS_DIR={snakemake.input.basecalls_dir}'
    ' BARCODES_DIR={snakemake.input.barcodes_dir}'
    ' LIBRARY_PARAMS={snakemake.input.library_params}'
    ' NUM_PROCESSORS={snakemake.threads}'
    ' {params}'
    ' {log}')
