"""Snakemake wrapper for Kraken classify."""
# Requires this patch to Kraken
# https://github.com/DerrickWood/kraken/pull/106
__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

import uuid

from pathlib import Path

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

def make_kraken_params(params):
    formatted_params = ''

    for key, value in params.items():
        if key == 'extra':
            continue

        if value is True:
            formatted_params += f' --{key.replace("_", "-")}'
        elif value is False:
            continue
        else:
            formatted_params += f' --{key.replace("_", "-")} {value}'
    return formatted_params

sam_to_fastq_params = make_picard_params(snakemake.params.get('sam_to_fastq', {}))
kraken_params = make_kraken_params(snakemake.params.get('kraken', {}))
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

fastq1 = Path(str(snakemake.input)).stem + str(uuid.uuid4())
fastq2 = Path(str(snakemake.input)).stem + str(uuid.uuid4())

shell(
    'mkfifo /tmp/{fastq1} && mkfifo /tmp/{fastq2} &&'
    ' picard SamToFastq'
    ' INPUT={snakemake.input}'
    ' FASTQ=/tmp/{fastq1}'
    ' SECOND_END_FASTQ=/tmp/{fastq2}'
    ' {sam_to_fastq_params}'
    ' {log}'
    ' | kraken'
    ' --threads {snakemake.threads}'
    ' {kraken_params}'
    ' /tmp/{fastq1}'
    ' /tmp/{fastq2}'
    ' > {snakemake.output}'
    ' {log}'
)
