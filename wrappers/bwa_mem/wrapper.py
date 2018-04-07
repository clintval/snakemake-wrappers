"""Snakemake wrapper for bwa mem ClipBam."""

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

def make_bwa_params(params):
    formatted_params = ''

    for key, value in params.items():
        if key == 'extra':
            continue

        if value is True:
            formatted_params += f' -{key}'
        elif value is False:
            continue
        else:
            formatted_params += f' -{key} {value}'
    return formatted_params

sam_to_fastq_params = make_picard_params(snakemake.params.get('sam_to_fastq', {}))
bwa_mem_params = make_bwa_params(snakemake.params.get('bwa_mem', {}))
merge_bam_alignment_params = make_picard_params(snakemake.params.get('merge_bam_alignment', {}))

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

shell(
    'picard SamToFastq'
    ' INPUT={snakemake.input.unmapped_bam}'
    ' FASTQ=/dev/stdout'
    ' {sam_to_fastq_params}'
    ' {log}'
    ' | bwa mem'
    ' -t {snakemake.threads}'
    ' {snakemake.input.reference}'
    ' /dev/stdin'
    ' {bwa_mem_params}'
    ' {log}'
    ' | picard MergeBamAlignment '
    ' UNMAPPED={snakemake.input.unmapped_bam}'
    ' REFERENCE_SEQUENCE={snakemake.input.reference}'
    ' ALIGNED=/dev/stdin'
    ' OUTPUT={snakemake.output}'
    ' {merge_bam_alignment_params}'
    ' {log}')
