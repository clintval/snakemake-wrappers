"""Snakemake wrapper for picard MergeSamFiles."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.utils import collect_jvm_resources
from snakescale.utils import collect_picard_style_jvm_resources

extra = snakemake.params.get('extra', '')
extra += collect_jvm_resources()
extra += collect_picard_style_jvm_resources()
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if isinstance(snakemake.input, dict):
    input_files = snakemake.input['bams']
    intervals = snakemake.input.pop('intervals', 'null')
else:
    input_files = snakemake.input
    intervals = 'null'

if isinstance(input_files, (list, tuple)):
    input_files = ''.join(f' INPUT={bam}' for bam in input_files)
else:
    input_files = f' INPUT={input_files}'


shell(
    'picard MergeSamFiles'
    ' {extra}'
    ' {input_files}'
    ' INTERVALS={intervals}'
    ' OUTPUT={snakemake.output}'
    ' {params}'
    ' {log}'
)
