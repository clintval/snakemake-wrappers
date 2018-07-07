"""Snakemake wrapper for picard FilterSamReads."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell
from snakemake_wrappers.utils import collect_jvm_resources
from snakemake_wrappers.utils import collect_picard_style_jvm_resources
from snakemake_wrappers.utils import make_fgbio_params

extra = snakemake.params.get('extra', '')
extra += collect_jvm_resources()
extra += collect_picard_style_jvm_resources()
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.input.get('bam'):
    params += f' INPUT={snakemake.input.bam}'

    if snakemake.input.get('read_list_file'):
        params += f' READ_LIST_FILE={snakemake.input.read_list_file}'

    if snakemake.input.get('interval_list'):
        params += f' INTERVAL_LIST={snakemake.input.interval_list}'

    if snakemake.input.get('javascript_file'):
        params += f' JAVASCRIPT_FILE={snakemake.input.javascript_file}'
else:
    params += f' INPUT={snakemake.input}'


shell(
    'picard FilterSamReads'
    ' {extra}'
    ' OUTPUT={snakemake.output}'
    ' {params}'
    ' {log}')
