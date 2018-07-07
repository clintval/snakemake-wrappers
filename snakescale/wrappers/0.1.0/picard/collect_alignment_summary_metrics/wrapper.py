"""Snakemake wrapper for picard CollectAlignmentSummaryMetrics."""

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

shell(
    'picard CollectAlignmentSummaryMetrics'
    ' {extra}'
    ' INPUT={snakemake.input.bam}'
    ' OUTPUT={snakemake.output}'
    ' REFERENCE_SEQUENCE={snakemake.input.reference}'
    ' {params}'
    ' {log}'
)
