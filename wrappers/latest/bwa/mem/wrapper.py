"""Snakemake wrapper for bwa mem ClipBam."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.utils import collect_jvm_resources
from snakescale.utils import collect_picard_style_jvm_resources
from snakescale.utils import make_bwa_params
from snakescale.utils import make_picard_params

extra = snakemake.params.get('extra', '')
extra += collect_jvm_resources()
extra += collect_picard_style_jvm_resources()
sam_to_fastq_params = make_picard_params(snakemake.params.get('sam_to_fastq', {}))
bwa_mem_params = make_bwa_params(snakemake.params.get('bwa_mem', {}))
merge_bam_alignment_params = make_picard_params(snakemake.params.get('merge_bam_alignment', {}))
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

shell(
    'picard {extra} SamToFastq'
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
    ' | picard {extra} MergeBamAlignment '
    ' UNMAPPED={snakemake.input.unmapped_bam}'
    ' REFERENCE_SEQUENCE={snakemake.input.reference}'
    ' ALIGNED=/dev/stdin'
    ' OUTPUT={snakemake.output}'
    ' {merge_bam_alignment_params}'
    ' {log}'
)
