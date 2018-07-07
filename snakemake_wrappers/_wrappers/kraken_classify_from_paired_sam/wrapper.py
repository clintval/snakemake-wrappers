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

sam_to_fastq_params = make_picard_params(snakemake.params.get('sam_to_fastq', {}))
kraken_params = make_kraken_params(snakemake.params.get('kraken', {}))
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

fastq1 = Path(str(snakemake.input)).stem + str(uuid.uuid4())
fastq2 = Path(str(snakemake.input)).stem + str(uuid.uuid4())

shell(
    'mkfifo /tmp/{fastq1} && mkfifo /tmp/{fastq2} &&'
    ' picard -Xmx5g SamToFastq'
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
