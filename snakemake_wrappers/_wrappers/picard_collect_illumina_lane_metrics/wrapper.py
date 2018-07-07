"""Snakemake wrapper for picard CollectIlluminaLaneMetrics."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from pathlib import Path
from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
params = make_picard_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not snakemake.params.get('output_directory'):
    params += f' OUTPUT_DIRECTORY={Path(snakemake.output.lane_metrics).parent}'
if not snakemake.params.get('output_prefix'):
    params += f' OUTPUT_PREFIX={Path(Path(snakemake.output.lane_metrics).stem).stem}'

shell(
    'picard CollectIlluminaLaneMetrics'
    ' {extra}'
    ' RUN_DIRECTORY={snakemake.input}'
    ' {params}'
    ' {log}')
