"""Snakemake wrapper for picard CollectIlluminaLaneMetrics."""

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

if not snakemake.params.get('output_directory'):
    params += f' OUTPUT_DIRECTORY={Path(snakemake.output.lane_metrics).parent}'
if not snakemake.params.get('output_prefix'):
    params += f' OUTPUT_PREFIX={Path(Path(snakemake.output.lane_metrics).stem).stem}'

shell(
    'picard CollectIlluminaLaneMetrics'
    ' {extra}'
    ' RUN_DIRECTORY={snakemake.input}'
    ' {params}'
    ' {log}'
)
