"""Snakemake wrapper for picard CollectIlluminaLaneMetrics."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output_directory = snakemake.params['output_directory']
output_prefix = snakemake.params['output_prefix']

extra = snakemake.params.get("extra", "")
read_structure = snakemake.params.get('read_structure', 'null')
file_extension = snakemake.params.get('file_extension', '.txt')
is_novaseq = snakemake.params.get('is_novaseq', 'false')
assert is_novaseq in ('true', 'false'), '`is_novaseq` must be "true" or "false"'

shell(
    "picard CollectIlluminaLaneMetrics"
    " {extra}"
    " RUN_DIRECTORY={snakemake.input}"
    " OUTPUT_DIRECTORY={output_directory}"
    " OUTPUT_PREFIX={output_prefix}"
    " READ_STRUCTURE={read_structure}"
    " FILE_EXTENSION={file_extension}"
    " IS_NOVASEQ={is_novaseq}"
    " {log}")
