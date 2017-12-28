"""Snakemake wrapper for picard BedToIntervalList."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

sort = snakemake.params.get('sort', 'true')
unique = snakemake.params.get('unique', 'false')

sequence_dictionary = snakemake.input.get('sequence_dictionary', 'null')

shell(
    "picard BedToIntervalList"
    " INPUT={snakemake.input[0]}"
    " OUTPUT={snakemake.output}"
    " SEQUENCE_DICTIONARY={sequence_dictionary}"
    " SORT={sort}"
    " UNIQUE={unique}")
