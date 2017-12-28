"""Snakemake wrapper for picard IlluminaBasecallsToSam."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

lane = snakemake.params["lane"]
run_barcode = snakemake.params["run_barcode"]
read_structure = snakemake.params["read_structure"]

extra = snakemake.params.get("extra", "")
barcodes_dir = snakemake.input.get("barcodes_dir", "null")
read_group_id = snakemake.params.get("read_group_id", "null")
sequencing_center = snakemake.params.get("sequencing_center", "null")
run_start_date = snakemake.params.get("run_start_date", "null")
platform = snakemake.params.get("platform", "null")
adapters_to_check = snakemake.params.get("adapters_to_check", "null")
five_prime_adapter = snakemake.params.get("five_prime_adapter", "null")
three_prime_adapter = snakemake.params.get("three_prime_adapter", "null")
first_tile = snakemake.params.get("first_tile", "null")
time_limit = snakemake.params.get("time_limit", "null")
force_gc = snakemake.params.get("force_gc", False)
apply_eamss_filter = snakemake.params.get("apply_eamss_filter", True)
max_reads_in_ram_per_tile = snakemake.params.get("max_reads_in_ram_per_tile", 1200000)
minimum_quality = snakemake.params.get("minimum_quality", 2)
include_non_pf_reads = snakemake.params.get("include_non_pf_reads", True)
ignore_unexpected_barcodes = snakemake.params.get("ignore_unexpected_barcodes", False)
molecular_index_tag = snakemake.params.get("molecular_index_tag", "RX")
molecular_index_base_quality_tag = snakemake.params.get("molecular_index_base_quality_tag", "QX")
tag_per_molecular_index = snakemake.params.get("tag_per_molecular_index", "null")


if force_gc != "null":
    force_gc = "true" if force_gc else "false"
if apply_eamss_filter != "null":
    apply_eamss_filter = "true" if apply_eamss_filter else "false"
if include_non_pf_reads != "null":
    include_non_pf_reads = "true" if include_non_pf_reads else "false"
if ignore_unexpected_barcodes != "null":
    ignore_unexpected_barcodes = "true" if ignore_unexpected_barcodes else "false"

command = (
    "picard IlluminaBasecallsToSam"
    " {extra}"
    " BASECALLS_DIR={snakemake.input.basecalls_dir}"
    " BARCODES_DIR={barcodes_dir}"
    " LANE={lane}"
    " RUN_BARCODE={run_barcode}"
    " READ_GROUP_ID={read_group_id}"
    " SEQUENCING_CENTER={sequencing_center}"
    " RUN_START_DATE={run_start_date}"
    " PLATFORM={platform}"
    " READ_STRUCTURE={read_structure}"
    " LIBRARY_PARAMS={snakemake.input.library_params}"
    " FIVE_PRIME_ADAPTER={five_prime_adapter}"
    " THREE_PRIME_ADAPTER={three_prime_adapter}"
    " NUM_PROCESSORS={snakemake.threads}"
    " FIRST_TILE={first_tile}"
    " TILE_LIMIT={time_limit}"
    " FORCE_GC={force_gc}"
    " APPLY_EAMSS_FILTER={apply_eamss_filter}"
    " MAX_READS_IN_RAM_PER_TILE={max_reads_in_ram_per_tile}"
    " MINIMUM_QUALITY={minimum_quality}"
    " INCLUDE_NON_PF_READS={include_non_pf_reads}"
    " IGNORE_UNEXPECTED_BARCODES={ignore_unexpected_barcodes}"
    " MOLECULAR_INDEX_TAG={molecular_index_tag}"
    " MOLECULAR_INDEX_BASE_QUALITY_TAG={molecular_index_base_quality_tag}")

if adapters_to_check != "null":
    command += "".join(f" ADAPTERS_TO_CHECK={adapter}" for adapter in adapters_to_check)
else:
    command += " ADAPTERS_TO_CHECK={adapters_to_check}"

if tag_per_molecular_index != "null":
    command += "".join(f" TAG_PER_MOLECULAR_INDEX={tag}" for tag in tag_per_molecular_index)
else:
    command += " TAG_PER_MOLECULAR_INDEX={tag_per_molecular_index}"
command += " {log}"

shell(command)
