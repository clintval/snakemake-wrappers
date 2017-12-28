"""Snakemake wrapper for bwa mem ClipBam."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

rg_tag = snakemake.params.get("rg_tag", "PU")
output_dir = snakemake.params.get("output_dir", "null")
re_reverse = snakemake.params.get("re_reverse", True)
interleave = snakemake.params.get("interleave", False)
include_non_pf_reads = snakemake.params.get("include_non_pf_reads", False)
clipping_attribute = snakemake.params.get("clipping_attribute", "null")
clipping_action = snakemake.params.get("clipping_action", "null")
clipping_min_length = snakemake.params.get("clipping_min_length", 0)
read1_trim = snakemake.params.get("read1_trim", 0)
read1_max_bases_to_write = snakemake.params.get("read1_max_bases_to_write", "null")
read2_trim = snakemake.params.get("read2_trim", 0)
read2_max_bases_to_write = snakemake.params.get("read2_max_bases_to_write", "null")
quality = snakemake.params.get("quality", "null")
include_non_primary_alignments = snakemake.params.get("include_non_primary_alignments", False)


if re_reverse != "null":
    re_reverse = "true" if re_reverse else "false"
if interleave != "null":
    interleave = "true" if interleave else "false"
if include_non_pf_reads != "null":
    include_non_pf_reads = "true" if include_non_pf_reads else "false"
if include_non_primary_alignments != "null":
    include_non_primary_alignments = "true" if include_non_primary_alignments else "false"

program_record_id = snakemake.params.get("program_record_id", "null")
program_group_version = snakemake.params.get("program_group_version", "null")
program_group_command_line = snakemake.params.get("program_group_command_line", "null")
program_group_name = snakemake.params.get("program_group_name", "null")
paired_run = snakemake.params.get("paired_run", True)
clip_adapters = snakemake.params.get("clip_adapters", True)
is_bisulfite_sequence = snakemake.params.get("is_bisulfite_sequence", False)
aligned_reads_only = snakemake.params.get("aligned_reads_only", False)
max_insertions_or_deletions = snakemake.params.get("max_insertions_or_deletions", 1)
attributes_to_retain = snakemake.params.get("attributes_to_retain", "null")
attributes_to_remove = snakemake.params.get("attributes_to_remove", "null")
attributes_to_reverse = snakemake.params.get("attributes_to_reverse", "null")
attributes_to_reverse_complement = snakemake.params.get("attributes_to_reverse_complement", "null")
expected_orientations = snakemake.params.get("expected_orientations", "null")
aligner_proper_pair_flags = snakemake.params.get("aligner_proper_pair_flags", False)
sort_order = snakemake.params.get("sort_order", "coordinate")
primary_alignment_strategy = snakemake.params.get("primary_alignment_strategy", "BestMapq")
clip_overlapping_reads = snakemake.params.get("clip_overlapping_reads", True)
include_secondary_alignments = snakemake.params.get("include_secondary_alignments", True)
add_mate_cigar = snakemake.params.get("add_mate_cigar", True)
unmap_contaminant_reads = snakemake.params.get("unmap_contaminant_reads", False)
min_unclipped_bases = snakemake.params.get("min_unclipped_bases", 32)
matching_dictionary_tags = snakemake.params.get("matching_dictionary_tags", ["M5", "LN"])
unmapped_read_strategy = snakemake.params.get("unmapped_read_strategy", "DO_NOT_CHANGE")

if paired_run != "null":
    paired_run = "true" if paired_run else "false"
if clip_adapters != "null":
    clip_adapters = "true" if clip_adapters else "false"
if is_bisulfite_sequence != "null":
    is_bisulfite_sequence = "true" if is_bisulfite_sequence else "false"
if aligned_reads_only != "null":
    aligned_reads_only = "true" if aligned_reads_only else "false"
if aligner_proper_pair_flags != "null":
    aligner_proper_pair_flags = "true" if aligner_proper_pair_flags else "false"
if clip_overlapping_reads != "null":
    clip_overlapping_reads = "true" if clip_overlapping_reads else "false"
if include_secondary_alignments != "null":
    include_secondary_alignments = "true" if include_secondary_alignments else "false"
if add_mate_cigar != "null":
    add_mate_cigar = "true" if add_mate_cigar else "false"
if unmap_contaminant_reads != "null":
    unmap_contaminant_reads = "true" if unmap_contaminant_reads else "false"

if isinstance(attributes_to_retain, (list, tuple)):
    attributes_to_retain = "".join(
        f" ATTRIBUTES_TO_RETAIN={attribute_to_retain}"
        for attribute_to_retain in attributes_to_retain)
else:
    attributes_to_retain = f" ATTRIBUTES_TO_RETAIN={attributes_to_retain}"
if isinstance(attributes_to_remove, (list, tuple)):
    attributes_to_remove = "".join(
        f" ATTRIBUTES_TO_REMOVE={attribute_to_remove}"
        for attribute_to_remove in attributes_to_remove)
else:
    attributes_to_remove = f" ATTRIBUTES_TO_REMOVE={attributes_to_remove}"
if isinstance(attributes_to_reverse, (list, tuple)):
    attributes_to_reverse = "".join(
        f" ATTRIBUTES_TO_REVERSE={attribute_to_reverse}"
        for attribute_to_reverse in attributes_to_reverse)
else:
    attributes_to_reverse = f" ATTRIBUTES_TO_REVERSE={attributes_to_reverse}"
if isinstance(attributes_to_reverse_complement, (list, tuple)):
    attributes_to_reverse_complement = "".join(
        f" ATTRIBUTES_TO_REVERSE_COMPLEMENT={attribute_to_reverse_complement}"
        for attribute_to_reverse_complement in attributes_to_reverse_complement)
else:
    attributes_to_reverse_complement = f" ATTRIBUTES_TO_REVERSE_COMPLEMENT={attributes_to_reverse_complement}"
if isinstance(expected_orientations, (list, tuple)):
    expected_orientations = "".join(
        f" EXPECTED_ORIENTATIONS={expected_orientation}"
        for expected_orientation in expected_orientations)
else:
    expected_orientations = f" EXPECTED_ORIENTATIONS={expected_orientations}"
if isinstance(matching_dictionary_tags, (list, tuple)):
    matching_dictionary_tags = "".join(
        f" MATCHING_DICTIONARY_TAGS={matching_dictionary_tag}"
        for matching_dictionary_tag in matching_dictionary_tags)
else:
    matching_dictionary_tags = f" MATCHING_DICTIONARY_TAGS={matching_dictionary_tags}"

bwa_verbosity = snakemake.params.get("bwa_verbosity", "3")

command = (
    "picard SamToFastq"
    " {extra}"
    " INPUT={snakemake.input.unmapped_bam}"
    " FASTQ=/dev/stdout"
    " RG_TAG={rg_tag}"
    " OUTPUT_DIR={output_dir}"
    " RE_REVERSE={re_reverse}"
    " INTERLEAVE={interleave}"
    " INCLUDE_NON_PF_READS={include_non_pf_reads}"
    " CLIPPING_ATTRIBUTE={clipping_attribute}"
    " CLIPPING_ACTION={clipping_action}"
    " CLIPPING_MIN_LENGTH={clipping_min_length}"
    " READ1_TRIM={read1_trim}"
    " READ1_MAX_BASES_TO_WRITE={read1_max_bases_to_write}"
    " READ2_TRIM={read2_trim}"
    " READ2_MAX_BASES_TO_WRITE={read2_max_bases_to_write}"
    " QUALITY={quality}"
    " INCLUDE_NON_PRIMARY_ALIGNMENTS={include_non_primary_alignments}"
    " {log} | "
    "bwa mem"
    " -t {snakemake.threads}"
    " -p "
    " -v {bwa_verbosity}"
    " {snakemake.input.reference}"
    " /dev/stdin"
    " {log}"
    " | "
    "picard MergeBamAlignment "
    " {extra}"
    " UNMAPPED={snakemake.input.unmapped_bam}"
    " REFERENCE_SEQUENCE={snakemake.input.reference}"
    " ALIGNED=/dev/stdin"
    " OUTPUT={snakemake.output}"
    " PROGRAM_RECORD_ID={program_record_id}"
    " PROGRAM_GROUP_VERSION={program_group_version}"
    " PROGRAM_GROUP_COMMAND_LINE={program_group_command_line}"
    " PROGRAM_GROUP_NAME={program_group_name}"
    " PAIRED_RUN={paired_run}"
    " CLIP_ADAPTERS=false"
    " IS_BISULFITE_SEQUENCE={is_bisulfite_sequence}"
    " ALIGNED_READS_ONLY={aligned_reads_only}"
    " MAX_INSERTIONS_OR_DELETIONS={max_insertions_or_deletions}"
    "{attributes_to_retain}"
    "{attributes_to_remove}"
    "{attributes_to_reverse}"
    "{attributes_to_reverse_complement}"
    " READ1_TRIM={read1_trim}"
    " READ2_TRIM={read2_trim}"
    "{expected_orientations}"
    " ALIGNER_PROPER_PAIR_FLAGS=false"
    " SORT_ORDER=queryname"
    " PRIMARY_ALIGNMENT_STRATEGY={primary_alignment_strategy}"
    " CLIP_OVERLAPPING_READS={clip_overlapping_reads}"
    " INCLUDE_SECONDARY_ALIGNMENTS={include_secondary_alignments}"
    " ADD_MATE_CIGAR={add_mate_cigar}"
    " UNMAP_CONTAMINANT_READS={unmap_contaminant_reads}"
    " MIN_UNCLIPPED_BASES={min_unclipped_bases}"
    "{matching_dictionary_tags}"
    " UNMAPPED_READ_STRATEGY={unmapped_read_strategy}"
    " {log}")

shell(command)
