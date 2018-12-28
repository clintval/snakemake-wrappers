"""Snakemake wrapper for VarDictJava."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from pathlib import Path

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

do_not_print_end_tag = snakemake.params.get('do_not_print_end_tag', True)
col_chromosome = snakemake.params.get('col_chromosome', 1)
col_end = snakemake.params.get('col_end', 3)
col_gene_name = snakemake.params.get('col_gene_name', 4)
col_start = snakemake.params.get('col_start', 2)
coordinates_are_zero_based = snakemake.params.get('coordinates_are_zero_based', True)
do_not_print_end_tag = snakemake.params.get('do_not_print_end_tag', True)
do_pileup = snakemake.params.get('do_pileup', True)
hex_filter = snakemake.params.get('hex_filter', '0x500')
include_Ns = snakemake.params.get('include_Ns', False)
max_mean_mismatches = snakemake.params.get('max_mean_mismatches', 5.25)
max_mismatches_per_read = snakemake.params.get('max_mismatches_per_read', 8)
min_allele_frequency = snakemake.params.get('min_allele_frequency', 0.05)
min_allele_frequency_for_homozygous = snakemake.params.get(
    'min_allele_frequency_for_homozygous', 0.25
)
min_base_quality_score = snakemake.params.get('min_base_quality_score', 20)
min_high_quality_variant_depth = snakemake.params.get('min_high_quality_variant_depth', 2)
min_mapping_quality = snakemake.params.get('min_mapping_quality', 10.0)
min_mean_base_quality = snakemake.params.get('min_mean_base_quality', 22.5)
min_total_depth = snakemake.params.get('min_total_depth', 3)
min_variant_reads = snakemake.params.get('min_variant_reads', 2)
minimum_mean_position = snakemake.params.get('minimum_mean_position', 8)
output_gvcf = snakemake.params.get('output_gvcf', False)
pstd_filter = snakemake.params.get('pstd_filter', 1)
sample_name = snakemake.wildcards.sample

prefix = snakemake.params.get('VarDictJavaRoot', None)

if prefix is None:
    vardict_bin = 'VarDict'
    teststrandbias_bin = 'teststrandbias.R'
    var2vcf_valid_bin = 'var2vcf_valid.pl'
else:
    prefix = Path(prefix)
    vardict_bin = str(prefix / 'build/install/VarDict/bin/VarDict')
    teststrandbias_bin = str(prefix / 'VarDict/teststrandbias.R')
    var2vcf_valid_bin = str(prefix / 'VarDict/var2vcf_valid.pl')

vardict_command = (
    "{vardict_bin}"
    " -b {snakemake.input.bam}"  # The indexed BAM file
    " -G {snakemake.input.reference}"  # The indexed reference FASTA
    " -N {sample_name}"  # Sample name to be used directly
    " -c {col_chromosome}"  # The column for chromosome
    " -S {col_start}"  # The column for region start, e.g. gene start
    " -E {col_end}"  # The column for region end, e.g. gene end
    " -g {col_gene_name}"  # The column for gene name, or segment annotation
    " -F {hex_filter}"  # The hexical to filter reads using samtools
    " -f {min_allele_frequency}"  # The threshold for allele frequency
    " -r {min_variant_reads}"  # The minimum number of variant reads
    " -q {min_base_quality_score} "  # The phred score for a base to be considered a good call
    " -m {max_mismatches_per_read}"  # Reads with mismatches greater than INT will be ignored
    " -th {snakemake.threads}"  # Threads count
)
vardict_command += " -p" if do_pileup else ""  # Do pileup regardless of frequency
vardict_command += " -K" if include_Ns else ""  # Include Ns in the total depth calculation
vardict_command += (
    " -z 1" if coordinates_are_zero_based else " -z 0"
)  # Whether coordinates are 0 or 1-based
vardict_command += " {snakemake.input.bed}"  # "region_info" file

filter_for_variants_only = (
    r'awk "{{if (\$6 != \$7) print}}"'  # Filter for non-reference matching alternate records
)

test_strand_bias_command = (
    "{teststrandbias_bin}"  # Performs fisher exact test on variants called on either strand
)

var2vcf_valid_command = (
    "{var2vcf_valid_bin}"  # Convert .variant to .vcf file
    " -m {max_mean_mismatches}"  # The maximum mean mismatches allowed, Mismatches do not include indels in the alignment
    " -p {minimum_mean_position}"  # The minimum mean position of variants in the read
    " -P {pstd_filter}"  # Whether to filter variants with pstd = 0
    " -q {min_mean_base_quality}"  # The minimum mean base quality
    " -Q {min_mapping_quality}"  # The minimum mapping quality
    " -d {min_total_depth}"  # The minimum total depth
    " -v {min_high_quality_variant_depth}"  # The minimum high quality variant depth
    " -f {min_allele_frequency}"  # The minimum allele frequency
    " -F {min_allele_frequency_for_homozygous}"  # The minimum allele frequency to consider to be homozygous
    " -N {sample_name}"  # The sample name to be used directly
)
var2vcf_valid_command += " -E" if do_not_print_end_tag else ""  # If set, do not print END tag

if output_gvcf:
    command = " | ".join([vardict_command, test_strand_bias_command, var2vcf_valid_command])
else:
    command = " | ".join(
        [
            vardict_command,
            filter_for_variants_only,
            test_strand_bias_command,
            var2vcf_valid_command,
        ]
    )

command += " > {snakemake.output} {log}"  # Pipe to output VCF and log stderr

shell(command)
