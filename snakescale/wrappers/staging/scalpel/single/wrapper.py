"""Snakemake wrapper for Scalpel."""

__author__ = 'fangyinlo'
__copyright__ = 'Copyright 2018, Fang Yin Lo'
__email__ = 'fangyinlot@gmail.com'
__license__ = 'MIT'

from pathlib import Path

from snakemake.shell import shell

kmer = snakemake.params.get('kmer', 25)
covthr = snakemake.params.get('covthr', 5)
lowcov = snakemake.params.get('lowcov', 2)
covratio = snakemake.params.get('covratio', 0.01)
radius = snakemake.params.get('radius', 100)
window = snakemake.params.get('window', 400)
maxregcov = snakemake.params.get('maxregcov', 10000)
step = snakemake.params.get('step', 100)
mapscore = snakemake.params.get('mapscore', 1)
pathlimit = snakemake.params.get('pathlimit', 1_000_000)
mismatches = snakemake.params.get('mismatches', 3)
outdir = snakemake.params.get('dir', './outdir')
numprocs = snakemake.params.get('numprocs', 1)
sample = snakemake.params.get('sample', 'ALL')
coords = snakemake.params.get('coords', 'null')
prefix = snakemake.input.get('scalpel_path')

if prefix is None:
    scalpel_path = 'scalpel-discovery'

else:
    prefix = Path(prefix)
    scalpel_path = str(prefix / 'scalpel-discovery')

scalpel_command = (
    "{scalpel_path}"
    " --single"
    " --bam {snakemake.input.bam}"
    " --bed {snakemake.input.bed}"
    " --ref {snakemake.input.ref}"
    " --kmer {kmer}"  # < int> k-mer size [default 25]
    " --covthr {covthr}"  # < int> threshold used to select source and sink [default 5]
    " --lowcov  {lowcov}"  # < int> threshold used to remove low-coverage nodes [default 2]
    " --covratio {covratio}"  # < float> minimum coverage ratio for sequencing errors (default: 0.01)
    " --radius {radius}"  # < int> left and right extension (in base-pairs) [default 100]
    " --window {window}"  # < int> window-size of the region to assemble (in base-pairs) [default 400]
    " --maxregcov {maxregcov}"  # < int> maximum average coverage allowed per region [default 10000]
    " --step {step}"  # < int> delta shift for the sliding window (in base-pairs) [default 100]
    " --mapscore {mapscore}"  # < int> minimum mapping quality for selecting reads to assemble [default 1]
    " --pathlimit {pathlimit}"  # < int> limit number of sequence paths to [default 1000000]
    " --mismatches {mismatches}"  # < int> max number of mismatches in near-perfect repeat detection [default 3]
    " --dir {outdir}"  # < directory> output directory [default ./outdir]
    " --numprocs {numprocs}"  # < int> number of parallel jobs (1 for no parallelization) [default 1]
    " --sample {sample}"  # < string> only process reads/fragments in sample [default ALL]
    " --coords {coords}"  # < file> file with list of selected locations to examine [defau
)

shell(scalpel_command)
