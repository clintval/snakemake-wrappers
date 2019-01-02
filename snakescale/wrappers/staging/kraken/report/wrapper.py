"""Snakemake wrapper for Kraken report."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.utils import make_kraken_params

params = make_kraken_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get('extra', '')

shell('kraken-report' ' {extra}' ' {params}' ' {snakemake.input}' ' > {snakemake.output}' ' {log}')
