"""Snakemake wrapper for dwgsim."""

__author__ = 'clintval'
__copyright__ = 'Copyright 2018, Clint Valentine'
__email__ = 'valentine.clint@gmail.com'
__license__ = 'MIT'

from snakemake.shell import shell

from snakescale.utils import make_dwgsim_params

output_prefix = snakemake.params["output_prefix"]

extra = snakemake.params.get('extra', '')
params = make_dwgsim_params(snakemake.params)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell('dwgsim' ' {extra}' ' {params}' ' {snakemake.input}' ' {output_prefix}' ' {log}')
