# snakescale

[![Testing Status](https://travis-ci.org/clintval/snakescale.svg?branch=master)](https://travis-ci.org/clintval/snakescale)
[![codecov](https://codecov.io/gh/clintval/snakescale/branch/master/graph/badge.svg)](https://codecov.io/gh/clintval/snakescale)
[![Documentation Build Status](https://readthedocs.org/projects/snakescale/badge/?version=latest)](https://snakescale.readthedocs.io/en/latest/?badge=latest)
[![PyPi Release](https://badge.fury.io/py/snakescale.svg)](https://badge.fury.io/py/snakescale)
[![Python Versions](https://img.shields.io/pypi/pyversions/snakescale.svg)](https://pypi.python.org/pypi/snakescale/)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Non-strict wrappers for the data pipelining language Snakemake.

```bash
❯ pip install snakescale
```

Features:

- Do the wrappers in the [official wrapper repository](https://bitbucket.org/snakemake/snakemake-wrappers) get you half of the way to writing rules in only Python syntax?
- Do want your rules fully parameterized with the `input`, `output`, `resources`, and `params` keys only?
- Do you want to use the builtin Python types as values to a rule?
- Do you want to use the Snakemake resource system for JVM resources?
- Do you want a Snakemake wrapper which hard-codes as little as possible besides the **style** of the CLI it's wrapping?
Read the documentation at: [snakescale.readthedocs.io](http://snakescale.readthedocs.io/)

This project aims to wrap bioinformatics utilities with style and variable type converters instead of strict, inflexible shell templates. The wrappers in this project are unaware of the command line flags of the tool the wrapper is wrapping!

### Quickstart

For example, [`Picard`](https://broadinstitute.github.io/picard/) tools accept a list in the form of repeating options:

```bash
❯ picard_tool \
    INPUT=infile \
    OUTPUT=outfile \
    INCLUDE_NON_PF_READS=false \
    ADAPTERS_TO_CHECK=option1 \
    ADAPTERS_TO_CHECK=option2
```

With these arguments as their respective Python types, this is how you would construct that Snakemake rule without a wrapper:

```python
include_non_pf_reads = False
adapters_to_check = ['option1', 'option2']

rule picard_tool:
    input: infile
    output: outfile
    params:
        include_non_pf_reads_option='true' if include_non_pf_reads else 'false'
        adapters_to_check_options=' '.join([f'ADAPTERS_TO_CHECK={s}' for s in adapters_to_check])
    shell: (
      'picard'
      ' INPUT={input} OUTPUT={output}'
      ' INCLUDE_NON_PF_READS={params.include_non_pf_reads_option}'
      ' {params.adapters_to_check_options}'
    )
```

Using the wrappers in this library, the rule takes the form:

```python
from snakescale import wrappers as run

rule picard_tool:
    input: infile
    output: outfile
    params:
        adapters_to_check=['option1', 'option2']
        include_non_pf_reads=False
    wrapper: run.picard_tool
```

#### Complete Example

```python
from snakescale import task

lanes = 4

rule illumina_basecalls_to_sam:
    """Demultiplexes an Illumina run folder of basecalls."""
    input:
        basecalls_dir=basecalls_dir,
        barcodes_dir=rules.extract_illumina_barcodes.output.barcodes_dir,
        library_params=rules.create_basecalling_params.output.library_params
    output: output_files
    threads: max(1, int(threads / len(lanes)))
    resources:
        gc_time_limit=50,
        gc_heap_free_limit=10,
        heap_size=lambda wildcards, attempt: int(attempt * 1.5 * 8192),
        samjdk_buffer_size=131072,
        use_async_io_read_samtools=1,
        use_async_io_write_samtools=1
    params:
        lane='{lane}',
        run_barcode=run_barcode,
        sequencing_center=sequencing_center,
        run_start_date=run_start_date,
        read_structure=read_structure,
        adapters_to_check=['INDEXED', 'NEXTERA_V2', 'FLUIDIGM'],
        include_non_pf_reads=False
    log: run_output / 'logs' / '{rule.name}.{lane}.log'
    wrapper: task('0.1.0', 'picard', 'illumina_basecalls_to_sam')
```

The above rule will execute the following shell command.

```bash
❯ picard \
      -XX:GCHeapFreeLimit=50 \
      -XX:GCTimeLimit=10 \
      -Xmx12288m \
      -Dsamjdk.buffer_size=131072 \
      -Dsamjdk.use_async_io_read_samtools=true \
      -Dsamjdk.use_async_io_write_samtools=true \
    IlluminaBasecallsToSam \
      BASECALLS_DIR="${basecalls_dir}" \
      BARCODES_DIR="${barcodes_dir}" \
      LIBRARY_PARAMS="${library_params}" \
      OUTPUT="${output_file}"
      LANE="${lane}" \
      RUN_BARCODE="${run_barcode}" \
      SEQUENCING_CENTER="${sequencing_center}" \
      RUN_START_DATE="${run_start_date}" \
      READ_STRUCTURE="${read_structure}" \
      ADAPTERS_TO_CHECK=INDEXED \
      ADAPTERS_TO_CHECK=NEXTERA_V2 \
      ADAPTERS_TO_CHECK=FLUIDIGM \
      INCLUDE_NON_PF_READS=false \
    > ".../logs/illumina_basecalls_to_sam.${lane}.log" 2>&1
```

#### Specifying Resources and Resource Constraints

Any toolkit with JVM resources also supports assignment of these key values to the [`resources`](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources) dictionary of the rule. This allows for CLI overriding of these resources! For instance, to turn off `htsjdk` asynchronous IO you can pass:

```bash
❯ snakemake --resources use_async_io_read_samtools 0 use_async_io_write_samtools 0
```
