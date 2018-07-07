<h1 align="center">snakemake-wrappers</h2>

<p align="center">A Collection of Awesome Snakemake Wrappers</p>

<p align="center">
  <a href="#installation"><strong>Installation</strong></a>
  ·
  <a href="#wrapper-documentation"><strong>Wrapper Documentation</strong></a>
</p>

<p align="center">
  <a href="https://travis-ci.org/clintval/snakemake-wrappers"><img src="https://travis-ci.org/clintval/snakemake-wrappers.svg?branch=master"></img></a>
  <a href="https://codecov.io/gh/clintval/snakemake-wrappers"><img src="https://codecov.io/gh/clintval/snakemake-wrappers/branch/master/graph/badge.svg"></img></a>
  <a href="https://badge.fury.io/py/snakemake-wrappers"><img src="https://badge.fury.io/py/snakemake-wrappers.svg" alt="PyPI version"></img></a>
  <a href="https://codeclimate.com/github/clintval/snakemake-wrappers/maintainability"><img src="https://api.codeclimate.com/v1/badges/7f6ce7780716a92c40b8/maintainability"></img></a>
  <a href="https://github.com/clintval/snakemake-wrappers/blob/master/LICENSE"><img src="https://img.shields.io/pypi/l/snakemake-wrappers.svg"></img></a>
</p>

<br>

- Do the wrappers in the [official wrapper repository](https://bitbucket.org/snakemake/snakemake-wrappers) get you half of the way to writing rules in only Python syntax?
- Do want your rules fully parameterized with the `input`, `output`, `resources`, and `params` keys only?
- Do you want to use the builtin Python types as values to a rule?
- Do you want to use the Snakemake resource system for JVM resources?
- Do you want a Snakemake wrapper which hard-codes as little as possible besides the **style** of the CLI it's wrapping?

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

include_non_pf_reads_option = 'true' if include_non_pf_reads else 'false'
adapters_to_check_options = ' '.join([f'ADAPTERS_TO_CHECK={s}' for s in adapters_to_check])

rule picard_tool:
    input: infile
    output: outfile
    shell: (
      'picard'
      ' INPUT={input} OUTPUT={output}'
      ' INCLUDE_NON_PF_READS={include_non_pf_reads_option}'
      ' {adapters_to_check_options}'
    )
```

Using the wrappers in this library, the rule takes the form:

```python
from snakemake_wrappers import wrappers as run

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
from snakemake_wrappers import wrappers as run

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
    log: f'{run_output}/logs/{{rule.name}}.{{lane}}.log'
    wrapper: run.picard_illumina_basecalls_to_sam
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

<br>

<h3 align="center">Installation</h3>

```
# Do not install yet, as we have yet to release a functional prototype!
❯ # pip install git+https://github.com/clintval/snakemake-wrappers
```

<br>

<h3 align="center">Wrapper Documentation</h3>

Every tool's complete example is provided here:

TODO

<br>

