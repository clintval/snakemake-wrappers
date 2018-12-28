# snakescale

[![Testing Status](https://travis-ci.org/clintval/snakescale.svg?branch=master)](https://travis-ci.org/clintval/snakescale)
[![codecov](https://codecov.io/gh/clintval/snakescale/branch/master/graph/badge.svg)](https://codecov.io/gh/clintval/snakescale)
[![Documentation Build Status](https://readthedocs.org/projects/snakescale/badge/?version=latest)](https://snakescale.readthedocs.io/en/latest/?badge=latest)
[![PyPi Release](https://badge.fury.io/py/snakescale.svg)](https://badge.fury.io/py/snakescale)
[![Anaconda-Server Badge](https://anaconda.org/clintval/snakescale/badges/version.svg)](https://anaconda.org/clintval/snakescale)
[![Python Versions](https://img.shields.io/pypi/pyversions/snakescale.svg)](https://pypi.python.org/pypi/snakescale/)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Non-strict wrappers for the data pipelining language Snakemake.

```bash
❯ pip install 'snakescale[full]'
```

Features:

- Do the wrappers in the [official wrapper repository](https://bitbucket.org/snakemake/snakemake-wrappers) get you half of the way to writing rules in only Python syntax?
- Do want your rules fully parameterized with the `input`, `output`, `resources`, and `params` keys only?
- Do you want to use the builtin Python types as values to a rule?
- Do you want to use the Snakemake resource system for JVM resources?
- Do you want a Snakemake wrapper which hard-codes as little as possible besides the **style** of the CLI it's wrapping?
Read the documentation at: [snakescale.readthedocs.io](http://snakescale.readthedocs.io/)

This project aims to wrap bioinformatics utilities with style and variable converters instead of strict, inflexible shell templates.

The wrappers in this project are unaware of the command line flags of the tool the wrapper is wrapping!
