import setuptools
import shutil
from setuptools import find_packages

from pathlib import Path

PACKAGE = 'snakescale'
VERSION = '0.8.0'

for path in Path('snakescale/wrappers/').rglob('.snakemake'):
    shutil.rmtree(path)

setuptools.setup(
    name=PACKAGE,
    version=VERSION,
    author='clintval',
    author_email='valentine.clint@gmail.com',
    description='Non-strict wrappers for the data pipelining language Snakemake.',
    url=f'https://github.com/clintval/{PACKAGE}',
    download_url=f'https://github.com/clintval/{PACKAGE}/archive/v{VERSION}.tar.gz',
    long_description=Path('README.md').read_text(),
    long_description_content_type='text/markdown',
    license='MIT',
    zip_safe=False,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    extras_require={'full': 'snakemake>=4.5.1'},
    keywords='snakemake pipeline wrapper',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    project_urls={
        'Documentation': f'https://{PACKAGE}.readthedocs.io',
        'Issue Tracker': f'https://github.com/clintval/{PACKAGE}/issues',
    },
)
