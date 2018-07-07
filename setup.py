from setuptools import setup

NAME = 'snakemake-wrappers'
PACKAGE = 'snakemake_wrappers'
VERSION = '0.1.0'

ARTIFACT = f'https://github.com/clintval/{NAME}/archive/v{VERSION}.tar.gz'
URL = f'https://github.com/clintval/{NAME}'

setup(
    name=PACKAGE,
    packages=[PACKAGE],
    version=VERSION,
    description='A collection of awesome snakemake wrappers',
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url=URL,
    download_url=ARTIFACT,
    py_modules=[PACKAGE],
    license='MIT',
    zip_safe=True,
    keywords='snakemake',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
    ]
)
