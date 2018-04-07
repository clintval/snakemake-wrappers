from setuptools import setup

setup(
    name='snakemake_wrappers',
    packages=['snakemake_wrappers'],
    version='0.1.0',
    description='A collection of awesome snakemake wrappers.',
    long_description=long_description,
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url='https://github.com/clintval/snakemake-wrappers',
    download_url='https://github.com/clintval/snakemake-wrappers/archive/v0.1.0.tar.gz',
    py_modules=['snakemake_wrappers'],
    install_requires=[],
    extras_require={
        'ci': ['nose', 'codecov'],
        'fancytest': ['nose', 'nose-progressive', 'coverage'],
    },
    license='MIT',
    zip_safe=True,
    keywords=[
        'signature',
        'mutation',
        'transition'
        'transversion',
        'spectra',
        'bioinformatics'
    ],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
    ]
)
