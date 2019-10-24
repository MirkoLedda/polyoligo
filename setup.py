"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from os.path import join
from setuptools import find_packages
from setuptools import setup

from src.polyoligo._version import __version__

# Get the long description from the README file
with open(join("README.md")) as f:
    long_description = f.read()

setup(
    name="polyoligo",
    version=__version__,
    description="Set of tools to design oligos with complex genomes.",
    long_description=long_description,
    url="https://github.com/MirkoLedda/polyoligo",
    author="Mirko Ledda",
    author_email="maledda@ucdavis.edu",
    license="BSD-2",
    keywords=[
        "genomics",
        "bioinformatics",
        "polyploid",
        "genotyping",
        "KASP",
        "CAPS",
        "PCR",
        "CRISPR",
    ],

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
    ],

    # Python requirements
    python_requires=">=3, <4",

    # Specify source packages
    packages=find_packages('src'),
    package_dir={'': 'src'},

    # Dependencies
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "pysam",
        "PyVCF",
        "PyYAML",
        "tqdm",
        "primer3-py",
    ],

    # Create the executable
    entry_points={
        'console_scripts': [
            'polyoligo = polyoligo.logo:main',
            'polyoligo-kasp = polyoligo.cli_kasp:main',
            'polyoligo-crispr = polyoligo.cli_crispr:main',
            'polyoligo-pcr = polyoligo.cli_pcr:main',
            'polyoligo-caps = polyoligo.cli_caps:main',
            'polyoligo-hrm = polyoligo.cli_hrm:main',
        ]
    },

    # Requirements for tests and coverage analysis
    setup_requires=["pytest-runner"],
    tests_require=[
        "pytest",
        "coverage",
        "pytest_cov",
        "memory_profiler",
        "matplotlib",
    ],

    package_data={
        'bin/linux_x64': ['bin/linux_x64/*'],
        'bin/macosx_x64': ['bin/macosx_x64/*'],
        'bin/win_x64': ['bin/win_x64/*'],
        'data': ['data/*'],
    },
    include_package_data=True,
    zip_safe=False,
)
