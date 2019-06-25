"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from os.path import join
from setuptools import find_packages
from setuptools import setup
# import os
# from setuptools.command.install_scripts import install_scripts
# from setuptools.command.install import install
# from distutils import log  # needed for outputting information messages

from src.polyoligo._version import __version__

# Get the long description from the README file
with open(join("README.md")) as f:
    long_description = f.read()

# class OverrideInstall(install):
#
#     def run(self):
#         uid, gid = 0, 0
#         mode = 0o771
#         install.run(
#             self)  # ensures the install proceed as usual
#         # Update the file permissions to make everything executable
#         for filepath in self.get_outputs():
#             if self.install_scripts in filepath:
#                 log.info("Overriding setuptools mode of scripts ...")
#                 log.info("Changing ownership of {} to uid:{} gid {}".format(filepath, uid, gid))
#                 os.chown(filepath, uid, gid)
#                 log.info("Changing permissions of {} to {}".format(filepath, oct(mode)))
#                 os.chmod(filepath, mode)


setup(
    name="polyoligo",
    version=__version__,
    description="Set of tools to design oligos with complex genomes.",
    long_description=long_description,
    url="",  # todo
    author="Mirko Ledda",
    author_email="maledda@ucdavis.edu",
    license="BSD-2",
    keywords=["genomics", "bioinformatics", "complex genomes", "KASP", "CRISPR"],

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
    ],

    # Python requirements
    python_requires=">=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4",

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
            'polyoligo-kasp = polyoligo.cli_kasp:main',
            'polyoligo-crispr = polyoligo.cli_crispr:main',
            'polyoligo-pcr = polyoligo.cli_pcr:main',
        ]
    },

    # Requirements for tests and coverage analysis
    setup_requires=["pytest-runner"],
    tests_require=["pytest",
                   "coverage",
                   "pytest_cov"],

    package_data={
        'bin/linux_x64': ['bin/linux_x64/*'],
        'bin/macosx_x64': ['bin/macosx_x64/*'],
        'bin/win_x64': ['bin/win_x64/*'],
    },
    include_package_data=True,
    zip_safe=False,
    # Run the custom install
    # cmdclass={'install': OverrideInstall},
)
