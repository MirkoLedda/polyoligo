# `polyoligo`

Design oligos for complex genomes.

## What Is It?

A cross-platform Python-based set of tools to design oligonucleotides on complex genomes.

Handles polyploid organims and genomic duplications, as well as known mutations within a selected population.

Threads are fully parallelized to benefit from multi-core system. 

## What Is Included?

* **`polyoligo-kasp`**: Design optimal primers for genotyping assays using Kompetitive allele specific PCR (KASP).
<!--* **`polyoligo-pcr`**: Design optimal primers for PCR.-->
<!--* **`polyoligo-snpseq`**: Design optimal primers for SNP-Seq assays.-->
* **`polyoligo-crispr`**: Design gRNAs for CRISPR/Cas9 assays.

## Getting Started

The software has a command-line interface which is accessed using your favorite terminal. Note that all positions use a **1-based indexing**.


### Prerequisites

**OS**: MacOSX or Linux

**Python v3 or newer**

To check if Python3 is installed on your system run the following command:

```
python3 -V
```

If the command failed install Python3 for your OS from the [official Python website](https://www.python.org/downloads/).

It is also recommended to use the lastest version of pip and setuptools, which can be installed by running the command:

```
python3 -m pip install -U pip setuptools --user
```

### Installation

Installation is done directly from source. For that, clone this repository using the commands:

```
git clone https://github.com/MirkoLedda/polyoligo.git
cd polyoligo
```

If you do have `sudo` permissions, such as on a personal laptop, then use the command:

```sudo python3 setup.py install```

Otherwise, run a local installation using the commands:

```
python3 setup.py install --user
echo 'export PATH="$PATH:~/.local/bin"' >> ~/.bashrc; source ~/.bashrc
```

### Instructions

<!--See this [wiki](https://github.com/MirkoLedda/polyoligo.git).-->

* [**`polyoligo-kasp`**](https://github.com/MirkoLedda/polyoligo/README_kasp.md)
* [**`polyoligo-crispr`**](https://github.com/MirkoLedda/polyoligo/README_crispr.md)
<!--* **`polyoligo-pcr`**: Design optimal primers for PCR.-->
<!--* **`polyoligo-snpseq`**: Design optimal primers for SNP-Seq assays.-->



## Citation

TBA

## Reporting bugs and requesting features

**`polyoligo`** is actively supported and all changes are listed in the [CHANGELOG](CHANGES.md). To report a bug open a ticket in the [issues tracker](https://github.com/MirkoLedda/polyoligo/issues). Features can be requested by opening a ticket in the [pull request](https://github.com/MirkoLedda/polyoligo/pulls).

<!-- ## Versioning

We use the [SemVer](http://semver.org/) convention for versioning. For the versions available, see the [tags on this repository](TBA/tags). -->

## Contributors

* [**Mirko Ledda**](https://mirkoledda.github.io/) - *Initial implementation and developer*
* [**Nicolas Cobo**](https://github.com/ncobo) - *Requirements and biological accuracy*

## License

This software is licensed under the BSD-2 License - see the [LICENSE](LICENSE.txt) file for details.
