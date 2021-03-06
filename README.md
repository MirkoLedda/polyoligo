# `PolyOligo`

Design oligos for complex genomes. (alpha version)

## What Is It?

A Python-based **multicore enabled** set of tools to **design oligonucleotides**. Supports **SNPs** and **indels** and can use known **mutations** across populations to refine designs.

A webapp version of this software is available [HERE](http://ec2-52-52-41-39.us-west-1.compute.amazonaws.com/) !

## Included binaries (click on the command for more information)

* [**`polyoligo-pcr`**](https://github.com/MirkoLedda/polyoligo/blob/master/README_pcr.md): Standard PCR or Sanger sequencing primers.
* [**`polyoligo-kasp`**](https://github.com/MirkoLedda/polyoligo/blob/master/README_kasp.md): Genotyping by Kompetitive allele specific PCR (KASP).
* [**`polyoligo-caps`**](https://github.com/MirkoLedda/polyoligo/blob/master/README_caps.md): Genotyping by Cleaved Amplified Polymorphic Sequences (CAPS).
* [**`polyoligo-hrm`**](https://github.com/MirkoLedda/polyoligo/blob/master/README_hrm.md): Genotyping by high-resolution melting (HRM) PCR.
* [**`polyoligo-crispr`**](https://github.com/MirkoLedda/polyoligo/blob/master/README_crispr.md): Guide RNAs for CRISPR/Cas9 assays.

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

## Running a test

To make sure the installation ran properly, run:

```
polyoligo
```

## Citation

Ledda M., Cobo N., Lorant A., Hardigan M.A. and Knapp S.J., PolyOligo: A Bioinformatic Platform for Identifying Target DNA Sequences for the Development of Sub-Genome Specific DNA Markers in Polyploid/Complex Genomes. Poster presented at: Annual Conference of the American Society of Horticultural Sciences; 2019 July 21-25; Las Vegas, NV, USA

## Reporting bugs and requesting features

This software is actively supported and all changes are listed in the [CHANGELOG](CHANGES.md). To report a bug open a ticket in the [issues tracker](https://github.com/MirkoLedda/polyoligo/issues). Features can be requested by opening a ticket in the [pull request](https://github.com/MirkoLedda/polyoligo/pulls).

## Contributors

* [**Mirko Ledda**](https://mirkoledda.github.io/) - *Initial implementation and developer*
* [**Nicolas Cobo**](https://github.com/ncobo) - *Requirements and biological accuracy*

## License

This software is licensed under the BSD-2 License - see the [LICENSE](https://github.com/MirkoLedda/polyoligo/blob/master/LICENSE) file for details.
