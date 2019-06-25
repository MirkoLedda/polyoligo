# `polyoligo-crispr`

### Running a test

To make sure the installation completed, run the command:

```
polyoligo-crispr -v
```

You can now run a full test by entering the command:

```
polyoligo-crispr Fvb2-4 1000 2000 out sample_data/blastdb
```

> This will create three output files called `out.txt`, `out.bed` and `out.log` in the current folder.


## Usage

### General usage and available options

```
polyoligo-crispr <CHROM/CONTIG> <START> <END> <OUTPUT> <FASTA/BLASTDB> <OPTIONS>
```

For a list of all available options and their descriptions, type:

```
polyoligo-crispr -h
```

Recommendations (when applicable) are given in the option caption. Note that switches, i.e. boolean options that do not need arguments, have defaults set to `False`.

### Inputs
The software requires three mandatory inputs:

**`CHROM/CONTIG`**: Target chromosome or contig.

**`START`**: Target region start.

**`END`**: Target region end.

**`<OUTPUT>`**: The base name of the output files.

**`<FASTA/BLASTDB>`**: Either a FASTA file or a BLAST database to use as the reference genome.


### Outputs
Three output files are produced:

**`<OUTPUT>.log`**: Log file.

**`<OUTPUT>.bed`**: BED file containing the list of designed gRNA.

**`<OUTPUT>.txt`** List of gRNAs with the following columns:

HR START END STRAND SEQ PAM N_PERFECT_OFFTARGETS N_OFFTARGETS_12MER N_OFFTARGETS_8MER 12_MER 8_MER TM TTTT

|Column|Description|
|---|---|
|`CHR`|Chromosome|
|`START`|Start|
|`END`|END|
|`STRAND`|Strand where the gRNA will bind|
|`SEQ`|gRNA sequence|
|`PAM`|PAM site|
|`N_PERFECT_OFFTARGETS`|Number of offtargets with perfect match to the gRNA sequence|
|`N_OFFTARGETS_12MER`|Number of offtargets with perfect match to 12-MER seed region of the gRNA and no more than 2 mutations elsewhere|
|`N_OFFTARGETS_8MER`|Number of offtargets with perfect match to 8-MER seed region of the gRNA and no more than 2 mutations elsewhere|
|`12_MER`|Number of offtargets with perfect match to 12-MER seed region of the gRNA|
|`8_MER`|Number of offtargets with perfect match to 8-MER seed region of the gRNA|
|`TM`|Melting temperature of the gRNA|
|`TTTT`|Does the gRNA contain TTTT repeats. 0=No or 1=Yes.|
