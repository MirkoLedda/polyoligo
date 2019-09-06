# `polyoligo-crispr`

## Usage

### General usage and available options

```
polyoligo-crispr <CHROM:START-END> <OUTPUT> <FASTA/BLASTDB> <OPTIONS>
```

For a list of all available options and their descriptions, type:

```
polyoligo-crispr -h
```

> Recommendations (when applicable) are given in the option caption. Note that switches, i.e. boolean options that do not need arguments, have defaults set to `False`.

### Example usage
In the following example, CRISPR/Cas9 gRNAs will be designed within the region of interest (ROI):

```
polyoligo-crispr Fvb2-4:4045000-4045200 out sample_data/blastdb
```

### Inputs
The software requires three mandatory inputs:

**`ROI`**: Target region of interest defined as CHR:START-END.

**`<OUTPUT>`**: The base name of the output files.

**`<FASTA/BLASTDB>`**: A FASTA file and/or a BLAST database to use as the reference genome. Both file types can be provided by using the same basename. If either is provided, then a conversion will automatically be made to obtain both file types.


### Outputs
Three output files are produced:

**`<OUTPUT>.log`**: Log file.

**`<OUTPUT>.bed`**: BED file containing the list of designed gRNA. Names are composites of `<N_PERFECT_OFFTARGETS>_<N_OFFTARGETS_12MER>_<N_OFFTARGETS_8MER>_<12_MER>_<8_MER>_<TM>_<TTTT>` (see below).

**`<OUTPUT>.txt`** List of gRNAs with the following columns:

|Column|Description|
|---|---|
|`CHR`|Chromosome|
|`START`|Start|
|`END`|END|
|`STRAND`|Strand where the gRNA will bind|
|`SEQ`|gRNA sequence|
|`PAM`|PAM site|
|`N_PERFECT_OFFTARGETS`|Number of off-targets with perfect match to the gRNA sequence|
|`N_OFFTARGETS_12MER`|Number of off-targets with perfect match to 12-MER seed region of the gRNA and no more than 2 mutations elsewhere|
|`N_OFFTARGETS_8MER`|Number of off-targets with perfect match to 8-MER seed region of the gRNA and no more than 2 mutations elsewhere|
|`12_MER`|Number of off-targets with perfect match to 12-MER seed region of the gRNA|
|`8_MER`|Number of off-targets with perfect match to 8-MER seed region of the gRNA|
|`TM`|Melting temperature of the gRNA|
|`TTTT`|Does the gRNA contain TTTT repeats. 0=No or 1=Yes.|
