# `polyoligo-pcr`

### Running a test

To make sure the installation completed, run the command:

```
polyoligo-pcr -v
```

You can now run a full test by entering the command:

```
polyoligo-pcr sample_data/pcr_targets.txt out sample_data/blastdb
```

> This will create two output files called `out.txt` and `out.log` in the current folder.


## Usage

### General usage and available options

```
polyoligo-pcr <ROI> <OUTPUT> <BLASTDB> <OPTIONS>
```

For a list of all available options and their descriptions, type:

```
polyoligo-pcr -h
```

Recommendations (when applicable) are given in the option caption. Note that switches, i.e. boolean options that do not need arguments, have defaults set to `False`.

### Inputs
The software requires three mandatory inputs:

**`<INPUT>`**: Regions of interest declared as CHR:START-END NAME(optional).

**`<OUTPUT>`**: The base name of the output files.

**`<FASTA/BLASTDB>`**: Either a FASTA file or a BLAST database to use as the reference genome.

Optional files include:

**`--vcf`**: A VCF file to design primers considering mutations (both SNPs and indels).

> Note that a tabix index file, created using the Samtools' [`tabix`](http://www.htslib.org/doc/tabix.html) binary is required. To create it, use the command `tabix -p vcf <VCF>.txt.gz`.

**`--vcf_include/--vcf_exclude`**: List of samples in a text file to include/exclude from the VCF. See this [example file](sample_data/vcf_include.txt).

**`--primer3`**: YAML configuration file for Primer3. All Primer3 arguments can be set here. See this [example file](sample_data/primer3_example.yaml.txt).

### Outputs
Two output files are produced:

**`<OUTPUT>.log`**: A log file which contain details on the number of valid primers found during each search and for each marker.

**`<OUTPUT>.bed`** Primer pairs reported in BED format for use with genome browsers. Names are composites of `<primer_id>-<goodness>` (see below).

**`<OUTPUT>.txt`** Primer pairs reported as a space-separated list with the following columns:

|Column|Description|
|---|---|
|`name`|Target name|
|`chr`|Chromosome|
|`start`|Primer start position in the genome|
|`end`|Primer end position in the genome|
|`direction`|Direction of the primer as F/R for forward/reverse, respectively|
|`assay_id`|ID of the primer pairs|
|`seq5_3`|Sequence of the primer in a 5'-3' direction|
|`seq_5_3_ambiguous`|Sequence of the primer in a 5'-3' direction with ambiguous nucleotides for mutations (no indels)|
|`primer_id`|Unique primer identification for each marker. Intended to ensure same primers are not purchased multiple time.|
|`goodness`|Heuristic goodness score based on multiple criteria. Maximum score is 10|
|`qcode`|Quality code containing warnings about the assay. Characters mean the following:<br>. =  No warnings <br>t = Bad TM<br>O = Off-targets<br>d = Heterodimerization<br>m/M = Mutations with allele frequencies >0/>0.1<br>i/I = Indels larger than 0/50 nucleotides|
|`length`|Primer length|
|`prod_size`|Expected PCR product size|
|`tm`|Predicted primer melting temperature (based on a NN thermodynamic model with SantaLucia et al, 1998 parameters)|
|`gc_content`|Percent GC in the primer sequence|
|`n_offtargets`|Number of possible genome-wide off-target PCR products. If larger than 5, then the number does not represent an exhaustive list.|
|`max_aaf`|Maximum alternative allele frequency across all mutations located in the primer|
|`indels`|Length of any indels located in the target PCR product|
|`offtargets`|Comma-separated list of expected off-target PCR products. If larger than 5, then this list is not exhaustive|
|`mutations`|Comma-separated list of mutations located in the primers and reported as [REF/ALT:AAF]|

### Optional output files
Optionally, a list of all subjects containing alternative alleles can be requested using the flag `--report_alts`, if a VCF file if provided. The list will be reported in `<OUTPUT>_altlist.txt`. Each mutation is listed in a pseudo-FASTA format as `>REFPOSALT`, similar to what is reported in the standard output but without allele frequencies. The first and second lines list all hets and homozygotes, respectively. This file can be used to investigate markers where no primers exempt of mutations exists.

### Example usage and tips
In the following example, primer pairs will be designed by considering both homologs and mutations within a selected subset population:

```
polyoligo-pcr sample_data/pcr_targets.txt out sample_data/blastdb --vcf sample_data/vcf.txt.gz --vcf_include sample_data/vcf_include.txt
```
