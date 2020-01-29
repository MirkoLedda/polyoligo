# `polyoligo-caps`

## Usage

### General usage and available options

```
polyoligo-caps <INPUT> <OUTPUT> <BLASTDB> <OPTIONS>
```

For a list of all available options and their descriptions, type:

```
polyoligo-caps -h
```

> Recommendations (when applicable) are given in the option caption. Note that switches, i.e. boolean options that do not need arguments, have defaults set to `False`.

### Example usage
In the following example, CAPS primers will be designed by considering both homologs and mutations within a selected subset population:

```
polyoligo-caps sample_data/markers.txt out sample_data/blastdb --vcf sample_data/vcf.txt.gz --vcf_include sample_data/vcf_include.txt
```

### Inputs
The software requires three mandatory inputs:

**`<INPUT>`**: A text file containing target markers as a list of [CHR POS NAME REF ALT]. See this [example file](sample_data/markers.txt). Deletions are denoted by `*` and `.` indicates no name . For example `Fvb1-1 10000 . G *` is a valid entry.

**`<OUTPUT>`**: The base name of the output files.

**`<FASTA/BLASTDB>`**: A FASTA file and/or a BLAST database to use as the reference genome. Both file types can be provided by using the same basename. If either is provided, then a conversion will automatically be made to obtain both file types.

Optional files include:

**`--vcf`**: A VCF file to design primers considering mutations (both SNPs and indels).

> Note that a tabix index file, created using the Samtools' [`tabix`](http://www.htslib.org/doc/tabix.html) binary is required. To create it, use the command `tabix -p vcf <VCF>.txt.gz`.

**`--vcf_include/--vcf_exclude`**: List of samples in a text file to include/exclude from the VCF. See this [example file](sample_data/vcf_include.txt).

**`--primer3`**: YAML configuration file for Primer3. All Primer3 arguments can be set here. See this [example file](sample_data/primer3_example.yaml).

**`--enzymes`**: List of restriction enzymes to consider. Can be useful to restrict the search to in stock enzymes only. See this [example file](sample_data/enzymes.txt).

### Outputs
Three output files are produced:

**`<OUTPUT>.log`**: A log file which contain details on the number of valid primers found during each search and for each marker.

**`<OUTPUT>.bed`** CAPS primers reported in BED format for use with genome browsers. Names are composites of `<primer_id>-<goodness>` (see below).

**`<OUTPUT>.txt`** CAPS primers reported as a space-separated list with the following columns:

|Column|Description|
|---|---|
|`marker`|SNP ID/label|
|`chr`|Chromosome|
|`pos`|Position|
|`ref`|Reference allele|
|`alt`|Alternative allele|
|`enzymes`|List of possible restriction enzymes|
|`enzyme_suppliers`|Restriction enzyme suppliers using the REBASE encoding. See this [file](src/polyoligo/data/REBASE_suppliers.txt) for the letter code|
|`start`|Primer start position in the genome|
|`end`|Primer end position in the genome|
|`direction`|Direction of the primer as F/R for forward/reverse, respectively|
|`assay_id`|ID of the assay|
|`seq5_3`|Sequence of the primer in a 5'-3' direction|
|`seq_5_3_ambiguous`|Sequence of the primer in a 5'-3' direction with ambiguous nucleotides for mutations (no indels)|
|`primer_id`|Unique primer identification for each marker. Intended to ensure same primers are not purchased multiple time.|
|`goodness`|Heuristic goodness score based on multiple criteria. Maximum score is 10|
|`qcode`|Quality code containing warnings about the assay. Characters mean the following:<br>. =  No warnings <br>t = Bad TM<br>O = Off-targets<br>d = Heterodimerization<br>m/M = Mutations with allele frequencies >0/>0.1<br>i/I = Indels larger than 0/50 nucleotides|
|`length`|Primer length|
|`prod_size`|Expected PCR product size|
|`fragment_left`|5'-end approximate fragment length (assumes the enzyme cuts at the marker site which isn't always true)|
|`fragment_right`|3'-end approximate fragment length (assumes the enzyme cuts at the marker site which isn't always true)|
|`tm`|Predicted primer melting temperature (based on a NN thermodynamic model with SantaLucia et al, 1998 parameters)|
|`gc_content`|Percent GC in the primer sequence|
|`n_offtargets`|Number of possible genome-wide off-target PCR products. If larger than 5, then the number does not represent an exhaustive list.|
|`max_aaf`|Maximum alternative allele frequency across all mutations located in the primer|
|`indels`|Length of any indels located in the target PCR product|
|`offtargets`|Comma-separated list of expected off-target PCR products. If larger than 5, then this list is not exhaustive|
|`mutations`|Comma-separated list of mutations located in the primers and reported as [REF/ALT:AAF]|
|`PCR_product`|Sequence of the entire PCR product with `x` masking the location of the REF allele|
