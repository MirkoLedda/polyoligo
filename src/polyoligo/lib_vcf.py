import vcf
import pandas as pd
import warnings
import numpy as np
import sys
from copy import deepcopy
import re


class VCF:
    def __init__(self, fp, fp_inc_samples="", fp_exc_samples=""):
        self.reader_fp = fp
        self.reader = vcf.Reader(filename=self.reader_fp)
        self.samples, self.samples_ix = self._parse_sample_files(fp_inc_samples, fp_exc_samples)

    def start_reader(self):
        self.reader = vcf.Reader(filename=self.reader_fp)

    def stop_reader(self):
        self.reader = None

    # noinspection PyProtectedMember
    def _parse_sample_files(self, fp_inc_samples, fp_exc_samples):

        if fp_inc_samples != "":
            inc_samples = read_list(fp_inc_samples)
        else:
            inc_samples = self.reader.samples

        if fp_exc_samples != "":
            exc_samples = read_list(fp_exc_samples)
        else:
            exc_samples = []

        # Exclude samples
        samples = []

        for s in inc_samples:
            if s not in exc_samples:
                samples.append(s)

        # Order the sample list and get sample indexes in the VCF file
        samples_ix = []

        for s in samples:
            if s in self.reader._sample_indexes:
                samples_ix.append(self.reader._sample_indexes[s])
            else:
                pass  # Skip sample if not in present in the VCF

        # Check we have at least 1 sample left
        if len(samples) < 1:
            raise NameError(
                "No samples selected in the VCF file. Please adjust --vcf_include and --vcf_exclude files.")

        six = np.argsort(samples_ix)
        samples = np.array(samples)[six]
        samples_ix = np.array(samples_ix)[six]

        return samples, samples_ix

    def fetch_kasp_markers(self, fp, chrom, **kwargs):

        start = None
        stop = None

        if "start" in kwargs.keys():
            start = kwargs["start"] - 1  # Position arguments is expected as 1-based

        if "stop" in kwargs.keys():
            stop = kwargs["stop"]

        if start < 0:
            start = 0

        try:
            mutations = self.reader.fetch(chrom, start=start, end=stop)
        except ValueError:
            sys.exit("ERROR: Region {}:{}-{} is not in the VCF file".format(chrom, start, stop))

        with open(fp, "a") as f:
            for mutation in mutations:
                if mutation.is_snp:
                    if len(mutation.ALT) == 1:
                        name = "{}@{}".format(mutation.CHROM, mutation.POS)
                        f.write("{} {} {} {} {}\n".format(mutation.CHROM,
                                                          mutation.POS,
                                                          name,
                                                          mutation.REF,
                                                          mutation.ALT[0],
                                                          ))

    # noinspection PyPep8Naming
    def fetch_genotypes(self, chrom, **kwargs):

        start = None
        stop = None

        if "start" in kwargs.keys():
            start = kwargs["start"] - 1  # Position arguments is expected as 1-based

        if "stop" in kwargs.keys():
            stop = kwargs["stop"]

        if start < 0:
            start = 0

        try:
            mutations = self.reader.fetch(chrom, start=start, end=stop)  # Position arguments is expected as 1-based
        except ValueError:
            sys.exit("ERROR: Region {}:{}-{} is not in the VCF file".format(chrom, start, stop))

        genotypes = {}  # Genotypes
        alleles = {}  # alleles

        for mutation in mutations:
            G = []
            alts = []  # List of possible alternative alleles

            sample_iter = iter(self.samples_ix)  # Iterator for the samples to include
            sid = next(sample_iter)

            for i, sample in enumerate(mutation.samples):

                if i == sid:
                    g = sample.gt_type

                    if g is None:
                        g = -1

                    if (g == 1) or (g == 2):
                        try:
                            alts.append(re.split(r'[/|]', sample.gt_bases)[1])
                        except IndexError:
                            sys.exit("Error: Unsupported VCF format. Please contact the developer.")

                    G.append(g)
                    try:
                        sid = next(sample_iter)
                    except StopIteration:
                        break

            if len(alts) > 0:  # Skip this mutation if none of the samples have the ALT allele
                genotypes[mutation.POS] = G
                genotypes[mutation.POS] = np.array(genotypes[mutation.POS], dtype=int)
                alleles[mutation.POS] = (mutation.REF, list(np.unique(alts)))

        genotypes = pd.DataFrame(genotypes, index=self.samples)

        # Compute allele frequencies
        aafs = compute_aaf(genotypes)

        # Create a list of mutations
        mutations = []

        i = 0
        for pos in genotypes:
            mutation = Mutation(
                chrom=chrom,
                pos=pos,
                ref=alleles[pos][0],
                alt=alleles[pos][1],
                aaf=aafs[i],
            )
            mutations.append(mutation)
            i += 1

        return genotypes, mutations

    def print_alternative_subjects(self, fp, chrom, start, stop):
        if fp != "":

            if start < 1:
                start = 1

            try:
                mutations = self.reader.fetch(chrom, start=start - 1, end=stop)
            except ValueError:
                sys.exit("ERROR: Region {}:{}-{} is not in the VCF file".format(chrom, start, stop))

            alleles = {}  # alleles
            alt_subjects = {}  # Subjects

            for mutation in mutations:
                alts = []  # List of possible alternative alleles
                alt_snames = {"het": [],
                              "hom": []}

                sample_iter = iter(self.samples_ix)  # Iterator for the samples to include
                sid = next(sample_iter)

                for i, sample in enumerate(mutation.samples):

                    if i == sid:
                        g = sample.gt_type

                        if g is None:
                            g = -1

                        if (g == 1) or (g == 2):
                            try:
                                alts.append(re.split(r'[/|]', sample.gt_bases)[1])

                                if g == 1:
                                    alt_snames["het"].append(sample.sample)
                                else:
                                    alt_snames["hom"].append(sample.sample)

                            except IndexError:
                                sys.exit("Error: Unsupported VCF format. Please contact the developer.")

                        try:
                            sid = next(sample_iter)
                        except StopIteration:
                            break

                if len(alts) > 0:  # Skip this mutation if none of the samples have the ALT allele
                    alleles[mutation.POS] = (mutation.REF, list(np.unique(alts)))
                    alt_subjects[mutation.POS] = deepcopy(alt_snames)

            # Print alternative subjects per variants
            with open(fp, "a") as f:
                for pos in alleles.keys():
                    mut_txt = "{}{:d}{}".format(alleles[pos][0], pos, "/".join(alleles[pos][1]))
                    f.write(">{}\n{}\n{}\n".format(
                        mut_txt,
                        "\t".join(alt_subjects[pos]["het"]),
                        "\t".join(alt_subjects[pos]["hom"])
                    ))


class Mutation:
    def __init__(self, chrom, pos, ref, alt, aaf=None, maf=None, indel_size=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = list(alt)

        # Get the max indel size for that SNP
        self.indel_size = indel_size
        if self.indel_size is None:
            self.indel_size = 0
            for a in self.alt:
                self.indel_size = np.max([self.indel_size, np.abs(len(ref) - len(a))])

        self.aaf = aaf
        self.maf = maf

        if aaf is not None:
            self.aaf = aaf
            self.maf = np.min([self.aaf, (1 - self.aaf)])


# noinspection PyPep8Naming
def compute_aaf(G):
    af = np.array([])
    if len(G) > 0:
        G = np.array(G, dtype=int)

        a1 = np.sum(G == 1, axis=0)
        a2 = np.sum(G == 2, axis=0)
        at = np.sum(G != -1, axis=0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            af = (a2 * 2 + a1) / (at * 2)
            af[np.isinf(af)] = 1  # In case no geno available for the SNP then we can consider it as not a variant

    return af


def read_list(fp):
    ls = []
    with open(fp, "r") as f:
        for line in f:
            ls.append(line.strip())

    return ls
