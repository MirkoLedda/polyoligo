from __future__ import print_function, division
from os.path import join
import numpy as np
import pandas as pd
import re
from copy import deepcopy
# noinspection PyPackageRequirements,PyPep8Naming
from Bio.SeqUtils import MeltingTemp as mt

from . import lib_utils, lib_blast


class ROI:
    def __init__(self, chrom=None, start=None, stop=None, blast_hook=None, malign_hook=None, vcf_hook=None, name=None,
                 marker=None, do_print_alt_subjects=False):
        self.name = None  # ROI name
        self.chrom = None  # Chromosome
        self.start = None  # Region start (1-based)
        self.stop = None  # Region end (1-based)
        self.blast_hook = None  # BLAST hook
        self.malign_hook = None  # Multiple alignment hook
        self.vcf_hook = None  # VCF hook
        self.seq = None  # Sequence of the region
        self.seq_alt = None  # Sequence of the region with the alternative allele if a marker is present
        self.seq_x = None  # Sequence with 'x' indicating the REF allele location
        self.n = None  # Length of the region
        self.G = None  # Possible genotypes within this region
        self.mutations = None  # List of vcf_lib.Mutation objects in the region
        self.has_marker = None  # Does the region contain a marker
        self.marker = None  # Target marker in the region
        self.p3_sequence_target = None  # For PRIMER3 design
        self.p3_sequence_included_maps = None  # For PRIMER3 design
        self.p3_sequence_excluded_regions = None  # For PRIMER3 design
        self.do_print_alt_subjects = do_print_alt_subjects  # Print alternative subjects?

        # Initialize attributes
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        self.blast_hook = blast_hook
        self.malign_hook = malign_hook
        self.vcf_hook = vcf_hook
        self.marker = marker

        if name is not None:
            self.name = name
        else:
            self.name = "{}:{}-{}".format(self.chrom, self.start, self.stop)

        if blast_hook is not None:
            query = [{
                "chr": self.chrom,
                "start": self.start,
                "stop": self.stop,
            }]

            seqs = self.blast_hook.fetch(query)
            self.seq = list(seqs.values())[0].upper()
            self.n = len(self.seq)

        if self.marker is not None:
            if (self.marker.pos1 >= self.start) and (self.marker.pos1 <= self.stop):
                self.seq_alt = list(self.seq)
                astart = self.marker.pos1 - self.start
                astop = self.marker.pos1 - self.start + self.marker.n
                lseq = self.seq_alt[:astart]
                rseq = self.seq_alt[astop:]
                self.seq_alt = lseq + list(self.marker.alt) + rseq
                self.seq_alt = "".join(self.seq_alt)
                self.seq_x = lseq + ["x"] * self.marker.n + rseq
                self.seq_x = "".join(self.seq_x)
            else:
                self.seq_alt = self.seq
                self.seq_x = self.seq

    def map_homologs(self, inner_len, outer_len, min_align_id, min_align_len):

        lp, rp = lib_utils.get_n_padding(x=self.n, n=inner_len)

        query = [{
            "chr": self.chrom,
            "start": int(self.start - lp),
            "stop": int(self.stop + rp),
        }]

        seqs = self.blast_hook.fetch(query)
        for k, seq in seqs.items():
            k = k
            nseq = {k: seq.upper()}
            break

        # Pad sequences with N's if needed
        exp_seq_len = query[0]["stop"] - query[0]["start"]
        # noinspection PyUnboundLocalVariable
        if len(nseq[k]) < exp_seq_len:
            if query[0]["start"] == 1:
                nseq[k] = lib_utils.padding_left(x=nseq[k],
                                                 n=exp_seq_len)  # left padding
            else:
                nseq[k] = lib_utils.padding_right(x=nseq[k],
                                                  n=exp_seq_len)  # right padding

        fp_query = join(self.blast_hook.temporary, "{}_homologs_blast_in.txt".format(self.name))
        lib_blast.write_fasta(nseq, fp_out=fp_query)

        # Run BLAST search
        fp_blast_out = join(self.blast_hook.temporary, "{}_homologs_blast_out.txt".format(self.name))
        self.blast_hook.blastn(
            fp_query=fp_query,
            fp_out=fp_blast_out,
            perc_identity=min_align_id,
            max_hsps=100,  # Max 100 HSPs returned as above 50 valid sequences we won't do multialignment
            dust="no",
            soft_masking="false",
            outfmt='"6 std qseq sseq slen"',
        )
        fp_blast_out2 = join(self.blast_hook.temporary, "{}_homologs_blast_out_pruned.txt".format(self.name))

        # Prune to keep only valid hit lengths
        self.blast_hook.prune_blastn_output(fp_blast_out,
                                            fp_blast_out2,
                                            min_align_len=min_align_len)

        # Retrieve the sequences of the homologs
        homologs_exceeded = False
        df = self.blast_hook.parse_blastn_output(fp_blast_out2)

        hseqs = []
        if df is not None:  # Can happen if the query is all N's
            for _, row in df.iterrows():
                is_target = False
                if row["schr"] == self.chrom:
                    if (self.start >= row["sstart"]) and (self.stop <= row["sstop"]):
                        is_target = True

                if not is_target:
                    if row["sstart"] > row["sstop"]:
                        strand = "minus"
                    else:
                        strand = "plus"

                    hseq = {
                        "chrom": row["schr"],
                        "start": np.min([row["sstart"], row["sstop"]]),
                        "stop": np.max([row["sstart"], row["sstop"]]),
                        "strand": strand
                    }
                    hseqs.append(hseq)

        if len(hseqs) > 50:
            homologs_exceeded = True

        # Retrieve outer sequences of the target and the homologs
        tname = None
        seqs = {}

        lp, rp = lib_utils.get_n_padding(x=self.n, n=outer_len)
        query = [{
            "chr": self.chrom,
            "start": int(self.start - lp),
            "stop": int(self.stop + rp),
        }]
        for k, seq in self.blast_hook.fetch(query).items():
            tname = k
            seqs[k] = seq.upper()
            break

        # For runtime optimization - If more than 50 homoloous sequences were found, then do not run fetching and
        # multialignment and instead make all nucleotides not specific (using twice the same target sequence)
        if not homologs_exceeded:
            for hseq in hseqs:
                n = np.abs(hseq["stop"] - hseq["start"])
                lp, rp = lib_utils.get_n_padding(x=n, n=outer_len)

                if hseq["strand"] == "plus":
                    query = [{
                        "chr": hseq["chrom"],
                        "start": int(hseq["start"] - lp),
                        "stop": int(hseq["stop"] + rp),
                        "strand": hseq["strand"],
                    }]
                else:
                    query = [{
                        "chr": hseq["chrom"],
                        "start": int(hseq["stop"] - rp),
                        "stop": int(hseq["start"] + lp),
                        "strand": hseq["strand"],
                    }]
                fetched_seqs = self.blast_hook.fetch(query)
                for k, seq in fetched_seqs.items():
                    seqs[k] = seq
                    break

        # Multiple alignment
        fp_malign_in = join(self.blast_hook.temporary, "{}.fa".format(self.name))
        fp_malign_out = join(self.blast_hook.temporary, "{}.afa".format(self.name))
        if homologs_exceeded:
            mock_seqs = {
                tname: seqs[tname],
                "mock": seqs[tname],
            }
            lib_blast.write_fasta(mock_seqs, fp_out=fp_malign_in)
        else:
            lib_blast.write_fasta(seqs, fp_out=fp_malign_in)

        self.malign_hook.align_fasta(fp=fp_malign_in, fp_out=fp_malign_out)

        # Map homologs
        maps = dict()

        seqs = lib_blast.read_fasta(fp_malign_out)
        tseq = seqs[tname]
        tlen = len(seqs[tname].replace("-", ""))
        del seqs[tname]

        n_homeo = len(seqs)
        match_arr = np.tile(0, (tlen, n_homeo))

        # Count mismatches
        homeo_id = 0
        for _, seq in seqs.items():
            j = 0  # Index in the original sequence without gaps
            for i, nuc in enumerate(tseq):
                if nuc == "-":
                    continue
                elif nuc == "N":  # Ns are considered matches
                    match_arr[j, homeo_id] += 1
                    j += 1
                elif nuc == seq[i]:
                    match_arr[j, homeo_id] += 1
                    j += 1
                else:
                    j += 1
            homeo_id += 1
            break

        # Compile the mismatch counts into maps
        if n_homeo == 0:
            maps["all"] = np.repeat(True, tlen)
            maps["partial_mism"] = np.repeat(True, tlen)
            maps["mism"] = np.repeat(True, tlen)
        else:
            match_cnts = np.sum(match_arr, axis=1)
            maps["all"] = np.repeat(True, tlen)
            maps["partial_mism"] = match_cnts != n_homeo
            maps["mism"] = match_cnts == 0

        chrom = tname.strip().split(":")[0]
        start = int(tname.strip().split(":")[1].split("-")[0])
        stop = int(tname.strip().split(":")[1].split("-")[1])

        hroi = ROI(
            chrom=chrom,
            start=start,
            stop=stop,
            blast_hook=self.blast_hook,
            malign_hook=self.malign_hook,
            vcf_hook=self.vcf_hook,
            name=self.name,
            marker=self.marker,
            do_print_alt_subjects=self.do_print_alt_subjects,
        )

        hroi.p3_sequence_included_maps = maps

        return hroi

    def print_alt_subjects(self, fp):
        if self.do_print_alt_subjects and (self.vcf_hook is not None):
            self.vcf_hook.print_alternative_subjects(
                fp=fp,
                chrom=self.chrom,
                start=self.start,
                stop=self.stop,
            )

    def upload_mutations(self):
        self.mutations = []
        if self.vcf_hook:
            self.G, self.mutations = self.vcf_hook.fetch_genotypes(
                chrom=self.chrom,
                start=self.start,
                stop=self.stop,
            )

            # If a marker is present, pop it from the list of mutations
            if self.marker is not None:
                i = np.where(np.array(self.G.columns) == self.marker.pos1)[0]
                if len(i) == 1:
                    _ = self.G.pop(self.marker.pos1)
                    _ = self.mutations.pop(i[0])
                else:
                    pass


class Marker:
    """Hold a single marker information."""

    def __init__(self,  name=None, chrom=None, pos=None, ref_allele=None, alt_allele=None):
        self.name = None  # Marker name
        self.chrom = None  # Chromosome
        self.pos = None  # Position 0-based
        self.pos1 = None  # Position 1-based
        self.ref = None  # Reference allele
        self.alt = None  # Alternative allele
        self.variant = None  # Variant notation
        self.n = None  # Maximum length of the mutation

        # Initialize attributes
        self.name = name
        self.chrom = chrom
        self.pos = int(pos)
        self.pos1 = int(pos) + 1
        self.ref = ref_allele
        self.alt = alt_allele

        # Handle indel '*' notation
        if self.ref == "*":
            self.ref = ""

        if self.alt == "*":
            self.alt = ""

        self.variant = "[{}/{}]".format(self.ref, self.alt)
        self.n = len(self.ref)

        if (name is None) or (name == "."):
            self.name = "{}:{}-{}".format(self.chrom, self.pos1, self.pos1)
        self.name = self.name.replace("_", "-")  # Underscores in the marker name will mess with the algorithm

    def assert_marker(self, blast_hook):
        if self.ref != "":
            query = [{
                "chr": self.chrom,
                "start": int(self.pos1),
                "stop": int(self.pos1 + self.n - 1),
            }]

            seq = "N"
            for _, seq in blast_hook.fetch(query).items():
                seq = seq.upper()
                break

            err_msg = None
            if self.ref != seq:
                err_msg = "REF allele in the marker file does not match the genomic REF allele.\n" \
                          "        SNP ID: {} | Marker {} vs Reference {}\n" \
                          "        Please double check your marker file.".format(self.name, self.ref, seq)

            return err_msg


class Region(ROI):
    def __init__(self, left_roi=None, right_roi=None, **kwargs):
        super().__init__(**kwargs)
        self.left_roi = None
        self.right_roi = None

        # Initialize attributes
        self.left_roi = left_roi
        self.right_roi = right_roi

    def merge_with_primers(self):

        # Spawn a new combined region
        region = Region(
            chrom=self.chrom,
            start=self.left_roi.start,
            stop=self.right_roi.stop,
            blast_hook=self.blast_hook,
            malign_hook=self.malign_hook,
            vcf_hook=self.vcf_hook,
            name=self.name,
            marker=self.marker,
            do_print_alt_subjects=self.do_print_alt_subjects,
            left_roi=self.left_roi,
            right_roi=self.right_roi,
        )

        region.p3_sequence_target = [[self.start - region.left_roi.start, self.n]]

        # Combine the inclusion maps
        ixs = np.arange(region.start, region.stop+1)
        ix_left = np.arange(region.left_roi.start, region.left_roi.stop+1)
        ix_center = np.arange(self.start, self.stop+1)
        ix_right = np.arange(region.right_roi.start, region.right_roi.stop+1)
        ix_left = np.intersect1d(ixs, ix_left, return_indices=True)[1]
        ix_center = np.intersect1d(ixs, ix_center, return_indices=True)[1]
        ix_right = np.intersect1d(ixs, ix_right, return_indices=True)[1]
        base_map = np.repeat(0, region.n)

        region.p3_sequence_included_maps = {}
        for k in region.left_roi.p3_sequence_included_maps.keys():
            mmap = deepcopy(base_map)
            imap = deepcopy(base_map)
            mmap[ix_left] += ~region.left_roi.p3_sequence_included_maps[k]
            imap[ix_left] += 1
            mmap[ix_center] += 1
            imap[ix_center] += 1
            mmap[ix_right] += ~region.right_roi.p3_sequence_included_maps[k]
            imap[ix_right] += 1
            mmap[imap == 0] += 1  # Exclude unmapped areas

            region.p3_sequence_included_maps[k] = mmap == 0

        if self.left_roi.mutations is not None:
            region.mutations = region.left_roi.mutations + region.right_roi.mutations

        return region


class PCRproduct:
    def __init__(self, roi, pp):
        self.roi = Region(
            chrom=roi.chrom,
            start=pp.primers["F"].start,
            stop=pp.primers["R"].stop,
            blast_hook=roi.blast_hook,
            malign_hook=roi.malign_hook,
            vcf_hook=roi.vcf_hook,
            name=roi.marker.name,
            marker=roi.marker,
            do_print_alt_subjects=roi.do_print_alt_subjects,
        )
        self.pp = pp

        self.enzymes = None  # For CAPS
        self.hrm = None  # For HRM

    def will_it_cut(self, enzymes, min_fragment_size):
        self.enzymes = []

        for enzyme in enzymes:
            m = re.findall(re.compile(enzyme["regex"]["F"]), self.roi.seq)
            mi = re.findall(re.compile(enzyme["regex"]["R"]), self.roi.seq)
            n_cut = np.max([len(m), len(mi)])

            ma = re.findall(re.compile(enzyme["regex"]["F"]), self.roi.seq_alt)
            mai = re.findall(re.compile(enzyme["regex"]["R"]), self.roi.seq_alt)
            n_cut_alt = np.max([len(ma), len(mai)])

            fragl = 0
            fragr = 0
            all_cut = ""

            if (n_cut == 1) and (n_cut_alt == 0):
                for m in re.finditer(re.compile(enzyme["regex"]["F"]), self.roi.seq):
                    cut_pos = m.span()[0] + int(np.ceil((m.span()[1]-m.span()[0]))/2)
                    fragl = cut_pos
                    fragr = self.roi.n - cut_pos
                    all_cut = "REF"

                for m in re.finditer(re.compile(enzyme["regex"]["R"]), self.roi.seq):
                    cut_pos = m.span()[0] + int(np.ceil((m.span()[1]-m.span()[0]))/2)
                    fragl = cut_pos
                    fragr = self.roi.n - cut_pos
                    all_cut = "REF"

            if (n_cut == 0) and (n_cut_alt == 1):
                for m in re.finditer(re.compile(enzyme["regex"]["F"]), self.roi.seq_alt):
                    cut_pos = m.span()[0] + int(np.ceil((m.span()[1] - m.span()[0])/2))
                    fragl = cut_pos
                    fragr = self.roi.n - cut_pos
                    all_cut = "ALT"

                for m in re.finditer(re.compile(enzyme["regex"]["R"]), self.roi.seq_alt):
                    cut_pos = m.span()[0] + int(np.ceil((m.span()[1] - m.span()[0])) / 2)
                    fragl = cut_pos
                    fragr = self.roi.n - cut_pos
                    all_cut = "ALT"

            if (fragl >= min_fragment_size) and (fragr >= min_fragment_size):
                enzyme["fragl"] = fragl
                enzyme["fragr"] = fragr
                enzyme["allele"] = all_cut
                self.enzymes.append(enzyme)

    def get_hrm_temps(self):
        # noinspection PyDictCreation
        self.hrm = {}
        self.hrm["ref"] = get_tm(self.roi.seq)
        self.hrm["alt"] = get_tm(self.roi.seq_alt)
        self.hrm["delta"] = np.abs(self.hrm["ref"] - self.hrm["alt"])

        for k, v in self.hrm.items():
            if k == "delta":
                self.hrm[k] = lib_utils.round_tidy(v, 2)
            else:
                self.hrm[k] = lib_utils.round_tidy(v, 4)


def read_markers(fp):
    """Reads a text file of input markers and returns a list of Marker objects.

    File should be tab delimited and no headers and contain the information columns:
    CHR POS NAME REF ALT

    """

    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    df.columns = ["chr", "pos", "name", "ref", "alt"]

    markers = []
    for _, row in df.iterrows():
        markers.append(Marker(
            chrom=row["chr"],
            pos=row["pos"] - 1,  # to 0-based indexing
            ref_allele=row["ref"],
            alt_allele=row["alt"],
            name=row["name"],
        ))

    return markers


def read_regions(fp):
    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    df.columns = ["region", "name"]

    regions = []
    for _, row in df.iterrows():
        chrom = row["region"].strip().split(":")[0]
        start = row["region"].strip().split(":")[1].split("-")[0]
        stop = row["region"].strip().split(":")[1].split("-")[1]

        regions.append(Region(
            name=row["name"],
            chrom=chrom,
            start=start,
            stop=stop,
        ))

    return regions


def get_tm(seq, **kwargs):
    return mt.Tm_NN(seq, strict=False, **kwargs)
