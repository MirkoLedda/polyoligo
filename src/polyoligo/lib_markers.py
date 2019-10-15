from __future__ import print_function, division
from os.path import join
import numpy as np
import pandas as pd

from . import lib_utils, lib_blast


class ROI:
    def __init__(self, chrom=None, start=None, stop=None, blast_hook=None, malign_hook=None, vcf_hook=None, name=None,
                 marker=None):
        self.name = None  # ROI name
        self.chrom = None  # Chromosome
        self.start = None  # Region start (1-based)
        self.stop = None  # Region end (1-based)
        self.blast_hook = None  # BLAST hook
        self.malign_hook = None  # Multiple alignment hook
        self.vcf_hook = None  # VCF hook
        self.seq = None  # Sequence of the region
        self.seq_alt = None  # Sequence of the region with the alternative allele if a marker is present
        self.n = None  # Length of the region
        self.G = None  # Possible genotypes within this region
        self.mutations = None  # List of vcf_lib.Mutation objects in the region
        self.has_marker = None  # Does the region contain a marker
        self.marker = None  # Target marker in the region

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

        query = [{
            "chr": self.chrom,
            "start": self.start,
            "stop": self.stop,
        }]

        seqs = self.blast_hook.fetch(query)
        self.seq = list(seqs.values())[0].upper()
        self.n = len(self.seq)

        if self.marker is not None:
            self.seq_alt = np.array(list(self.seq))
            self.seq_alt[
            (self.marker.pos1 - self.start):(self.marker.pos1 - self.start + len(self.marker.ref))
            ] = self.marker.alt
            self.seq_alt = "".join(self.seq_alt)

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

        # Retrieve the center of the sequence of the homologs
        df = self.blast_hook.parse_blastn_output(fp_blast_out2)

        hseqs = []
        for _, row in df.iterrows():
            is_target = False
            if row["schr"] == self.chrom:
                if (row["sstart"] >= self.start) and (row["sstop"] <= self.stop):
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

        # Multiple alignment
        fp_malign_in = join(self.blast_hook.temporary, "{}.fa".format(self.name))
        fp_malign_out = join(self.blast_hook.temporary, "{}.afa".format(self.name))
        if len(seqs) >= 50:
            # For runtime optimization - If more than 50 homoloous sequences were found, then do not run
            # multialignment and instead make all nucleotides not specific (using twice the same target sequence)
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

        mseqs = lib_blast.read_fasta(fp_malign_out)
        tseq = seqs[tname]
        tlen = len(seqs[tname])
        del seqs[tname]

        n_homeo = len(mseqs)
        match_arr = np.tile(np.nan, (tlen, n_homeo))

        # Count mismatches
        homeo_id = 0
        for chrom, seq in mseqs.items():
            j = 0  # Index in the original sequence without gaps
            for i, nuc in enumerate(tseq):
                if nuc == "-":
                    continue
                elif nuc == "N":  # Ns are considered matches
                    match_arr[j, homeo_id] = 1
                    j += 1
                elif nuc == seq[i]:
                    match_arr[j, homeo_id] = 1
                    j += 1
                else:
                    match_arr[j, homeo_id] = 0
                    j += 1
            homeo_id += 1

        # Compile the mismatch counts into maps
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
        )

        return hroi, maps

    def get_alt_seq(self):
        pass  # TODO

    def upload_mutations(self):
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
                    # logger.warning("Marker {} ({}:{:d}) is not present in the VCF file.".format(self.name,
                    #                                                                             self.chrom,
                    #                                                                             self.pos1))
        else:
            self.mutations = []


class Marker:
    """Hold a single marker information."""

    def __init__(self, chrom=None, pos=None, ref_allele=None, alt_allele=None, name=None):
        self.name = None  # Marker name
        self.chrom = None  # Chromosome
        self.pos = None  # Position 0-based
        self.pos1 = None  # Position 1-based
        self.ref = None  # Reference allele
        self.alt = None  # Alternative allele
        self.variant = None  # Variant notation

        # Initialize attributes
        self.chrom = chrom
        self.pos = int(pos)
        self.pos1 = int(pos) + 1
        self.ref = ref_allele
        self.alt = alt_allele
        self.variant = "[{}/{}]".format(self.ref, self.alt)

        if (name is None) or (name == "."):
            self.name = "{}:{}-{}".format(self.chrom, self.pos1, self.pos1)
        self.name = self.name.replace("_", "-")  # Underscores in the marker name will mess with the algorithm


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
