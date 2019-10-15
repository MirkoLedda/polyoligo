from __future__ import print_function, division
import pandas as pd
import logging
from os.path import join, exists
import sys
import numpy as np
import os

from . import lib_blast, lib_utils, _lib_kasp

logger = logging.getLogger(__name__)  # Initialize the logger


class Marker:
    """Hold informations for a Marker."""


    #
    # def upload_mutations(self, vcf_obj):
    #     if vcf_obj:
    #         self.G, self.mutations = vcf_obj.fetch_genotypes(
    #             chrom=self.chrom,
    #             start=self.start,
    #             stop=self.stop,
    #         )
    #         # Pop the mutation corresponding to the target
    #         i = np.where(np.array(self.G.columns) == self.pos1)[0]
    #         if len(i) == 1:
    #             _ = self.G.pop(self.pos1)
    #             _ = self.mutations.pop(i[0])
    #         else:
    #             pass
    #             # logger.warning("Marker {} ({}:{:d}) is not present in the VCF file.".format(self.name,
    #             #                                                                             self.chrom,
    #             #                                                                             self.pos1))
    #     else:
    #         self.mutations = []

    # def print_alt_subjects(self, vcf_obj, fp):
    #     if vcf_obj:
    #         vcf_obj.print_alternative_subjects(
    #             fp=fp,
    #             chrom=self.chrom,
    #             start=self.start,
    #             stop=self.stop,
    #         )

    # def upload_sequence(self, seq):
    #     self.seq = seq  # Sequence of the region
    #     self.n = len(self.seq)


class Markers:
    """Holds a list of Marker objects."""

    def __init__(self, blast_db=None, **kwargs):
        self.markers = None
        self.blast_db = blast_db
        self.MARKER_FLANKING_N = None
        self.MIN_ALIGN_LEN = None
        self.MIN_ALIGN_ID = None
        self.HOMOLOG_FLANKING_N = None

        # Set attributes from kwargs
        lib_utils.kwargs2attr_deep(self, kwargs)

    def __len__(self):
        return len(self.markers)

    def read_markers(self, fp):
        """Reads a text file of input markers and returns a list of Marker objects.

        File should be tab delimited and no headers and contain the information columns:
        CHR POS NAME REF ALT

        """

        df = pd.read_csv(fp, delim_whitespace=True, header=None)
        df.columns = ["chr", "pos", "name", "ref", "alt"]

        self.markers = []
        for _, row in df.iterrows():
            self.markers.append(Marker(
                chrom=row["chr"],
                pos=row["pos"] - 1,  # to 0-based indexing
                ref_allele=row["ref"],
                alt_allele=row["alt"],
                name=row["name"],
                HOMOLOG_FLANKING_N=self.HOMOLOG_FLANKING_N,
            ))

    def upload_fasta_ids(self, fasta_ids):
        for marker in self.markers:
            fasta_name = fasta_ids[marker.name]
            marker.upload_fasta_name(fasta_name)

    def upload_mutations(self, *args, **kwargs):
        for marker in self.markers:
            marker.upload_mutations(*args, **kwargs)

    def print_alt_subjects(self, *args, **kwargs):
        for marker in self.markers:
            marker.print_alt_subjects(*args, **kwargs)

    def get_queries_for_homolog_flanking(self, fp, snp_pos, padding):
        queries_dict = {}
        targets = {}

        for marker in self.markers:
            fpm = join(self.blast_db.temporary, marker.name) + ".hfa"
            if exists(fpm):
                os.remove(fpm)

        # Subset the input BLAST file into marker files
        with open(fp, "r") as f:
            for line in f:
                marker_name = line.split("_")[0]
                fpm = join(self.blast_db.temporary, marker_name) + ".hfa"

                with open(fpm, "a") as f_out:
                    f_out.write(line)

        for marker in self.markers:
            marker_found = False
            try:
                # Parse the BLAST output file and return a pandas.DataFrame
                fpm = join(self.blast_db.temporary, marker.name) + ".hfa"
                df = self.blast_db.parse_blastn_output(fpm)
                os.remove(fpm)

                # Add snp names and chromosomes to df
                df["snp"] = df.qname.str.split("_").str[0].tolist()
                df["qchr"] = df.qname.str.split("_").str[1].tolist()
                df["qpos"] = pd.to_numeric(df.qname.str.split("_").str[2].tolist())

                queries = []
                target = None

                df_sub = df[(df["qstart"] <= snp_pos) &
                            (df["qstop"] >= snp_pos)]

                for _, row in df_sub.iterrows():

                    # Real position of the SNP in the query accounting for gaps
                    qpos = snp_pos - row["qstart"]
                    for j, nuc in enumerate(row["qseq"]):
                        if nuc == "-":
                            qpos += 1
                        if qpos == j:
                            break

                    # Number of gaps until the SNP in the subject
                    sgap = row["sseq"][0:qpos].count("-")

                    if row["sstart"] > row["sstop"]:
                        strand = "minus"
                        spos = row["sstart"] - (qpos - sgap)
                    else:
                        strand = "plus"
                        spos = row["sstart"] + (qpos - sgap)

                    q = {
                        "name": marker.name,
                        "chr": row["schr"],
                        "start": spos - padding,
                        "stop": spos + padding,
                        "strand": strand,
                    }

                    queries.append(q)

                    # Determine which sequence correspond to the marker
                    if (spos == row["qpos"]) and (str(row["schr"]) == str(marker.chrom)):
                        target = "{}:{:d}-{:d}".format(q["chr"], q["start"], q["stop"])
                        marker_found = True

                    # Limit to 50 homologous sequences as we won't do multialignment for more sequences
                    if marker_found and (len(queries) >= 51):
                        break

                queries_dict[marker.name] = queries
                targets[marker.name] = target
            except FileNotFoundError:
                queries_dict[marker.name] = []
                targets[marker.name] = None

            if not marker_found:  # In rare instances BLAST did not return the target sequence
                q = {
                    "name": marker.name,
                    "chr": marker.chrom,
                    "start": marker.pos1 - padding,
                    "stop": marker.pos1 + padding,
                    "strand": "plus",
                }
                target = "{}:{:d}-{:d}".format(q["chr"], q["start"], q["stop"])
                targets[marker.name] = target
                queries_dict[marker.name].append(q)

        return queries_dict, targets

    def get_marker_flanks(self):

        # Make the queries for BLASTDBCMD
        queries = {}

        for marker in self.markers:
            queries[marker.name] = {
                "chr": marker.chrom,
                "start": marker.pos1 - self.MARKER_FLANKING_N,
                "stop": marker.pos1 + self.MARKER_FLANKING_N,
            }

        # Fetch sequences flanking the markers
        exp_seq_len = self.MARKER_FLANKING_N * 2 + 1
        seqs = self.blast_db.fetch(list(queries.values()))
        seqsr = {}
        for i, marker in enumerate(self.markers):
            q = queries[marker.name]
            qname = "{}:{}-{}".format(q["chr"], q["start"], q["stop"])
            seqsr[marker.blast_name] = seqs[qname]
            # Pad sequences with N's if needed
            if len(seqsr[marker.blast_name]) < exp_seq_len:
                if q["start"] == 1:
                    seqsr[marker.blast_name] = lib_utils.padding_left(x=seqsr[marker.blast_name],
                                                                      n=exp_seq_len)  # left padding
                else:
                    seqsr[marker.blast_name] = lib_utils.padding_right(x=seqsr[marker.blast_name],
                                                                       n=exp_seq_len)  # right padding

            # # Assert that the REF alleles in the genomic reference matches the ones provided as input
            # if (seqsr[marker.blast_name][self.MARKER_FLANKING_N] != marker.ref) and (marker.ref != "X"):
            #     logger.error("REF allele in the marker file does not match the genomic REF allele.\n"
            #                  "SNP ID: {} | Marker {} vs Reference {}. "
            #                  "Please double check your marker file.".format(marker.name, marker.ref,
            #                                                                 seqs[qname][self.MARKER_FLANKING_N]))
            #     sys.exit()

        return seqsr

    def find_homologs(self, seqs):

        # Write a file with queries for BLAST
        fp_query = join(self.blast_db.temporary, "homologs_blast_in.txt")
        lib_blast.write_fasta(seqs=seqs, fp_out=fp_query)

        # Run BLAST
        fp_blast_out = join(self.blast_db.temporary, "homologs_blast_out.txt")
        self.blast_db.blastn(
            fp_query=fp_query,
            fp_out=fp_blast_out,
            perc_identity=self.MIN_ALIGN_ID,
            max_hsps=100,  # Max 100 HSPs returned as above 50 valid sequences we won't do multialignment
            dust="no",
            soft_masking="false",
            outfmt='"6 std qseq sseq slen"',
        )
        fp_blast_out2 = join(self.blast_db.temporary, "homologs_blast_out_pruned.txt")

        # Prune to keep only valid hit lengths
        self.blast_db.prune_blastn_output(fp_blast_out,
                                          fp_blast_out2,
                                          min_align_len=self.MIN_ALIGN_LEN)

        # Retrieve homolog sequences and write them to a file
        queries_dict, fasta_ids = self.get_queries_for_homolog_flanking(
            fp_blast_out2,
            snp_pos=self.MARKER_FLANKING_N + 1,  # 1-based
            padding=self.HOMOLOG_FLANKING_N,
        )

        for snp_id, queries in queries_dict.items():
            _ = self.blast_db.fetch(queries, fp_out=join(self.blast_db.temporary, snp_id + ".fa"))

        self.upload_fasta_ids(fasta_ids)  # Store FASTA ID of the markers in each FASTA file

    def write_kasp_reports(self, fp_base_out):
        _lib_kasp.print_report_header(fp_base_out + ".txt")
        with open(fp_base_out + ".txt", "a") as f:
            for marker in self.markers:
                with open(join(self.blast_db.temporary, marker.name + ".txt"), "r") as f_marker:
                    for line in f_marker:
                        f.write(line)

        # BED file
        line_memory = []
        with open(fp_base_out + ".bed", "w") as f:
            for marker in self.markers:
                with open(join(self.blast_db.temporary, marker.name + ".bed"), "r") as f_marker:
                    for line in f_marker:
                        if line not in line_memory:
                            line_memory.append(line)
                            f.write(line)


if __name__ == "__main__":
    pass
