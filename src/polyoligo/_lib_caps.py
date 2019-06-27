from __future__ import print_function, division
# noinspection PyPackageRequirements
from Bio.Seq import Seq
import numpy as np
# noinspection PyUnresolvedReferences
from primer3.thermoanalysis import ThermoAnalysis
from copy import deepcopy
import logging
from os.path import join
from os import path
import os
import sys
import json
import re

from . import lib_blast, lib_utils, lib_primer3, _lib_pcr, _lib_kasp

ENZYME_FILENAME = join(os.path.dirname(__file__), "data/type2.json")

# Initialize the logger
logger = logging.getLogger(__name__)


class CAPS:
    def __init__(self, marker, included_enzymes=None):

        self.marker = marker
        self.marker_pos = None
        self.seqs = None
        self.substrings = None
        self.enzymes = []
        self.valid_enzymes = None
        self.valid_enzymes_list = None

        enzymes = self.load_enzymes()

        if included_enzymes is not None:
            for enzyme in enzymes:
                if enzyme["name"] in included_enzymes:
                    self.enzymes.append(enzyme)
        else:
            self.enzymes = enzymes

        self._get_substrings()

    @staticmethod
    def load_enzymes():
        # Load restriction enzymes
        with open(ENZYME_FILENAME, "r") as f:
            enzymes = json.load(f)

        return enzymes

    def _get_recognition_site_max_length(self):
        max_n = 0
        for enzyme in self.enzymes:
            max_n = max([max_n, enzyme["n"]])
        return max_n

    def _get_substrings(self):
        # Adjust the sequence to make sure the marker is in the center
        mseq = deepcopy(self.marker.seq)

        # Pad the sequence if needed to ensure the SNP is always in the center
        exp_seq_len = self.marker.HOMOLOG_FLANKING_N * 2 + 1
        act_seq_len = len(mseq)
        if act_seq_len < exp_seq_len:
            if int(self.marker.fasta_name.split(":")[1].split("-")[0]) == 1:
                # Left side is not full
                mseq = lib_utils.padding_left(
                    x=mseq,
                    n=exp_seq_len,
                )
            else:
                # Right side is not full
                mseq = lib_utils.padding_right(
                    x=mseq,
                    n=exp_seq_len,
                )

        self.marker_pos = int((len(mseq) - 1) / 2)  # SNP is now always in the center as we padded the sequence

        mseql = list(mseq)
        mseql[self.marker_pos] = self.marker.alt
        self.seqs = {
            "REF": mseq,
            "ALT": "".join(mseql),
        }

        self.substrings = {"REF": {}, "ALT": {}}
        for i in np.flip(range(1, self._get_recognition_site_max_length() + 1)):
            self.substrings["REF"][i] = self.seqs["REF"][(self.marker_pos - (i - 1)):(self.marker_pos + (i + 1) + 1)]
            self.substrings["ALT"][i] = self.seqs["ALT"][(self.marker_pos - (i - 1)):(self.marker_pos + (i + 1) + 1)]

    def find_valid_enzymes(self):
        self.valid_enzymes = []
        self.valid_enzymes_list = []

        for enzyme in self.enzymes:
            flag_ref = self.will_it_single_cut(self.substrings["REF"][enzyme["n"]], enzyme)
            flag_alt = self.will_it_cut(self.substrings["ALT"][enzyme["n"]], enzyme)

            if flag_ref and flag_alt:
                self.valid_enzymes.append(enzyme)
                self.valid_enzymes_list.append(enzyme["name"])

    @staticmethod
    def will_it_cut(seq, enzyme):
        m = re.findall(re.compile(enzyme["regex"]["F"]), seq)
        mi = re.findall(re.compile(enzyme["regex"]["R"]), seq)

        if (len(m) > 0) or (len(mi) > 0):
            return True
        else:
            return False

    @staticmethod
    def will_it_single_cut(seq, enzyme):
        m = re.findall(re.compile(enzyme["regex"]["F"]), seq)
        mi = re.findall(re.compile(enzyme["regex"]["R"]), seq)

        if (len(m) == 1) or (len(mi) == 1):
            return True
        else:
            return False


class PCR(_lib_pcr.PCR):
    def __init__(self, snp_id, chrom, pos, ref, alt):
        super().__init__(chrom=chrom)
        self.snp_id = snp_id
        self.ref = ref
        self.alt = alt
        self.chrom = chrom
        self.pos = pos


def print_report_header(fp, delimiter="\t"):
    header = delimiter.join([
        "marker",
        "chr",
        "pos",
        "ref",
        "alt",
        "enzymes",
        "band_left",
        "band_right"
        "start",
        "end",
        "direction",
        "assay_id",
        "seq_5_3",
        "primer_id",
        "goodness",
        "qcode",
        "length",
        "prod_size",
        "tm",
        "gc_content",
        "n_offtargets",
        "max_aaf",
        "indels",
        "offtargets",
        "mutations",
    ])

    with open(fp, "w") as f:
        f.write("{}\n".format(header))


def print_report(pcr, caps, fp, delimiter="\t"):
    primer_type_ordering = ["F", "R"]
    sorted_scores = np.sort(np.unique(list(pcr.pps_pruned.keys())))[::-1]
    ppid = 0

    seq_ids = {}
    for pt in primer_type_ordering:
        seq_ids[pt] = []

    enzyme_str = ",".join(caps.valid_enzymes_list)
    if enzyme_str == "":
        enzyme_str = "NA"

    with open(fp, "w") as f:
        for i in sorted_scores:
            for pp in pcr.pps_pruned[i]:
                pp.id = ppid

                seqs = {
                    "F": pp.primers["F"].sequence,
                    "R": pp.primers["R"].sequence,
                }

                curr_seq_ids = {}

                for d in primer_type_ordering:
                    if seqs[d] not in seq_ids[d]:
                        curr_seq_ids[d] = d + "_" + str(len(seq_ids[d]))
                        seq_ids[d].append(seqs[d])
                    else:
                        curr_seq_ids[d] = d + "_" + str(seq_ids[d].index(seqs[d]))

                # Shared properties
                offtargets = ",".join([x for x in pp.offtargets])
                if offtargets == "":
                    offtargets = "NA"

                # Write primers parameters
                for d in primer_type_ordering:
                    mutations = ",".join([x for x in pp.primers[d].mutations])

                    if mutations == "":
                        mutations = "NA"
                    if pp.qcode == "":
                        pp.qcode = "."

                    fields = [
                        pcr.snp_id,
                        pcr.chrom,
                        int(pcr.pos),
                        pcr.ref,
                        pcr.alt,
                        enzyme_str,
                        int(pcr.pos) - int(pp.primers["F"].start),
                        int(pp.primers["R"].stop) - int(pcr.pos),
                        int(pp.primers[d].start),
                        int(pp.primers[d].stop),
                        d,
                        pp.id,
                        seqs[d],
                        curr_seq_ids[d],
                        pp.goodness,
                        pp.qcode,
                        pp.primers[d].length,
                        pp.product_size,
                        lib_utils.round_tidy(pp.primers[d].tm, 3),
                        lib_utils.round_tidy(pp.primers[d].gc_content, 3),
                        len(pp.offtargets),
                        lib_utils.round_tidy(pp.primers[d].max_aaf, 2),
                        pp.max_indel_size,
                        offtargets,
                        mutations,
                    ]
                    f.write("{}\n".format(delimiter.join([str(x) for x in fields])))
                f.write("\n")
                ppid += 1

        if ppid == 0:
            fields = [
                pcr.snp_id,
                pcr.chrom,
                int(pcr.pos),
                pcr.ref,
                pcr.alt,
            ]
            fields = [str(x) for x in fields]
            f.write(delimiter.join(fields + 20 * ["NA"]) + "\n")
            f.write("\n")

        f.write("\n")


def main(kwarg_dict):
    # kwargs to variables
    fp_fasta = kwarg_dict["fp_fasta"]
    marker = kwarg_dict["marker"]
    fp_out = kwarg_dict["fp_out"]
    blast_db = kwarg_dict["blast_db"]
    muscle = kwarg_dict["muscle"]
    n_primers = kwarg_dict["n_primers"]
    p3_search_depth = kwarg_dict["p3_search_depth"]
    primer_seed = kwarg_dict["primer_seed"]
    tm_delta = kwarg_dict["tm_delta"]
    offtarget_size = kwarg_dict["offtarget_size"]
    primer3_configs = kwarg_dict["primer3_configs"]
    included_enzymes = kwarg_dict["included_enzymes"]
    debug = kwarg_dict["debug"]

    # Set primer3 globals
    if len(primer3_configs) > 0:
        lib_primer3.set_globals(**primer3_configs)

    lib_primer3.set_globals(
        PRIMER_NUM_RETURN=int(np.ceil(n_primers * p3_search_depth)),
    )
    lib_primer3.set_tm_delta(tm_delta)
    lib_primer3.set_primer_seed(primer_seed)

    # Set lib_primer offtarget sizes
    lib_primer3.set_offtarget_size(offtarget_size[0], offtarget_size[1])

    # Set a logger message that will be printed at the end (to be threadsafe)
    header = "Primer search results for {}".format(marker.name)
    sepline = "=" * (len(header) + 2)
    logger_msg = "\n{}\n{}\n{}\n".format(sepline, header, sepline)

    # Read sequences
    seqs = lib_blast.read_fasta(fp_fasta)
    marker.upload_sequence(seqs[marker.fasta_name])
    del seqs[marker.fasta_name]

    # Find all restriction enzymes that may work on the marker
    caps = CAPS(
        marker=marker,
        included_enzymes=included_enzymes,
    )

    caps.find_valid_enzymes()
    logger_msg += "Number of possible restriction enzymes: {}\n".format(caps.valid_enzymes)

    # Set target in PRIMER3
    lib_primer3.set_globals(
        SEQUENCE_TARGET=[max([1, caps.marker_pos - 49]), 101],
    )

    # Align homologs to the target sequence
    if len(seqs) > 49:
        # If more than 50 homoloous sequences were found, then do not run multialignment and
        # instead make all nucleotides not specific (using twice the same marker sequence) for runtime optimization
        mock_seqs = {
            marker.fasta_name: marker.seq,
            "mock": marker.seq,
        }
        lib_blast.write_fasta(mock_seqs, fp_out=fp_fasta)

    fp_aligned = join(blast_db.temporary, marker.name + ".afa")
    muscle.align_fasta(fp=fp_fasta, fp_out=fp_aligned)

    # Map mismatches in the homeologs
    maps = dict()
    maps["all"] = np.repeat(True, marker.n)
    maps["partial_mism"], maps["mism"] = _lib_kasp.map_homologs(fp_aligned, marker.fasta_name, marker.n)

    # Cleanup temporary files
    if not debug:
        os.remove(fp_fasta)
        os.remove(fp_aligned)

    # Get valid regions for the design of the primers
    ivs = {}
    for k, mmap in maps.items():
        ivs[k] = {}
        k_mut = k + "_mut"
        ivs[k] = _lib_pcr.get_valid_primer_regions(mmap)

        if len(marker.mutations) > 0:
            ivs[k_mut] = _lib_pcr.get_valid_primer_regions(mmap,
                                                           hard_exclude=marker.mutations)

    # Design primers
    pcr = PCR(marker.name, marker.chrom, marker.pos1, marker.ref, marker.alt)

    # Only design primers if a restriction enzyme exists
    if len(caps.valid_enzymes) > 0:
        map_names = {
            "mism_mut": "Variants excluded | Homologs specific",
            "partial_mism_mut": "Variants excluded | Homologs partial spec",
            "all_mut": "Variants excluded | Unspecific",
            "mism": "Variants included | Homologs specific",
            "partial_mism": "Variants included | Homologs partial spec",
            "all": "Variants included | Unspecific",
        }
        # Set the search mask ordering
        if len(marker.mutations) > 0:
            search_types = ["mism_mut", "mism", "partial_mism_mut", "partial_mism", "all_mut", "all"]
        else:
            search_types = ["mism", "partial_mism", "all"]

        for search_type in search_types:
            if len(pcr.pps) < lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN']:
                n_before = len(pcr.pps)
                pcr.pps = _lib_pcr.design_primers(
                    pps_repo=pcr.pps,  # Primer pair repository
                    target_seq=marker.seq,
                    target_chrom=marker.chrom,
                    target_start=marker.start,
                    ivs=ivs[search_type],
                    n_primers=lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN'],
                )
                n_new = len(pcr.pps) - n_before
                logger_msg += "{:42}: {:3d} pairs\n".format(map_names[search_type], n_new)

        hit_cnts, valid_hit_cnts = pcr.check_offtargeting(blast_db, debug=debug)  # Check for off-targets

        # Report hit counts
        logger_msg += "Offtarget hits across the genome\n"
        logger_msg += "      Forward : {:d}\n".format(valid_hit_cnts["F"])
        logger_msg += "      Reverse : {:d}\n".format(valid_hit_cnts["R"])

    pcr.add_mutations(marker.mutations)  # List mutations in primers
    pcr.classify()  # Classify primers by scores using a heuristic "goodness" score
    n = pcr.prune(n_primers)  # Retain only the top n primers

    # Print to logger
    logger_msg += "Returned top {:d} primer pairs\n".format(n)
    logger.debug(logger_msg)
    print_report(pcr, caps, fp_out)


if __name__ == "__main__":
    pass
