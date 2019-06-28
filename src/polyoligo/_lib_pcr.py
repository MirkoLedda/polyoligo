from __future__ import print_function, division
# noinspection PyPackageRequirements
from Bio.Seq import Seq
import numpy as np
# noinspection PyUnresolvedReferences
from primer3.thermoanalysis import ThermoAnalysis
from copy import deepcopy
import logging
from os.path import join
import os
import sys
import pandas as pd

from . import lib_blast, lib_utils, lib_primer3

# Initialize the logger
logger = logging.getLogger(__name__)


class ROI:
    def __init__(self, roi, blast_db):
        self.chrom = roi.strip().split(":")[0]
        self.start = int(roi.strip().split(":")[1].split("-")[0])
        self.end = int(roi.strip().split(":")[1].split("-")[1])
        self.blast_db = blast_db

        self.seq = None
        self.n = None
        self.G = None  # Genotypes within the region
        self.mutations = []  # List of vcf_lib.Mutation objects
        self.fasta_names = {"winleft": None, "winright": None}
        self.win_seqs = {"winleft": None, "winright": None}

    def fetch_roi(self):
        query = [{
            "chr": self.chrom,
            "start": self.start,
            "stop": self.end,
        }]

        seqs = self.blast_db.fetch(query)
        self.seq = list(seqs.values())[0].upper()
        self.n = len(self.seq)

    def upload_mutations(self, vcf_obj, start, stop):
        if vcf_obj:
            self.G, self.mutations = vcf_obj.fetch_genotypes(
                chrom=self.chrom,
                start=start,
                stop=stop,
            )
        else:
            self.mutations = []

    def print_alt_subjects(self, vcf_obj, fp):
        if vcf_obj:
            vcf_obj.print_alternative_subjects(
                fp=fp,
                chrom=self.chrom,
                start=self.start,
                stop=self.end,
            )


class PCR(lib_primer3.PCR):
    def __init__(self, chrom):
        super().__init__([])
        self.pps_classified = {}
        self.pps_pruned = {}
        self.chrom = chrom

    def add_mutations(self, mutations):
        if len(mutations) > 0:
            for pp in self.pps:
                pp.add_mutations(mutations)

    def classify(self):
        for pp in self.pps:
            pp.score()
            if pp.goodness not in self.pps_classified.keys():
                self.pps_classified[pp.goodness] = []

            self.pps_classified[pp.goodness].append(pp)

    def prune(self, n):
        nprim = 0
        score_cats = np.sort(list(self.pps_classified.keys()))[::-1]

        for score in score_cats:
            self.pps_pruned[score] = []
            for pp in self.pps_classified[score]:
                if nprim < n:
                    self.pps_pruned[score].append(pp)
                    nprim += 1
                else:
                    break

        return nprim


def map_homologs(fp_aligned, target_name, target_len):
    seqs = lib_blast.read_fasta(fp_aligned)
    target_seq = seqs[target_name]
    del seqs[target_name]

    n_homeo = len(seqs)
    match_arr = np.tile(np.nan, (target_len, n_homeo))

    # Count mismatches
    homeo_id = 0
    for chrom, seq in seqs.items():
        j = 0  # Index in the original sequence without gaps
        for i, nuc in enumerate(target_seq):
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
    partial_mismatch = match_cnts != n_homeo
    full_mismatch = match_cnts == 0

    return partial_mismatch, full_mismatch


def print_report_header(fp, delimiter="\t"):
    header = delimiter.join([
        "chr",
        "start",
        "end",
        "direction",
        "assay_id",
        "seq_5_3",
        "seq_5_3_ambiguous",
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


def print_report(pcr, fp, delimiter="\t"):
    primer_type_ordering = ["F", "R"]
    sorted_scores = np.sort(np.unique(list(pcr.pps_pruned.keys())))[::-1]
    ppid = 0

    seq_ids = {}
    for pt in primer_type_ordering:
        seq_ids[pt] = []

    with open(fp, "a") as f:
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
                        pcr.chrom,
                        int(pp.primers[d].start),
                        int(pp.primers[d].stop),
                        d,
                        pp.id,
                        seqs[d],
                        pp.primers[d].sequence_ambiguous,
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
                pcr.chrom,
            ]
            fields = [str(x) for x in fields]
            f.write(delimiter.join(fields + 18 * ["NA"]) + "\n")
            f.write("\n")

        f.write("\n")


def design_primers(pps_repo, target_seq, target_chrom, target_start, ivs, n_primers=10):
    primer3_seq_args = {
        'SEQUENCE_TEMPLATE': target_seq,
        'SEQUENCE_EXCLUDED_REGION': ivs,
    }

    p3_repo = lib_primer3.get_primers(primer3_seq_args, target_start=target_start)

    # Find valid primer pairs by checking top hits for each marker primers iteratively
    flag_continue = True
    while flag_continue:

        if len(p3_repo) == 0:
            flag_continue = False

        for pp in p3_repo:

            if len(pps_repo) == n_primers:  # We have enough primers
                flag_continue = False
                break

            pp.chrom = target_chrom
            is_pp_valid = True
            for d in pp.primers.keys():
                is_valid = pp.primers[d].is_valid()
                if not is_valid:
                    is_pp_valid = False

            if is_pp_valid:
                pps_repo.append(pp)

    return pps_repo


def main(kwarg_dict):
    # kwargs to variables
    roi = kwarg_dict["roi"]
    blast_db = kwarg_dict["blast_db"]
    muscle = kwarg_dict["muscle"]
    n_primers = kwarg_dict["n_primers"]
    p3_search_depth = kwarg_dict["p3_search_depth"]
    primer_seed = kwarg_dict["primer_seed"]
    tm_delta = kwarg_dict["tm_delta"]
    offtarget_size = kwarg_dict["offtarget_size"]
    primer3_configs = kwarg_dict["primer3_configs"]
    debug = kwarg_dict["debug"]
    fp_out = kwarg_dict["fp_out"]

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
    header = "Primer search results"
    sepline = "=" * (len(header) + 2)
    logger_msg = "\n{}\n{}\n{}\n".format(sepline, header, sepline)

    # Read and align sequences for both the left and right window
    win_maps = {"winleft": {}, "winright": {}}
    for win_name in ["winleft", "winright"]:
        fp_fasta = join(blast_db.temporary, win_name + ".fa")
        seqs = lib_blast.read_fasta(fp_fasta)

        roi.win_seqs[win_name] = seqs[roi.fasta_names[win_name]]
        del seqs[roi.fasta_names[win_name]]

        # Align homologs to the target sequence
        if len(seqs) > 49:
            # If more than 50 homoloous sequences were found, then do not run multialignment and
            # instead make all nucleotides not specific (using twice the same marker sequence) for runtime optimization
            mock_seqs = {
                roi.fasta_name: roi.seq,
                "mock": roi.seq,
            }
            lib_blast.write_fasta(mock_seqs, fp_out=fp_fasta)

        fp_aligned = join(blast_db.temporary, win_name + "_malign.afa")
        muscle.fast_align_fasta(fp=fp_fasta, fp_out=fp_aligned)

        # Map mismatches in the homeologs
        win_maps[win_name]["all"] = np.repeat(True, len(roi.win_seqs[win_name]))
        win_maps[win_name]["partial_mism"], win_maps[win_name]["mism"] = map_homologs(
            fp_aligned,
            roi.fasta_names[win_name],
            len(roi.win_seqs[win_name]))

        # Cleanup temporary files
        if not debug:
            os.remove(fp_fasta)
            os.remove(fp_aligned)

    # Merge all maps
    maps = dict()
    maps["all"] = np.concatenate(
        [win_maps["winleft"]["all"][:-1], np.repeat(False, roi.n), win_maps["winright"]["all"][1:]])
    maps["partial_mism"] = np.concatenate(
        [win_maps["winleft"]["partial_mism"][:-1], np.repeat(False, roi.n), win_maps["winright"]["partial_mism"][1:]])
    maps["mism"] = np.concatenate(
        [win_maps["winleft"]["mism"][:-1], np.repeat(False, roi.n), win_maps["winright"]["mism"][1:]])

    # Get valid regions for the design of the primers
    ivs = {}
    roi.start = int(roi.fasta_names["winleft"].split(":")[1].split("-")[0])
    roi.end = int(roi.fasta_names["winright"].split(":")[1].split("-")[1])
    for k, mmap in maps.items():
        ivs[k] = {}
        k_mut = k + "_mut"
        ivs[k] = lib_primer3.get_exclusion_zone(mmap)  # Mutations not excluded

        if len(roi.mutations) > 0:
            # Make a list of mutation positions based on mutation for exclusion
            mut_ixs = []
            for mutation in roi.mutations:
                mut_ixs.append(mutation.pos - roi.start)
            ivs[k_mut] = lib_primer3.get_exclusion_zone(mmap, hard_exclude=mut_ixs)

    # Set the product size range
    lib_primer3.set_globals(
        PRIMER_PRODUCT_SIZE_RANGE=[roi.n, len(maps["all"])]
    )

    # Update roi to be the entire sequence including the side windows
    roi.fetch_roi()

    # Design primers
    pcr = PCR(roi.chrom)
    map_names = {
        "mism_mut": "Variants excluded | Homologs specific",
        "partial_mism_mut": "Variants excluded | Homologs partial spec",
        "all_mut": "Variants excluded | Unspecific",
        "mism": "Variants included | Homologs specific",
        "partial_mism": "Variants included | Homologs partial spec",
        "all": "Variants included | Unspecific",
    }
    # Set the search mask ordering
    if len(roi.mutations) > 0:
        search_types = ["mism_mut", "mism", "partial_mism_mut", "partial_mism", "all_mut", "all"]
    else:
        search_types = ["mism", "partial_mism", "all"]

    print(ivs["all"])

    for search_type in search_types:
        if len(pcr.pps) < lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN']:
            n_before = len(pcr.pps)
            pcr.pps = design_primers(
                pps_repo=pcr.pps,  # Primer pair repository
                target_seq=roi.seq,
                target_chrom=roi.chrom,
                target_start=roi.start,
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

    pcr.add_mutations(roi.mutations)  # List mutations in primers
    pcr.classify()  # Classify primers by scores using a heuristic "goodness" score
    n = pcr.prune(n_primers)  # Retain only the top n primers

    # Print to logger
    logger_msg += "Returned top {:d} primer pairs\n".format(n)
    logger.debug(logger_msg)
    print_report_header(fp_out)
    print_report(pcr, fp_out)


if __name__ == "__main__":
    pass
