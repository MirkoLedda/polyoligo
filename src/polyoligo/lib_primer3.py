import numpy as np
# noinspection PyUnresolvedReferences
from primer3.thermoanalysis import ThermoAnalysis
from copy import deepcopy
# noinspection PyPackageRequirements
from Bio import SeqUtils
import primer3
import sys
from os.path import join
import os

from . import lib_blast

# Default Primer3 globals
PRIMER3_GLOBALS = {
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 55.0,
    'PRIMER_MAX_TM': 65.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE': 1,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_DIFF_TM': 6.0,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PRODUCT_SIZE_RANGE': [[90, 150], [100, 200], [50, 250]],
    'PRIMER_FIRST_BASE_INDEX': 0,
}

PRIMER_SEED = 12
MIN_OFFTARGET_SIZE = 0  # Minimum size of a PCR product in the offtargets
MAX_OFFTARGET_SIZE = 1000  # Maximum size of a PCR product in the offtargets
TM_DELTA = 5  # Minimum difference in melting temp for unwanted structures in the primer


class Primer:
    def __init__(self, **kwargs):
        self.id = None  # Primer ID
        self.is_marker = None  # If this primer contains the target SNP
        self.marker_pos = None  # Position of the marker in the primer
        self.start = None  # Primer start position
        self.stop = None  # Primer stop position
        self.length = None  # Primer length
        self.direction = None  # Primer direction
        self.hairpin_tm = None  # Hairpin melting temp
        self.homodimer_tm = None  # Homodimer melting temp
        self.gc_content = None  # % GC content
        self.has_thermals = False  # If the primers has thermal properties
        self.tm_delta = TM_DELTA  # Minimum difference in melting temp for unwanted structures in the primer
        self.tm = None  # tm from PRIMER3 is excluded from automated import
        self.mutations = []  # List of mutations within the primer
        self.max_aaf = 0  # Max alternative allele frequency for mutations within the primer
        self.nN = 0  # Number of N's in the primer sequence
        self.seed = None  # Sequence of the seed region

        # Attributes shared with PRIMER3
        self.sequence = None
        self.end_stability = None
        self.hairpin_th = None
        self.penalty = None
        self.self_any_th = None
        self.self_end_th = None

        # Set attributes from kwargs
        kwargs2attr_deep(self, kwargs)

    def calc_thermals(self):
        if not self.has_thermals:
            thals_obj = ThermoAnalysis()  # Initialize a primer3 thal object
            self.tm = thals_obj.calcTm(self.sequence)
            self.hairpin_tm = thals_obj.calcHairpin(self.sequence).tm
            self.homodimer_tm = thals_obj.calcHomodimer(self.sequence).tm
            self.length = len(self.sequence)
            self.gc_content = SeqUtils.GC(self.sequence)
            self.nN = self.sequence.count("N")

    def is_valid(self):
        if not self.has_thermals:
            self.calc_thermals()

        passed_checks = True
        # Check for N in the sequence
        if self.nN > 0:
            passed_checks = False

        # Check length
        if (self.length < PRIMER3_GLOBALS["PRIMER_MIN_SIZE"]) or (self.length > PRIMER3_GLOBALS["PRIMER_MAX_SIZE"]):
            passed_checks = False

        # Check GC content
        if (self.gc_content < PRIMER3_GLOBALS["PRIMER_MIN_GC"]) or (
                self.gc_content > PRIMER3_GLOBALS["PRIMER_MAX_GC"]):
            passed_checks = False

        # Check tm
        if (self.tm < PRIMER3_GLOBALS["PRIMER_MIN_TM"]) or (self.tm > PRIMER3_GLOBALS["PRIMER_MAX_TM"]):
            passed_checks = False

        # Check hairpins
        if self.hairpin_tm >= (self.tm - self.tm_delta):
            passed_checks = False

        # Check homodimerization
        if self.homodimer_tm >= (self.tm - self.tm_delta):
            passed_checks = False

        return passed_checks

    def add_primer3_attributes(self, p3):
        self.id = p3.id
        self.start = p3.start
        self.stop = p3.stop
        self.sequence = p3.sequence
        self.end_stability = p3.end_stability
        self.hairpin_th = p3.hairpin_th
        self.penalty = p3.penalty
        self.self_any_th = p3.self_any_th
        self.self_end_th = p3.self_end_th


class PrimerPair:
    def __init__(self):
        self.tm_delta = TM_DELTA  # Minimum difference in melting temp for unwanted structures in the primer
        self.id = None  # ID of the primer pair
        self.primers = {
            "F": Primer(),  # Forward primer
            "R": Primer(),  # Reverse primer
        }
        self.offtargets = []  # List of offtarget sites
        self.max_indel_size = 0  # Maximum indel size within the PCR product
        self.chrom = None
        self.goodness = None  # Heuristic goodness score
        self.qcode = None  # Heuristic quality code, gives context to the goodness score

        # Attributes from PRIMER3
        self.compl_any_th = None
        self.compl_end_th = None
        self.product_size = None
        self.penalty = None

    def add_primer3_attributes(self, pp3):
        for d in ["F", "R"]:
            self.primers[d].add_primer3_attributes(pp3.primers[d])
        self.compl_any_th = pp3.compl_any_th
        self.compl_end_th = pp3.compl_end_th
        self.product_size = pp3.product_size
        self.penalty = pp3.penalty

    def add_mutations(self, mutations):

        # Check each primer individually
        for d in self.primers.keys():
            window = np.arange(self.primers[d].start, self.primers[d].stop + 1)

            for m in mutations:
                if m.pos in window:
                    mut_txt = "{}{:d}{}:{}".format(m.ref, m.pos, "/".join(m.alt), round_tidy(m.aaf, 2))
                    self.primers[d].mutations.append(mut_txt)
                    self.primers[d].max_aaf = np.max([self.primers[d].max_aaf, m.aaf])

        # Check indels in the PCR product
        window = np.arange(self.primers["F"].start, self.primers["R"].stop + 1)

        for m in mutations:
            if m.pos in window:
                self.max_indel_size = np.max([self.max_indel_size, m.indel_size])

    def score(self):
        score = 0
        qcode = ""

        tms = np.array([self.primers[d].tm for d in self.primers.keys()])
        tms_l1_norm = np.sum(np.abs(tms - np.mean(tms)))
        if tms_l1_norm <= 3:
            score += 1
        else:
            qcode += "t"

        if len(self.offtargets) == 0:
            score += 3
        else:
            qcode += "O"

        # todo heterodimers
        # if (not self.ref_dimer) and (not self.alt_dimer):
        #     score += 1
        # else:
        #     qcode += "d"
        score += 1

        max_aafs = np.array([self.primers[d].max_aaf for d in self.primers.keys()])
        if np.all(max_aafs < 0.1):
            score += 2

            if np.all(max_aafs == 0):
                score += 1
            else:
                qcode += "m"

        else:
            qcode += "M"

        if self.max_indel_size < 50:
            score += 1

            if self.max_indel_size == 0:
                score += 1
            else:
                qcode += "i"

        else:
            qcode += "I"

        self.goodness = score
        self.qcode = qcode


class PCR:
    def __init__(self, primer_pairs=None):
        self.pps = primer_pairs

        if self.pps is None:
            self.pps = []

    def get_seeds(self):
        for pp in self.pps:

            for k in pp.primers.keys():
                pp.primers[k].seed = pp.primers[k].sequence[(len(pp.primers[k].sequence) - PRIMER_SEED):]

    def get_unique_seed_sequences(self):
        self.get_seeds()
        seqs_dict = {
            "F": [],
            "R": [],
        }

        # List all sequences
        for pp in self.pps:
            seqs_dict["F"].append(pp.primers["F"].seed)
            seqs_dict["R"].append(pp.primers["R"].seed)

        # Get unique sequences
        for d, seqs in seqs_dict.items():
            seqs_dict[d] = list(np.unique(seqs))

        return seqs_dict

    # noinspection PyPep8Naming
    def check_offtargeting(self, blast_db, debug=False):

        # Init counters for logging
        hits_cnts = {"F": 0, "R": 0}
        valid_hits_cnts = {"F": 0, "R": 0}

        if self.pps:
            pseqs = self.get_unique_seed_sequences()

            hits = {}

            # Lookup forward primers
            fp_query_F = join(blast_db.temporary, "{}_blast_F.fa".format(blast_db.job_id))
            fp_out_F = join(blast_db.temporary, "{}_blast_F.json".format(blast_db.job_id))
            lib_blast.write_fasta(pseqs["F"], fp_query_F)
            blast_db.blastn(
                fp_query=fp_query_F,
                fp_out=fp_out_F,
                word_size=np.min([PRIMER3_GLOBALS["PRIMER_MIN_SIZE"], PRIMER_SEED]),
                strand="plus",
                evalue=10000,
                max_hsps=10000,
                perc_identity=100,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
                ungapped="",
            )

            hits["F"] = blast_db.parse_blastn_json(fp_out_F, min_identity=PRIMER_SEED)

            # Lookup reverse primers
            fp_query_R = join(blast_db.temporary, "{}_blast_R.fa".format(blast_db.job_id))
            fp_out_R = join(blast_db.temporary, "{}_blast_R.json".format(blast_db.job_id))
            lib_blast.write_fasta(pseqs["R"], fp_query_R)
            blast_db.blastn(
                fp_query=fp_query_R,
                fp_out=fp_out_R,
                word_size=np.min([PRIMER3_GLOBALS["PRIMER_MIN_SIZE"], PRIMER_SEED]),
                strand="minus",
                evalue=10000,
                max_hsps=10000,
                perc_identity=100,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
                ungapped="",
            )
            hits["R"] = blast_db.parse_blastn_json(fp_out_R, min_identity=PRIMER_SEED)

            # Check if primers could bind
            valid_hits = {}
            for d in ["F", "R"]:
                valid_hits[d] = {}

                for k, hit in hits[d].items():
                    pseq = pseqs[d][int(k)]  # Sequence of the primer/blast query
                    valid_hits[d][pseq] = {}

                    for chrom, shit in hit.items():
                        if chrom not in valid_hits[d][pseq].keys():
                            valid_hits[d][pseq][chrom] = []

                        for sshit in shit.values():
                            hits_cnts[d] += 1
                            if will_it_bind(sshit["midline"]):
                                valid_hits[d][pseq][chrom].append(sshit)
                                valid_hits_cnts[d] += 1

            # Check for each primer pair if PCR products are made at offtarget sites
            for pp in self.pps:
                pp_hits = {
                    "F": deepcopy(valid_hits["F"][pp.primers["F"].seed]),
                    "R": deepcopy(valid_hits["R"][pp.primers["R"].seed]),
                }

                offtargets = self.list_offtargets(pp_hits)

                # Remove the entry for the on target
                left_padding = len(pp.primers["F"].sequence) - len(pp.primers["F"].seed)
                right_padding = len(pp.primers["R"].sequence) - len(pp.primers["R"].seed)

                ontarget = "{}:{}-{}".format(
                    pp.chrom,
                    pp.primers["F"].start + left_padding,
                    pp.primers["R"].stop - right_padding,
                )
                if ontarget in offtargets:
                    _ = offtargets.pop(offtargets.index(ontarget))
                else:
                    sys.exit("ERROR: Oups, something went wrong. Please contact the developer.")

                # Rename offtargets to adjust for seed padding
                offtargets_parsed = []
                for offt in offtargets:
                    fields = offt.split(":")
                    loc = fields[1].split("-")
                    offtargets_parsed.append("{}:{}-{}".format(
                        fields[0],
                        int(loc[0]) - left_padding,
                        int(loc[1]) + right_padding,
                    ))

                pp.offtargets = offtargets_parsed

            # Cleanup temporary files
            if not debug:
                os.remove(fp_query_F)
                os.remove(fp_query_R)
                os.remove(fp_out_F)
                os.remove(fp_out_R)

        return hits_cnts, valid_hits_cnts

    @staticmethod
    def list_offtargets(pp_hits):
        offtargets = []
        for chrom in pp_hits["F"].keys():
            if chrom in pp_hits["R"].keys():
                starts = []
                stops = []

                for hit in pp_hits["F"][chrom]:
                    starts.append(hit["sstart"])

                for hit in pp_hits["R"][chrom]:
                    stops.append(hit["sstart"])

                starts = np.sort(np.unique(starts))
                stops = np.sort(np.unique(stops))

                for i in starts:
                    for j in stops:
                        if j <= i:
                            pass
                        else:
                            delta = j-i
                            if (delta > MIN_OFFTARGET_SIZE) and (delta < MAX_OFFTARGET_SIZE):
                                offtargets.append("{}:{}-{}".format(chrom, i, j))
                            else:
                                break
        return offtargets


def get_primers(primer3_seq_args, target_start=1):
    pps = {}

    if "SEQUENCE_PRIMER" in PRIMER3_GLOBALS:
        primer3_seq_args["SEQUENCE_PRIMER"] = PRIMER3_GLOBALS["SEQUENCE_PRIMER"]

    if "SEQUENCE_PRIMER_REVCOMP" in PRIMER3_GLOBALS:
        primer3_seq_args["SEQUENCE_PRIMER_REVCOMP"] = PRIMER3_GLOBALS["SEQUENCE_PRIMER_REVCOMP"]

    p3_primers = primer3.bindings.designPrimers(primer3_seq_args)

    # Note that primers are returned ordered by "quality" (i.e. lower penalties first)
    for k, v in p3_primers.items():
        if (k == "PRIMER_WARNING") or (k == "PRIMER_ERROR"):
            pps = {}  # Remove all stored primers for this design run as there was an issue
            break

        fields = k.split("_")

        try:
            pid = int(fields[2])
        except ValueError:
            continue

        field_name = "_".join(fields[3:])
        field_name = field_name.lower()

        if field_name not in ["tm"]:  # Exclude fields from Primer3
            # Check if we initialized objects for this primer ID
            if pid not in pps.keys():
                pps[pid] = PrimerPair()
                for d in ["F", "R"]:
                    pps[pid].primers[d].id = pid

            # Check if the curr parameter relates to a primer or a pair and assign it to the corresponding attribute
            if fields[1] == "LEFT":
                d = "F"
                if field_name:
                    pps[pid].primers[d].__dict__[field_name] = v
                else:
                    pps[pid].primers[d].start = target_start + v[0]
                    pps[pid].primers[d].stop = pps[pid].primers[d].start + (v[1] - 1)

            elif fields[1] == "RIGHT":
                d = "R"
                if field_name:
                    pps[pid].primers[d].__dict__[field_name] = v
                else:
                    pps[pid].primers[d].stop = target_start + v[0]
                    pps[pid].primers[d].start = pps[pid].primers[d].stop - (v[1] - 1)

            elif fields[1] == "PAIR":
                pps[pid].__dict__[field_name] = v

    # dict to ordered list
    pps_list = []
    for i in np.sort(list(pps.keys())):
        pps_list.append(pps[i])

    return pps_list


def will_it_bind(midline):
    if "X" in midline:
        return False
    else:
        return True


def set_globals(**kwargs):
    global PRIMER3_GLOBALS
    for k, v in kwargs.items():
        PRIMER3_GLOBALS[k] = v

    primer3.setP3Globals(PRIMER3_GLOBALS)


def set_offtarget_size(min_size, max_size):
    global MIN_OFFTARGET_SIZE
    MIN_OFFTARGET_SIZE = min_size
    global MAX_OFFTARGET_SIZE
    MAX_OFFTARGET_SIZE = max_size


def set_tm_delta(v):
    global TM_DELTA
    TM_DELTA = v


def set_primer_seed(v):
    global PRIMER_SEED
    PRIMER_SEED = v


def kwargs2attr_deep(obj, kwargs):
    """Deep copy attributes values based on keyword arguments. Assign only if the attribute already exists in obj."""
    if kwargs is None:
        return

    if not isinstance(kwargs, dict):
        raise TypeError("Attempted to set arguments not using a dictionary.")

    attr = obj.__dict__.keys()
    for key in kwargs.keys():
        if key in attr:
            setattr(obj, key, deepcopy(kwargs[key]))


def round_tidy(x, n):
    """Return 'x' rounded to 'n' significant digits."""
    y = np.abs(x)
    if y <= sys.float_info.min:
        return 0.0
    else:
        return np.round(x, int(n - np.ceil(np.log10(y))))
