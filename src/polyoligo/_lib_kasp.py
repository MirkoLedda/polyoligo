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

from . import lib_blast, lib_utils, lib_primer3

# Initialize the logger
logger = logging.getLogger(__name__)


class Primer(lib_primer3.Primer):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.reporter = ""  # Sequence of the reporter dye


class PrimerPair(lib_primer3.PrimerPair):
    def __init__(self):
        super().__init__()
        self.snp_id = None
        self.chr = None
        self.primers = {
            "F": lib_primer3.Primer(),  # Forward primer
            "R": lib_primer3.Primer(),  # Reverse primer
            "A": lib_primer3.Primer(),  # Alternative primer
        }
        self.dir = None  # The direction of the primer where the marker is located (F or R)
        self.ref_dimer = None  # Dimerization for the primer with the ref allele and the reverse primer
        self.alt_dimer = None  # Dimerization for the primer with the alt allele and the reverse primer
        self.goodness = None  # Heuristic goodness score
        self.qcode = None  # Heuristic quality code, gives context to the goodness score

    # noinspection PyUnresolvedReferences
    def add_reporter_dyes(self, reporters):

        if self.dir == "F":
            self.primers["F"].reporter = reporters[0]
            self.primers["A"].reporter = reporters[1]
            self.primers["R"].reporter = ""
        else:
            self.primers["F"].reporter = ""
            self.primers["R"].reporter = reporters[0]
            self.primers["A"].reporter = reporters[1]
            # self.primers["A"].reporter = str(Seq(reporters[1]).reverse_complement())
            # self.primers["R"].reporter = str(Seq(reporters[0]).reverse_complement())

    # noinspection PyUnresolvedReferences
    def check_heterodimerization(self):

        for d in self.primers.keys():
            if not self.primers[d].has_thermals:
                self.primers[d].calc_thermals()

        seqs = {}
        for d in self.primers.keys():
            seqs[d] = self.primers[d].reporter + self.primers[d].sequence

        thals_obj = ThermoAnalysis()  # Initialize a primer3 thal object
        tm_ref = thals_obj.calcHeterodimer(seqs["F"], seqs["R"]).tm - self.tm_delta

        if self.dir == "F":
            tm_alt = thals_obj.calcHeterodimer(seqs["A"], seqs["R"]).tm - self.tm_delta
        else:
            tm_alt = thals_obj.calcHeterodimer(seqs["F"], seqs["A"]).tm - self.tm_delta

        # Primer with the reference allele
        if (tm_ref <= (self.primers["F"].tm - self.tm_delta)) and (tm_ref <= (self.primers["R"].tm - self.tm_delta)):
            self.ref_dimer = False
        else:
            self.ref_dimer = True

        # Primer with the alternative allele
        if (tm_alt <= (self.primers["F"].tm - self.tm_delta)) and (tm_alt <= (self.primers["R"].tm - self.tm_delta)):
            self.alt_dimer = False
        else:
            self.alt_dimer = True

    def score(self):
        score = 0
        qcode = ""

        tms = np.array([self.primers[d].tm for d in self.primers.keys()])
        tms_l1_norm = np.sum(np.abs(tms - np.mean(tms)))
        if tms_l1_norm <= 5:
            score += 1
        else:
            qcode += "t"

        if len(self.offtargets) == 0:
            score += 3
        else:
            qcode += "O"

        if (not self.ref_dimer) and (not self.alt_dimer):
            score += 1
        else:
            qcode += "d"

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


# Inherited
# noinspection PyPep8Naming
class PCR(lib_primer3.PCR):
    def __init__(self, snp_id, chrom, pos, ref, alt):
        super().__init__([])
        self.pps_classified = {}
        self.pps_pruned = {}
        self.snp_id = snp_id
        self.ref = ref
        self.alt = alt
        self.chrom = chrom
        self.pos = pos

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
            seqs_dict[pp.dir].append(pp.primers["A"].seed)

        # Get unique sequences
        for d, seqs in seqs_dict.items():
            seqs_dict[d] = list(np.unique(seqs))

        return seqs_dict

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
                word_size=np.min([lib_primer3.PRIMER3_GLOBALS["PRIMER_MIN_SIZE"], lib_primer3.PRIMER_SEED]),
                strand="plus",
                evalue=10000,
                perc_identity=100,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
                ungapped="",
            )
            hits["F"] = blast_db.parse_blastn_json(fp_out_F, min_identity=lib_primer3.PRIMER_SEED)

            # Lookup reverse primers
            fp_query_R = join(blast_db.temporary, "{}_blast_R.fa".format(blast_db.job_id))
            fp_out_R = join(blast_db.temporary, "{}_blast_R.json".format(blast_db.job_id))
            lib_blast.write_fasta(pseqs["R"], fp_query_R)
            blast_db.blastn(
                fp_query=fp_query_R,
                fp_out=fp_out_R,
                word_size=np.min([lib_primer3.PRIMER3_GLOBALS["PRIMER_MIN_SIZE"], lib_primer3.PRIMER_SEED]),
                strand="minus",
                evalue=10000,
                perc_identity=100,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
                ungapped="",
            )
            hits["R"] = blast_db.parse_blastn_json(fp_out_R, min_identity=lib_primer3.PRIMER_SEED)

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
                            if lib_primer3.will_it_bind(sshit["midline"]):
                                valid_hits[d][pseq][chrom].append(sshit)
                                valid_hits_cnts[d] += 1

            # Check for each primer pair if PCR products are made at offtarget sites
            for pp in self.pps:
                pp_hits = {
                    "F": deepcopy(valid_hits["F"][pp.primers["F"].seed]),
                    "R": deepcopy(valid_hits["R"][pp.primers["R"].seed]),
                }

                # Merge the hits from the alternative marker allele
                for chrom, ahits in valid_hits[pp.dir][pp.primers["A"].seed].items():
                    if chrom in pp_hits[pp.dir].keys():
                        pp_hits[pp.dir][chrom] += ahits
                    else:
                        pp_hits[pp.dir][chrom] = ahits

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
                    pass

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

    def add_reporter_dyes(self, reporters):
        for pp in self.pps:
            pp.add_reporter_dyes(reporters)

    def check_heterodimerization(self):
        for pp in self.pps:
            pp.check_heterodimerization()

    def add_mutations(self, mutations):
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


# noinspection PyDefaultArgument
def get_marker_primers(marker, allowed_end_pos=[0]):
    mseq = deepcopy(marker.seq)

    # Pad the sequence if needed to ensure the SNP is always in the center
    exp_seq_len = marker.HOMOLOG_FLANKING_N * 2 + 1
    act_seq_len = len(mseq)
    if act_seq_len < exp_seq_len:
        if int(marker.fasta_name.split(":")[1].split("-")[0]) == 1:
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

    marker_pos = int((len(mseq) - 1) / 2)  # SNP is always in the center as we padded the sequence

    seq = Seq(mseq)
    seq_rc = seq.reverse_complement()
    alt_a = Seq(marker.alt)
    alt_a_rc = alt_a.reverse_complement()

    pps = []
    for padding in range(lib_primer3.PRIMER3_GLOBALS["PRIMER_MIN_SIZE"],
                         lib_primer3.PRIMER3_GLOBALS["PRIMER_MAX_SIZE"] + 1):
        for i in allowed_end_pos:
            x2 = marker_pos + i
            x1 = x2 - padding

            for direction in ["F", "R"]:
                if direction == "F":
                    pseq = str(seq[x1:(x2 + 1)])
                    pos_a = len(pseq) - 1 - i  # Position of the allele in the primer
                    pseq_alt = list(pseq)
                    pseq_alt[pos_a] = alt_a[0]
                    pseq_alt = "".join(pseq_alt)
                else:
                    pseq = str(seq_rc[x1:(x2 + 1)])
                    pos_a = len(pseq) - 1 - i  # Position of the allele in the primer
                    pseq_alt = list(pseq)
                    pseq_alt[pos_a] = alt_a_rc[0]
                    pseq_alt = "".join(pseq_alt)

                ref_p = lib_primer3.Primer(sequence=pseq)
                alt_p = lib_primer3.Primer(sequence=pseq_alt)

                # Compute primer properties and thermals
                ref_p.calc_thermals()
                alt_p.calc_thermals()

                # Series of checks to determine if the primer and its alternative version are valid
                if ref_p.is_valid() and alt_p.is_valid():
                    # Store primers
                    pp = PrimerPair()
                    pp.chrom = marker.chrom
                    pp.snp_id = marker.name

                    pp.primers[direction] = ref_p
                    pp.primers["A"] = alt_p
                    pp.dir = direction

                    pps.append(pp)
            padding += 1

    return pps


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
        "marker",
        "chr",
        "pos",
        "ref",
        "alt",
        "start",
        "end",
        "direction",
        "type",
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
        "will_dimerize",
        "n_offtargets",
        "max_aaf",
        "indels",
        "offtargets",
        "mutations",
    ])

    with open(fp, "w") as f:
        f.write("{}\n".format(header))


def print_report(pcr, fp, delimiter="\t"):
    primer_type_ordering = ["REF", "ALT", "COM"]
    sorted_scores = np.sort(np.unique(list(pcr.pps_pruned.keys())))[::-1]
    ppid = 0

    seq_ids = {}
    for pt in primer_type_ordering:
        seq_ids[pt] = []

    with open(fp + ".txt", "w") as f, open(fp + ".bed", "w") as f_bed:
        for i in sorted_scores:
            for pp in pcr.pps_pruned[i]:
                pp.id = ppid

                dyes = {}
                seqs = {}
                seqs_amb = {}

                for d in pp.primers.keys():
                    dyes[d] = pp.primers[d].reporter

                if pp.dir == "F":
                    seqs["REF"] = pp.primers["F"].sequence[:-1].lower() + pp.primers["F"].sequence[-1]
                    seqs_amb["REF"] = pp.primers["F"].sequence_ambiguous[:-1].lower() + \
                                      pp.primers["F"].sequence_ambiguous[-1]
                    seqs["COM"] = pp.primers["R"].sequence.lower()
                    seqs_amb["COM"] = pp.primers["R"].sequence_ambiguous[:-1].lower() + \
                                      pp.primers["R"].sequence_ambiguous[-1]
                else:
                    seqs["COM"] = pp.primers["F"].sequence.lower()
                    seqs_amb["COM"] = pp.primers["F"].sequence_ambiguous.lower()
                    seqs["REF"] = pp.primers["R"].sequence[:-1].lower() + pp.primers["R"].sequence[-1]
                    seqs_amb["REF"] = pp.primers["R"].sequence_ambiguous[:-1].lower() + \
                                      pp.primers["R"].sequence_ambiguous[-1]

                seqs["ALT"] = pp.primers["A"].sequence[:-1].lower() + pp.primers["A"].sequence[-1]
                seqs_amb["ALT"] = pp.primers["A"].sequence_ambiguous[:-1].lower() + \
                                  pp.primers["A"].sequence_ambiguous[-1]

                curr_seq_ids = {}

                for d in primer_type_ordering:
                    if seqs[d] not in seq_ids[d]:
                        curr_seq_ids[d] = d + "_" + str(len(seq_ids[d]))
                        seq_ids[d].append(seqs[d])
                    else:
                        curr_seq_ids[d] = d + "_" + str(seq_ids[d].index(seqs[d]))

                    if d == pp.dir:
                        if seqs["ALT"] not in seq_ids[d]:
                            curr_seq_ids["ALT"] = d + "_" + str(len(seq_ids[d]))
                            seq_ids[d].append(seqs[d])
                        else:
                            curr_seq_ids["ALT"] = d + "_" + str(seq_ids[d].index(seqs[d]))

                # Shared properties
                offtargets = ",".join([x for x in pp.offtargets])
                if offtargets == "":
                    offtargets = "NA"

                # Find the order of the reported primers so that they are reported as REF/ALT/COM
                if pp.dir == "F":
                    type_ordering = ["F", "A", "R"]
                else:
                    type_ordering = ["R", "A", "F"]

                # Write primers parameters
                for ptype, d in zip(primer_type_ordering, type_ordering):
                    mutations = ",".join([x for x in pp.primers[d].mutations])

                    if d in ["A"]:
                        will_dimerize = pp.alt_dimer
                        direction = pp.dir
                    else:
                        will_dimerize = pp.ref_dimer
                        direction = d

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
                        int(pp.primers[d].start),
                        int(pp.primers[d].stop),
                        direction,
                        ptype,
                        pp.id,
                        dyes[d] + seqs[ptype],
                        dyes[d] + seqs_amb[ptype],
                        pcr.snp_id + "_" + curr_seq_ids[ptype],
                        pp.goodness,
                        pp.qcode,
                        pp.primers[d].length,
                        pp.product_size,
                        lib_utils.round_tidy(pp.primers[d].tm, 3),
                        lib_utils.round_tidy(pp.primers[d].gc_content, 3),
                        will_dimerize,
                        len(pp.offtargets),
                        lib_utils.round_tidy(pp.primers[d].max_aaf, 2),
                        pp.max_indel_size,
                        offtargets,
                        mutations,
                    ]
                    f.write("{}\n".format(delimiter.join([str(x) for x in fields])))

                    # BED file
                    if direction == "F":
                        direction = "+"
                    else:
                        direction = "-"

                    fields = [
                        pcr.chrom,
                        int(pp.primers[d].start),
                        int(pp.primers[d].stop),
                        "{}_{}-{}".format(pcr.snp_id, curr_seq_ids[ptype], pp.goodness),
                        "0",
                        direction,
                    ]
                    f_bed.write("{}\n".format("\t".join([str(x) for x in fields])))

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


def merge_primers(mpp, p3_primers):
    pps = []

    for pid, pp3 in enumerate(p3_primers):
        pp = deepcopy(mpp)
        pp.add_primer3_attributes(pp3)

        # Add start and stop to the alternative primer
        pp.primers["A"].start = pp.primers[pp.dir].start
        pp.primers["A"].stop = pp.primers[pp.dir].stop
        pp.primers["A"].id = pp.primers[pp.dir].id

        pps.append(pp)

    return pps


def design_primers(pps_repo, mpps, target_seq, target_start, ivs, n_primers=10):
    p3_repo = {}  # Primer pairs returned by PRIMER3

    for pid, mpp in enumerate(mpps):
        primer3_seq_args = {
            'SEQUENCE_TEMPLATE': target_seq,
            'SEQUENCE_EXCLUDED_REGION': ivs,
        }

        if mpp.dir == "F":
            primer3_seq_args['SEQUENCE_PRIMER'] = mpp.primers[mpp.dir].sequence
        else:
            primer3_seq_args['SEQUENCE_PRIMER_REVCOMP'] = mpp.primers[mpp.dir].sequence

        p3_primers = lib_primer3.get_primers(primer3_seq_args, target_start=target_start)
        p3_repo[pid] = merge_primers(mpp=mpp, p3_primers=p3_primers)  # Merge the already designed marker primer

    # Find valid primer pairs by checking top hits for each marker primers iteratively
    flag_continue = True
    while flag_continue:

        # Remove marker primers entries where no additional primer pairs are possible
        mpp_ids = list(p3_repo.keys())
        for mpp_id in mpp_ids:
            if len(p3_repo[mpp_id]) == 0:
                del p3_repo[mpp_id]

        if len(p3_repo) == 0:
            flag_continue = False

        for mpp_id, p3_pps in p3_repo.items():
            pp = p3_pps.pop(0)

            # Check if already present in the main repository
            is_new = True
            for rpp in pps_repo:
                cnt = 0
                for d in ["F", "R", "A"]:
                    if rpp.primers[d].sequence == pp.primers[d].sequence:
                        cnt += 1

                if cnt == 3:
                    is_new = False
                    break

            if is_new:
                # Check the complementary primer validity
                for d in ["F", "R"]:
                    if d != pp.dir:
                        is_valid = pp.primers[d].is_valid()
                        if is_valid:
                            pps_repo.append(pp)

                if len(pps_repo) == n_primers:  # We have enough primers
                    flag_continue = False
                    break

    return pps_repo


def main(kwarg_dict):
    # kwargs to variables
    fp_fasta = kwarg_dict["fp_fasta"]
    marker = kwarg_dict["marker"]
    fp_base_out = kwarg_dict["fp_base_out"]
    blast_db = kwarg_dict["blast_db"]
    muscle = kwarg_dict["muscle"]
    n_primers = kwarg_dict["n_primers"]
    p3_search_depth = kwarg_dict["p3_search_depth"]
    primer_seed = kwarg_dict["primer_seed"]
    tm_delta = kwarg_dict["tm_delta"]
    offtarget_size = kwarg_dict["offtarget_size"]
    primer3_configs = kwarg_dict["primer3_configs"]
    reporters = kwarg_dict["reporters"]
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

    # Find all marker primers
    mpps = get_marker_primers(marker=marker)
    logger_msg += "{:d} possible KASP marker primers found\n".format(len(mpps))

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
    maps["partial_mism"], maps["mism"] = map_homologs(fp_aligned, marker.fasta_name, marker.n)

    # Cleanup temporary files
    if not debug:
        os.remove(fp_fasta)
        os.remove(fp_aligned)

    # Get valid regions for the design of complementary primers
    ivs = {}
    for k, mmap in maps.items():
        k_mut = k + "_mut"

        ivs[k] = lib_primer3.get_exclusion_zone(mmap)  # Mutations not excluded

        if len(marker.mutations) > 0:
            # Make a list of mutation positions based on mutation for exclusion
            mut_ixs = []

            reg_start = marker.start
            for mutation in marker.mutations:
                mut_ixs.append(mutation.pos - reg_start)

            ivs[k_mut] = lib_primer3.get_exclusion_zone(mmap, hard_exclude=mut_ixs)

    # Loop across marker primers and design valid complementary primers using PRIMER3
    # PCR assay including multiple primer pairs
    pcr = PCR(marker.name, marker.chrom, marker.pos1, marker.ref, marker.alt)
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
            pcr.pps = design_primers(
                pps_repo=pcr.pps,  # Primer pair repository
                mpps=mpps,  # Marker primer pairs
                target_seq=marker.seq,
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

    pcr.add_reporter_dyes(reporters)  # Add the reporter dye sequence to the primers
    pcr.check_heterodimerization()  # Check for heterodimerization
    pcr.add_mutations(marker.mutations)  # List mutations in primers
    pcr.classify()  # Classify primers by scores using a heuristic "goodness" score
    n = pcr.prune(n_primers)  # Retain only the top n primers

    # Print to logger
    logger_msg += "Returned top {:d} primer pairs\n".format(n)
    logger.debug(logger_msg)
    print_report(pcr, join(fp_base_out))


if __name__ == "__main__":
    pass
