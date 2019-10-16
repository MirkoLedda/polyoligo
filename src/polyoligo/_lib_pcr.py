from __future__ import print_function, division
# noinspection PyPackageRequirements
import numpy as np
# noinspection PyUnresolvedReferences
from primer3.thermoanalysis import ThermoAnalysis
import logging
from os.path import join

from . import lib_blast, lib_utils, lib_primer3, lib_markers

# GLOBALS
INNER_LEN = 100  # Length of the region used to find homologs
OUTER_LEN = 500  # Length of the region were homologs will be mapped
MIN_ALIGN_ID = 88  # Minimum alignment identity to declare homologs
MIN_ALIGN_LEN = 50  # Minimum alignment to declare homologs
DELIMITER = "\t"  # Delimiter for output files
HEADER = [
    "target",
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
]

# Initialize the logger
logger = logging.getLogger(__name__)


def write_final_reports(fp_base_out, rois):
    print_report_header(fp_base_out + ".txt")
    with open(fp_base_out + ".txt", "a") as f:
        for roi in rois:
            with open(join(roi.blast_hook.temporary, roi.name + ".txt"), "r") as f_marker:
                for line in f_marker:
                    f.write(line)

    # BED file
    line_memory = []
    with open(fp_base_out + ".bed", "w") as f:
        for roi in rois:
            with open(join(roi.blast_hook.temporary, roi.name + ".bed"), "r") as f_marker:
                for line in f_marker:
                    if line not in line_memory:
                        line_memory.append(line)
                        f.write(line)


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


def print_report_header(fp):
    header = DELIMITER.join(HEADER)

    with open(fp, "w") as f:
        f.write("{}\n".format(header))


def print_report(pcr, fp):
    primer_type_ordering = ["F", "R"]
    sorted_scores = np.sort(np.unique(list(pcr.pps_pruned.keys())))[::-1]
    ppid = 0

    seq_ids = {}
    for pt in primer_type_ordering:
        seq_ids[pt] = []

    with open(fp + ".txt", "w") as f, open(fp + ".bed", "w") as f_bed:
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
                        pcr.name,
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
                    f.write("{}\n".format(DELIMITER.join([str(x) for x in fields])))

                    # BED file
                    if d == "F":
                        direction = "+"
                    else:
                        direction = "-"

                    fields = [
                        pcr.chrom,
                        int(pp.primers[d].start),
                        int(pp.primers[d].stop),
                        "{}_{}-{}".format(pcr.name, curr_seq_ids[d], pp.goodness),
                        "0",
                        direction,
                    ]
                    f_bed.write("{}\n".format("\t".join([str(x) for x in fields])))

                f.write("\n")
                ppid += 1

        if ppid == 0:
            fields = [
                pcr.name,
                pcr.chrom,
            ]
            fields = [str(x) for x in fields]
            f.write(DELIMITER.join(fields + (len(HEADER)-2) * ["NA"]) + "\n")
            f.write("\n")

        f.write("\n")


def design_primers(pps_repo, roi, sequence_target=None,
                   sequence_excluded_region=None, n_primers=10):
    primer3_seq_args = {'SEQUENCE_TEMPLATE': roi.seq}

    if sequence_target is not None:
        primer3_seq_args['SEQUENCE_TARGET'] = sequence_target

    if sequence_excluded_region is not None:
        primer3_seq_args['SEQUENCE_EXCLUDED_REGION'] = sequence_excluded_region

    p3_repo = lib_primer3.get_primers(primer3_seq_args, target_start=roi.start)

    # Find valid primer pairs by checking top hits for each marker primers iteratively
    flag_continue = True
    while flag_continue:

        if len(p3_repo) == 0:
            flag_continue = False

        for pp in p3_repo:

            if len(pps_repo) == n_primers:  # We have enough primers
                flag_continue = False
                break

            pp.chrom = roi.chrom
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
    region = kwarg_dict["region"]
    n_primers = kwarg_dict["n_primers"]
    p3_search_depth = kwarg_dict["p3_search_depth"]
    primer_seed = kwarg_dict["primer_seed"]
    tm_delta = kwarg_dict["tm_delta"]
    offtarget_size = kwarg_dict["offtarget_size"]
    primer3_configs = kwarg_dict["primer3_configs"]
    debug = kwarg_dict["debug"]
    fp_base_out = kwarg_dict["fp_base_out"]

    if region.vcf_hook:
        region.vcf_hook.start_reader()

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

    # Define two rois around the target where the primers will be built
    region.left_roi = lib_markers.ROI(
        chrom=region.chrom,
        start=region.start,
        stop=region.start,
        blast_hook=region.blast_hook,
        malign_hook=region.malign_hook,
        vcf_hook=region.vcf_hook,
        name=region.name,
        marker=region.marker,
        do_print_alt_subjects=region.do_print_alt_subjects,
    )
    region.right_roi = lib_markers.ROI(
        chrom=region.chrom,
        start=region.stop,
        stop=region.stop,
        blast_hook=region.blast_hook,
        malign_hook=region.malign_hook,
        vcf_hook=region.vcf_hook,
        name=region.name,
        marker=region.marker,
        do_print_alt_subjects=region.do_print_alt_subjects,
    )

    # Find homologs
    region.left_roi = region.left_roi.map_homologs(
        inner_len=INNER_LEN,
        outer_len=OUTER_LEN,
        min_align_id=MIN_ALIGN_ID,
        min_align_len=MIN_ALIGN_LEN,
    )
    region.right_roi = region.right_roi.map_homologs(
        inner_len=INNER_LEN,
        outer_len=OUTER_LEN,
        min_align_id=MIN_ALIGN_ID,
        min_align_len=MIN_ALIGN_LEN,
    )

    # Find mutations in the region
    region.left_roi.upload_mutations()
    lib_primer3.include_mut_in_included_maps(region.left_roi)
    region.right_roi.upload_mutations()
    lib_primer3.include_mut_in_included_maps(region.right_roi)

    # Combine the maps with the target region
    region.merge_with_primers()

    # Build exclusion maps
    lib_primer3.get_sequence_excluded_regions(region)

    # Print alternative subjects if required TODO

    # Design primers
    pcr = lib_primer3.PCR(chrom=region.chrom, name=region.name)
    map_names = {
        "mism_mut": "Variants excluded | Homologs specific",
        "partial_mism_mut": "Variants excluded | Homologs partial spec",
        "all_mut": "Variants excluded | Unspecific",
        "mism": "Variants included | Homologs specific",
        "partial_mism": "Variants included | Homologs partial spec",
        "all": "Variants included | Unspecific",
    }

    # Set the search mask ordering
    if "mism_mut" in region.p3_sequence_included_maps.keys():
        search_types = ["mism_mut", "mism", "partial_mism_mut", "partial_mism", "all_mut", "all"]
    else:
        search_types = ["mism", "partial_mism", "all"]

    for search_type in search_types:
        if len(pcr.pps) < lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN']:
            n_before = len(pcr.pps)
            pcr.pps = design_primers(
                pps_repo=pcr.pps,  # Primer pair repository
                roi=region,
                sequence_excluded_region=region.p3_sequence_excluded_regions[search_type],
                n_primers=lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN'],
            )
            n_new = len(pcr.pps) - n_before
            logger_msg += "{:42}: {:3d} pairs\n".format(map_names[search_type], n_new)

    hit_cnts, valid_hit_cnts = pcr.check_offtargeting(region.blast_hook, debug=debug)  # Check for off-targets

    # Report hit counts
    logger_msg += "Offtarget hits across the genome\n"
    logger_msg += "      Forward : {:d}\n".format(valid_hit_cnts["F"])
    logger_msg += "      Reverse : {:d}\n".format(valid_hit_cnts["R"])

    pcr.check_heterodimerization()  # Check for heterodimerization
    pcr.add_mutations(region.mutations)  # List mutations in primers
    pcr.classify()  # Classify primers by scores using a heuristic "goodness" score
    n = pcr.prune(n_primers)  # Retain only the top n primers

    # Print to logger
    logger_msg += "Returned top {:d} primer pairs\n".format(n)
    logger.debug(logger_msg)

    print_report(pcr, fp_base_out)


if __name__ == "__main__":
    pass
