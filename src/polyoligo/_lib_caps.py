from __future__ import print_function, division
import numpy as np
# noinspection PyUnresolvedReferences
from primer3.thermoanalysis import ThermoAnalysis
import logging
from os.path import join
import os
import json

from . import lib_utils, lib_primer3, lib_markers, _lib_pcr

# GLOBALS
ENZYME_FILENAME = join(os.path.dirname(__file__), "data/type2.json")
INNER_LEN = 100  # Length of the region used to find homologs
OUTER_LEN = 1000  # Length of the region were homologs will be mapped
MIN_ALIGN_ID = 88  # Minimum alignment identity to declare homologs
MIN_ALIGN_LEN = 50  # Minimum alignment to declare homologs
DELIMITER = "\t"  # Delimiter for output files
HEADER = [
    "marker",
    "chr",
    "pos",
    "ref",
    "alt",
    "enzymes",
    "enzyme_suppliers",
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
    "allele_cut",
    "fragment_1",
    "fragment_2",
    "tm",
    "gc_content",
    "n_offtargets",
    "max_aaf",
    "indels",
    "offtarget_products",
    "mutations",
    "PCR_product",
]

# Initialize the logger
logger = logging.getLogger(__name__)


def get_enzymes(included_enzymes):
    # Load restriction enzymes
    with open(ENZYME_FILENAME, "r") as f:
        enzymes = json.load(f)

    # Filter enzymes
    valid_enzymes = []
    if included_enzymes is not None:
        for enzyme in enzymes:
            if enzyme["name"] in included_enzymes:
                valid_enzymes.append(enzyme)
    else:
        valid_enzymes = enzymes
    return valid_enzymes


class PCR(lib_primer3.PCR):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


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

                # Loop across restriction enzymes
                for enzyme in pp.enzymes:

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
                            enzyme["name"],
                            enzyme["supplier"],
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
                            enzyme["allele"],
                            enzyme["fragl"],
                            enzyme["fragr"],
                            lib_utils.round_tidy(pp.primers[d].tm, 3),
                            lib_utils.round_tidy(pp.primers[d].gc_content, 3),
                            len(pp.offtargets),
                            lib_utils.round_tidy(pp.primers[d].max_aaf, 2),
                            pp.max_indel_size,
                            offtargets,
                            mutations,
                            pp.seq_x,
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
                            "{}_{}-{}".format(pcr.snp_id, curr_seq_ids[d], pp.goodness),
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
            f.write(DELIMITER.join(fields + 24 * ["NA"]) + "\n")
            f.write("\n")

        f.write("\n")


def write_final_reports(fp_base_out, rois):
    print_report_header(fp_base_out + ".txt")
    with open(fp_base_out + ".txt", "a") as f:
        for roi in rois:
            with open(join(roi.blast_hook.temporary, roi.marker.name + ".txt"), "r") as f_marker:
                for line in f_marker:
                    f.write(line)

    # BED file
    line_memory = []
    with open(fp_base_out + ".bed", "w") as f:
        for roi in rois:
            with open(join(roi.blast_hook.temporary, roi.marker.name + ".bed"), "r") as f_marker:
                for line in f_marker:
                    if line not in line_memory:
                        line_memory.append(line)
                        f.write(line)


def main(kwarg_dict):
    # kwargs to variables
    region = kwarg_dict["region"]
    fp_base_out = kwarg_dict["fp_base_out"]
    n_primers = kwarg_dict["n_primers"]
    p3_search_depth = kwarg_dict["p3_search_depth"]
    primer_seed = kwarg_dict["primer_seed"]
    tm_delta = kwarg_dict["tm_delta"]
    offtarget_size = kwarg_dict["offtarget_size"]
    primer3_configs = kwarg_dict["primer3_configs"]
    included_enzymes = kwarg_dict["included_enzymes"]
    debug = kwarg_dict["debug"]
    fragment_min_size = kwarg_dict["fragment_min_size"]

    if region.vcf_hook:
        region.vcf_hook.start_reader()

    # Set primer3 globals
    if len(primer3_configs) > 0:
        lib_primer3.set_globals(**primer3_configs)

    lib_primer3.set_globals(
        PRIMER_NUM_RETURN=p3_search_depth,
    )
    lib_primer3.set_tm_delta(tm_delta)
    lib_primer3.set_primer_seed(primer_seed)

    # Set lib_primer offtarget sizes
    lib_primer3.set_offtarget_size(offtarget_size[0], offtarget_size[1])

    # Set a logger message that will be printed at the end (to be threadsafe)
    header = "Primer search results for {}".format(region.name)
    sepline = "=" * (len(header) + 2)
    logger_msg = "\n{}\n{}\n{}\n".format(sepline, header, sepline)

    # List all considered restriction enzymes
    enzymes = get_enzymes(included_enzymes=included_enzymes)
    logger_msg += "Number of considered restriction enzymes: {}\n".format(len(enzymes))

    # Define two rois around the target where the primers will be built
    region.left_roi = lib_markers.ROI(
        chrom=region.chrom,
        start=int(region.start-(OUTER_LEN/2)),
        stop=int(region.start-(OUTER_LEN/2)),
        blast_hook=region.blast_hook,
        malign_hook=region.malign_hook,
        vcf_hook=region.vcf_hook,
        name=region.name,
        marker=region.marker,
        do_print_alt_subjects=region.do_print_alt_subjects,
    )
    region.right_roi = lib_markers.ROI(
        chrom=region.chrom,
        start=int(region.start+(OUTER_LEN/2)),
        stop=int(region.start+(OUTER_LEN/2)),
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

    # Ensure we have room around the marker
    region = lib_markers.Region(
        chrom=region.chrom,
        start=region.marker.pos1-fragment_min_size,
        stop=region.marker.pos1+fragment_min_size,
        blast_hook=region.blast_hook,
        malign_hook=region.malign_hook,
        vcf_hook=region.vcf_hook,
        name=region.name,
        marker=region.marker,
        do_print_alt_subjects=region.do_print_alt_subjects,
        left_roi=region.left_roi,
        right_roi=region.right_roi,
    )

    # Combine the maps with the target region
    region = region.merge_with_primers()

    # Build exclusion maps
    lib_primer3.get_sequence_excluded_regions(region)

    # Design primers
    pcr = PCR(
        snp_id=region.marker.name,
        chrom=region.marker.chrom,
        pos=region.marker.pos1,
        ref=region.marker.ref,
        alt=region.marker.alt,
    )
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

    lib_primer3.set_globals(PRIMER_PRODUCT_SIZE_RANGE=[[region.p3_sequence_target[0][1], region.n]])

    for search_type in search_types:
        if len(pcr.pps) < lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN']:
            n_before = len(pcr.pps)
            pcr.pps = _lib_pcr.design_primers(
                pps_repo=pcr.pps,  # Primer pair repository
                roi=region,
                sequence_target=region.p3_sequence_target,
                sequence_excluded_region=region.p3_sequence_excluded_regions[search_type],
                n_primers=lib_primer3.PRIMER3_GLOBALS['PRIMER_NUM_RETURN'],
                max_unique=2,  # This ensure diversity in the output
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

    # Prune only valid pcr_products
    pruned_pps = []
    for pp in pcr.pps:
        pcr_product = lib_markers.PCRproduct(roi=region, pp=pp)
        pcr_product.will_it_cut(enzymes=enzymes, min_fragment_size=fragment_min_size)

        if len(pcr_product.enzymes) > 0:
            pp.enzymes = pcr_product.enzymes
            pp.seq_x = pcr_product.roi.seq_x
            pruned_pps.append(pp)
    pcr.pps = pruned_pps
    pcr.classify(assay_name="CAPS")  # Classify primers by scores using a heuristic "goodness" score
    n = pcr.prune(n_primers)  # Retain only the top n primers

    # Print to logger
    logger_msg += "Returned top {:d} primer pairs\n".format(n)
    logger.debug(logger_msg)

    print_report(
        pcr=pcr,
        fp=fp_base_out,
    )


if __name__ == "__main__":
    pass
