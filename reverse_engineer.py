import os
from os.path import join, exists
import argparse
import sys
import multiprocessing as mp
import logging
import pandas as pd
import numpy as np
import yaml

# noinspection PyProtectedMember
from src.polyoligo import lib_blast, _lib_kasp, lib_markers, lib_utils, lib_vcf, _logger_config, lib_primer3

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}

PRIMER3_DEFAULTS = join(os.path.dirname(__file__), "src/polyoligo/data/PRIMER3_KASP.yaml")


def process_markers(fp, blast_hook, vcf_hook, reporters, fp_out):
    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    df.columns = ["name", "chr", "pos", "ref", "alt", "start", "stop", "dir", "type", "id", "seq"]

    assays = {}
    pcrs = []
    rois = []

    for _, row in df.iterrows():

        if pd.isna(row["type"]):
            continue

        n = "{}-{}".format(row["name"], int(row["id"]))

        if n not in assays.keys():
            assays[n] = {}

        assays[n]["chr"] = row["chr"]
        assays[n]["pos"] = row["pos"]
        assays[n]["ref"] = row["ref"]
        assays[n]["alt"] = row["alt"]

        # Remove the current reporter dye
        seq = np.array(list(row["seq"]))

        seqp = []
        for i, s in enumerate(seq):
            if s.islower():
                seqp = seq[i:]
                break
        seq = "".join(list(seqp))

        assays[n][row["type"]] = {
            "start": int(row["start"]),
            "stop": int(row["stop"]),
            "dir": row["dir"],
            "seq": seq.upper(),
        }

    for n, assay in assays.items():
        if ("REF" in assay.keys()) and ("ALT" in assay.keys()) and ("COM" in assay.keys()):

            region = [
                np.min([assay["REF"]["start"], assay["ALT"]["start"], assay["COM"]["start"]]),
                np.max([assay["REF"]["stop"], assay["ALT"]["stop"], assay["COM"]["stop"]]),
            ]

            marker = lib_markers.Marker(
                name=n,
                chrom=assay["chr"],
                pos=assay["pos"]-1,
                ref_allele=assay["ref"],
                alt_allele=assay["alt"],
            )

            roi = lib_markers.ROI(
                name=n,
                chrom=assay["chr"],
                start=region[0],
                stop=region[1],
                marker=marker,
                blast_hook=blast_hook,
                vcf_hook=vcf_hook,
            )

            roi.upload_mutations()

            pp = _lib_kasp.PrimerPair()
            if assay["REF"]["dir"] == "F":
                mapping = {"F": "REF", "R": "COM", "A": "ALT"}
            else:
                mapping = {"F": "COM", "R": "REF", "A": "ALT"}

            for d in ["F", "R", "A"]:
                pp.primers[d].start = assay[mapping[d]]["start"]
                pp.primers[d].stop = assay[mapping[d]]["stop"]
                pp.primers[d].sequence = assay[mapping[d]]["seq"]
            pp.dir = assay["REF"]["dir"]
            pp.chrom = roi.chrom

            pcr = _lib_kasp.PCR(
                snp_id=roi.marker.name,
                chrom=roi.chrom,
                pos=roi.marker.pos1,
                ref=roi.marker.ref,
                alt=roi.marker.alt,
                primer_pairs=[pp],
            )

            _, _ = pcr.check_offtargeting(blast_hook, debug=False)  # Check for off-targets
            pcr.add_mutations(roi.mutations)
            pcr.add_reporter_dyes(reporters)  # Add the reporter dye sequence to the primers
            pcr.check_heterodimerization()  # Check for heterodimerization
            pcr.classify(assay_name="KASP")  # Classify primers by scores using a heuristic "goodness" score
            pcr.prune(np.inf)
            pcrs.append(pcr)
            rois.append(roi)

    return rois, pcrs


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="polyoligo",
                                     description="Design primers for Kompetitive allele specific PCR (KASP) assays",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s 0.0",
    )
    parser.add_argument(
        "markers",
        metavar="MARKERS",
        type=str,
        help="File containing a list of input primers in space-separated format and with no header as: "
             "NAME CHR POS REF ALT START END DIR TYPE ID SEQ",
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        type=str,
        help="Output filename (no extension).",
    )
    parser.add_argument(
        "refgenome",
        metavar="FASTA/BLASTDB",
        type=str,
        help="Either a FASTA file or a BLAST database to use as the reference genome.\n"
             "FASTA: Both raw and compressed files are supported natively (see sample_data/blastdb.fa.gz).\n"
             "BLASTDB: Extensions (e.g. .nsq/.nin/...) are not required, just enter the basename of the "
             "database (see sample_data/blastdb).",
    )
    parser.add_argument(
        "--dye1",
        metavar="<TXT>",
        type=str,
        default="GAAGGTCGGAGTCAACGGATT",
        help="DNA sequence of the reporter dye for the reference allele in a 5'-3' orientation. Default: VIC",
    )
    parser.add_argument(
        "--dye2",
        metavar="<TXT>",
        type=str,
        default="GAAGGTGACCAAGTTCATGCT",
        help="DNA sequence of the reporter dye for the alternative allele in a 5'-3' orientation. Default: FAM",
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Silent mode, will not print to STDOUT.",
    )
    parser.add_argument(
        "--vcf",
        metavar="<VCF>",
        type=str,
        default="",
        help="VCF file. Note that a tabix index file (.tbi) associated with the VCF is required.",
    )
    parser.add_argument(
        "--vcf_include",
        metavar="<TXT>",
        type=str,
        default="",
        help="List of samples (as listed in the VCF file) to consider during primer design.",
    )
    parser.add_argument(
        "--vcf_exclude",
        metavar="<TXT>",
        type=str,
        default="",
        help="List of samples (as listed in the VCF file) to exclude during primer design. "
             "As priority over '--vcf_include'.",
    )
    parser.add_argument(
        "--seed",
        metavar="<FLOAT>",
        type=int,
        default=12,
        help="Length of the primer 3'-end seed that will be considered when searching potential offtargets.",
    )
    parser.add_argument(
        "--primer3",
        metavar="<YAML>",
        type=str,
        default="",
        help="Configuration file for PRIMER3 in YAML format. Uses the same convention as PRIMER3. "
             "See sample_data/PRIMER3.yaml for a list of default values.",
    )
    parser.add_argument(
        "--tm_delta",
        metavar="<FLOAT>",
        type=float,
        default=5.0,
        help="The minimum difference in melting temperature between the primer tm and the various "
             "structures that could form (homo/heterodimer and hairpins).",
    )
    parser.add_argument(
        "-nt", "--n-tasks",
        metavar="<INT>",
        type=int,
        default=-1,
        help="Number of parallel threads. If negative, then all but '-nt' CPUs will be used with a minimum of 1 CPU.",
    )
    parser.add_argument(
        "--offtarget_min_size",
        metavar="<INT>",
        type=int,
        default=0,
        help="Minimum size of offtarget PCR products.",
    )
    parser.add_argument(
        "--offtarget_max_size",
        metavar="<INT>",
        type=int,
        default=10000,
        help="Maximum size of offtarget PCR products.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        # help="Debug mode. FOR DEVS ONLY.",
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args(inputargs)  # Parse input args

    return args


def score_marker(kwargs_dict):
    marker = kwargs_dict["marker"]
    pcr = marker.score_primers(kwargs_dict["fp_reporters"], kwargs_dict["blast_db"])
    _lib_kasp.print_report(pcr, kwargs_dict["fp_out"])


def main():
    main_time = lib_utils.timer_start()  # Set main timer
    args = parse_args(sys.argv[1:])

    # Set the number of CPUs
    if args.n_tasks < 1:
        args.n_tasks = mp.cpu_count() + args.n_tasks

    if args.n_tasks < 1:  # To make sure at least 1 CPU is used
        args.n_tasks = 1

    # Prepare directories
    out_path = os.path.dirname(os.path.abspath(args.output))
    if out_path == "":
        out_path = os.getcwd()
    args.output = os.path.basename(args.output)
    temp_path = join(out_path, "temporary")

    if not exists(out_path):
        os.makedirs(out_path)

    if not exists(temp_path):
        os.makedirs(temp_path)

    # Init the logger
    _logger_config.setup_console_logging(verbose=not args.silent)
    logger = logging.getLogger(__name__)

    # Detect the os and point to respective binaries
    curr_os = lib_utils.get_os()
    if curr_os is None:
        logger.error("OS not supported or not detected properly.")
        sys.exit()
    bin_path = BINARIES[curr_os]

    # Initialize hooks
    blast_hook = lib_blast.BlastDB(
        path_db=args.refgenome,
        path_temporary=temp_path,
        path_bin=bin_path,
        job_id="main",
        n_cpu=1,
    )

    if args.vcf == "":
        vcf_hook = None
    else:
        logger.info("Loading VCF information ...")
        vcf_hook = lib_vcf.VCF(fp=args.vcf, fp_inc_samples=args.vcf_include, fp_exc_samples=args.vcf_exclude)

    # Read primer3 configs
    primer3_configs = {}
    if args.primer3 != "":
        with open(args.primer3, "r") as f:
            primer3_configs = yaml.safe_load(f)

    # Set primer3 default values for unset values
    with open(PRIMER3_DEFAULTS, "r") as f:
        primer3_defaults = yaml.safe_load(f)
        for k, v in primer3_defaults.items():
            if k not in primer3_configs.keys():  # overwrite only if not set
                primer3_configs[k] = v

    # Set primer3 globals
    if len(primer3_configs) > 0:
        lib_primer3.set_primer_seed(args.seed)
    lib_primer3.set_tm_delta(args.tm_delta)
    lib_primer3.set_globals(**primer3_configs)
    lib_primer3.set_offtarget_size(args.offtarget_min_size, args.offtarget_max_size)

    # Read reporter dyes
    reporters = [args.dye1, args.dye2]

    for reporter in reporters:
        if not lib_utils.is_dna(reporter):
            logger.error("Reporter dyes DNA sequence incorrect: {}".format(reporter))
            sys.exit(1)

    rois, pcrs = process_markers(
        fp=args.markers,
        blast_hook=blast_hook,
        vcf_hook=vcf_hook,
        reporters=reporters,
        fp_out=join(out_path, args.output),
    )

    for pcr in pcrs:
        _lib_kasp.print_report(pcr, join(temp_path, pcr.snp_id))
        break

    _lib_kasp.write_final_reports(join(out_path, args.output), rois)

    logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
    logger.info("Report written to -> {} [{}, {}]".format(join(out_path, args.output) + ".txt", ".bed", ".log"))


if __name__ == "__main__":
    main()
