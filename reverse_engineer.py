import os
from os.path import join, exists
import argparse
import sys
import shutil
import multiprocessing as mp
import logging
import pandas as pd
import numpy as np
import yaml
from copy import deepcopy
import tqdm

# noinspection PyProtectedMember
from src.polyoligo import lib_blast, lib_primer3, lib_markers, lib_vcf, lib_utils, _logger_config, _lib_pcr

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}

PRIMER3_DEFAULTS = join(os.path.dirname(__file__), "src/polyoligo/data/PRIMER3_KASP.yaml")


def read_assays(fp):
    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    df.columns = ["name", "chr", "start", "end", "dir", "id", "seq"]

    assays = {}

    for _, row in df.iterrows():

        if pd.isna(row["start"]):
            continue

        n = "{}-{}".format(row["name"], int(row["id"]))

        if n not in assays.keys():
            assays[n] = {}

        assays[n]["chrom"] = row["chr"]
        assays[n]["name"] = n
        assays[n][row["dir"]] = {
            "start": int(row["start"]),
            "end": int(row["end"]),
            "dir": row["dir"],
            "seq": row["seq"].upper(),
        }

    # Make sure we have primer pairs
    for n, assay in assays.items():
        if ("F" in assay.keys()) and ("R" in assay.keys()):
            continue
        else:
            del assays[n]

    return assays


def process_assay(kwargs):
    assay = kwargs["assay"]
    blast_hook = kwargs["blast_hook"]
    vcf_hook = kwargs["vcf_hook"]
    path_out = kwargs["path_out"]

    region = [
        np.min([assay["F"]["start"], assay["R"]["start"], assay["F"]["end"], assay["R"]["end"]]),
        np.max([assay["F"]["start"], assay["R"]["start"], assay["F"]["end"], assay["R"]["end"]]),
    ]

    roi = lib_markers.ROI(
        chrom=assay["chrom"],
        start=region[0],
        stop=region[1],
        blast_hook=blast_hook,
        vcf_hook=vcf_hook,
        name=assay["name"],
    )

    roi.upload_mutations()

    # Build primer pairs
    pp = lib_primer3.PrimerPair()

    for d in ["F", "R"]:
        pp.primers[d].start = assay[d]["start"]
        pp.primers[d].stop = assay[d]["end"]
        pp.primers[d].sequence = assay[d]["seq"]
    pp.chrom = roi.chrom
    pp.product_size = region[1] - region[0] + 1

    # Build PCR assay
    pcr = lib_primer3.PCR(
        name=assay["name"],
        chrom=roi.chrom,
        primer_pairs=[pp],
        snp_id=assay["name"],
    )

    # Check attributes
    _, _ = pcr.check_offtargeting(blast_hook, debug=False)  # Check for off-targets
    pcr.add_mutations(roi.mutations)
    pcr.check_heterodimerization()  # Check for heterodimerization
    pcr.classify(assay_name="PCR")  # Classify primers by scores using a heuristic "goodness" score
    pcr.prune(np.inf)

    # Write sub-report
    _lib_pcr.print_report(pcr, join(path_out, pcr.snp_id))


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="polyoligo-reveng",
                                     description="Reverse-engineer primers designed externally.",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format("0.1"),
    )
    parser.add_argument(
        "markers",
        metavar="MARKERS",
        type=str,
        help="File containing a list of input primers in space-separated format and with no header as: "
             "NAME CHR START END DIR ID SEQ",
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        type=str,
        help="Output filename.",
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
    _lib_pcr.print_report(pcr, kwargs_dict["fp_out"])


def write_final_reports(fp_base_out, assays, blast_hook):
    _lib_pcr.print_report_header(fp_base_out + ".txt")
    with open(fp_base_out + ".txt", "a") as f:
        for assay in assays.values():
            with open(join(blast_hook.temporary, assay["name"] + ".txt"), "r") as f_marker:
                for line in f_marker:
                    f.write(line)


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
        n_cpu=args.n_tasks,
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

    assays = read_assays(args.markers)

    # Stop VCF hook reader for serialization
    if vcf_hook:
        vcf_hook.stop_reader()

    # Parallelized loop
    # Initialize a list of kwargs for the worker and a job queue
    kwargs_worker = []

    for assay in assays.values():
        blast_hook.job_id = assay["name"]  # To make sure files are not overwritten during multiprocessing

        kwarg_dict = {
            "assay": assay,
            "blast_hook": blast_hook,
            "vcf_hook": vcf_hook,
            "path_out": temp_path,
        }

        kwargs_worker.append(deepcopy(kwarg_dict))

    # Run the job queue
    n_jobs = len(kwargs_worker)
    with mp.Pool(processes=args.n_tasks, maxtasksperchild=1) as p:
        if args.debug:
            _ = list(tqdm.tqdm(
                p.imap_unordered(process_assay, kwargs_worker),
                total=n_jobs,
            ))
        else:
            _ = list(tqdm.tqdm(
                p.imap_unordered(process_assay, kwargs_worker),
                total=n_jobs,
                disable=args.silent,
            ))

    write_final_reports(join(out_path, args.output), assays, blast_hook)

    if not args.debug:
        shutil.rmtree(temp_path)

    logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
    logger.info("Report written to -> {}".format(join(out_path, args.output) + ".txt"))


if __name__ == "__main__":
    main()
