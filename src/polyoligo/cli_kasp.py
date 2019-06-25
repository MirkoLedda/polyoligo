from __future__ import print_function, division
import sys
import argparse
import logging
import os
from os.path import join, exists
import shutil
import multiprocessing as mp
import tqdm
from copy import deepcopy
import yaml
import cProfile

from . import lib_blast, _lib_kasp, _lib_markers, _logger_config, lib_utils, lib_vcf, _version

__version__ = _version.__version__

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "bin/win_x64"),
}

MARKER_FLANKING_N = 50  # Number of nucleotides on each sides when retrieving the sequence flanking the marker
MIN_ALIGN_LEN = 50  # Minimum alignment to declare homologs
MIN_ALIGN_ID = 88  # Minimum alignment identity to declare homologs
HOMOLOG_FLANKING_N = 250  # Number of nucleotides on each sides of the SNP when retrieving homolog sequences
WEBAPP_MAX_N = 100  # Limit for the number of markers to be designed when using the webapp mode


def cprofile_worker(kwargs):
    """Wrapper to cProfile subprocesses."""
    cProfile.runctx('_lib_kasp.main(kwargs)', globals(), locals(), 'profile_{}.out'.format(kwargs["marker"].name))


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="polyoligo-kasp",
                                     description="Design primers for Kompetitive allele specific PCR (KASP) assays",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "markers",
        metavar="MARKERS",
        type=str,
        help="File containing a list of input markers in tab-separated format and with no header as: "
             "CHR POS NAME REF ALT",
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
        "--fast",
        action="store_true",
        help="Run in fast mode. This mode loads the entire reference genome in memory making the design "
             "of large number of probes (> 1000) faster at the expense of heavy RAM consumption.",
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Silent mode, will not print to STDOUT.",
    )
    parser.add_argument(
        "-n", "--n_primers",
        metavar="<INT>",
        type=int,
        default=10,
        help="Maximum number of KASP primers to return for each marker.",
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
        "--report_alts",
        action="store_true",
        help="Report alternative subjects for each mutations in the VCF file? (a VCF file is needed for that option).",
    )
    parser.add_argument(
        "--depth",
        metavar="<FLOAT>",
        type=float,
        default=10.0,
        help="This parameter controls the exhaustiveness of the primer pair search, which is given by "
             "'n' * 'depth'. By increasing this value, more primer pairs will be considered but the process will "
             "be computationally heavier.",
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
        metavar="<TXT>",
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
        default=1000,
        help="Maximum size of offtarget PCR products.",
    )
    parser.add_argument(
        "-nt", "--n-tasks",
        metavar="<INT>",
        type=int,
        default=-1,
        help="Number of parallel threads. If negative, then all but '-nt' CPUs will be used with a minimum of 1 CPU.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        # help="Debug mode. FOR DEVS ONLY.",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--webapp",
        action="store_true",
        # help="Webapp mode. FOR DEVS ONLY.",
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args(inputargs)  # Parse input args

    return args


def main(strcmd=None):
    main_time = lib_utils.timer_start()  # Set main timer

    # Input arguments handling
    if strcmd:  # Means we are running the script using a string of arguments (e.g. for testing)
        testcmd = lib_utils.absolute_paths(strcmd)  # Make paths absolute
        args = parse_args(testcmd.split()[1:])
    else:
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
    _logger_config.setup_logging(log_fp=join(out_path, args.output + ".log"), verbose=not args.silent)
    logger = logging.getLogger(__name__)

    # Detect the os and point to respective binaries
    curr_os = lib_utils.get_os()
    if curr_os is None:
        logger.error("OS not supported or not detected properly.")
        sys.exit(1)
    bin_path = BINARIES[curr_os]

    # Init the BlastDB
    blast_db = lib_blast.BlastDB(
        path_db=args.refgenome,
        path_temporary=temp_path,
        path_bin=bin_path,
        job_id="main",
        n_cpu=args.n_tasks,
    )

    # Build the BlastDB if a fasta was provided
    if not blast_db.has_db:
        logger.info("Converting the input reference genome to BlastDB ...")
        blast_db.fasta2db()

    # Make a FASTA reference genome if fast mode is activated
    if args.fast:
        if not blast_db.has_fasta:
            logger.info("Converting the input reference genome to FASTA ...")
            blast_db.db2fasta()

        logger.info("Fast mode On: Loading the input reference genome in memory ...")
        blast_db.load_fasta()

    else:
        blast_db.has_fasta = False
        blast_db.purge()

    # Init the MUSCLE aligner
    muscle = lib_blast.Muscle(path_temporary=blast_db.temporary, exe=bin_path)

    # Read primer3 configs
    primer3_configs = {}
    if args.primer3 != "":
        with open(args.primer3, "r") as f:
            primer3_configs = yaml.safe_load(f)

    # Read reporter dyes
    reporters = [args.dye1, args.dye2]

    for reporter in reporters:
        if not lib_utils.is_dna(reporter):
            logger.error("Reporter dyes DNA sequence incorrect: {}".format(reporter))
            sys.exit(1)

    # Init Markers object
    markers = _lib_markers.Markers(
        blast_db=blast_db,
        MARKER_FLANKING_N=MARKER_FLANKING_N,
        MIN_ALIGN_LEN=MIN_ALIGN_LEN,
        MIN_ALIGN_ID=MIN_ALIGN_ID,
        HOMOLOG_FLANKING_N=HOMOLOG_FLANKING_N,
    )

    # Read markers
    try:
        markers.read_markers(args.markers)
    except:
        logger.error("Failed to read input markers. Please double check the format is CHR POS NAME REF ALT")
        sys.exit(1)

    if args.webapp:
        if len(markers) > WEBAPP_MAX_N:
            markers.markers = markers.markers[0:WEBAPP_MAX_N]
            logger.warning("Number of target markers exceeds the limit of {}. "
                           "Only designing first {} markers".format(WEBAPP_MAX_N, WEBAPP_MAX_N))

    logger.info("Number of target markers = {}".format(len(markers)))

    # Init a VCF object if a file is provided
    if args.vcf == "":
        vcf_obj = None
    else:
        logger.info("Loading VCF information ...")
        vcf_obj = lib_vcf.VCF(fp=args.vcf, fp_inc_samples=args.vcf_include, fp_exc_samples=args.vcf_exclude)

    # Get flanking sequences around the markers
    logger.info("Retrieving flanking sequence around markers ...")
    seqs = markers.get_marker_flanks()

    # Find homologs for each marker and write them to a unique fasta file per marker
    logger.info("Finding homeologs/duplications by sequence homology ...")
    markers.find_homologs(seqs)

    if vcf_obj is not None:
        logger.info("Uploading the VCF file ...")
        markers.upload_mutations(vcf_obj)  # Upload mutations for each markers from a VCF file (if provided)

        # Write subjects containing alternative alleles for each mutations
        if args.report_alts:
            args.report_alts = join(out_path, args.output + "_altlist.txt")
            logger.info("Writing subjects with alternative alleles -> {}".format(args.report_alts))

            if os.path.exists(args.report_alts):
                os.remove(args.report_alts)

            markers.print_alt_subjects(
                vcf_obj=vcf_obj,
                fp=args.report_alts,
            )

    # Start the search for candidate KASP primers
    if not args.webapp:
        logger.info("Searching for KASP candidates using {} parallel processes ...".format(args.n_tasks))
    else:
        logger.info("Searching for KASP candidates ...")

    # Purge sequences from BlastDB
    blast_db.purge()

    # Initialize a list of kwargs for the worker and a job queue
    kwargs_worker = []

    for marker in markers.markers:
        blast_db.job_id = marker.name  # To make sure files are not overwritten during multiprocessing
        blast_db.n_cpu = 1  # Blast should now only use 1 thread as we will be spawning jobs

        kwarg_dict = {
            "fp_fasta": join(temp_path, marker.name + ".fa"),
            "marker": marker,
            "fp_out": join(temp_path, marker.name + ".txt"),
            "blast_db": blast_db,
            "muscle": muscle,
            "reporters": reporters,
            "n_primers": args.n_primers,
            "p3_search_depth": args.depth,
            "tm_delta": args.tm_delta,
            "offtarget_size": [args.offtarget_min_size, args.offtarget_max_size],
            "primer_seed": args.seed,
            "primer3_configs": primer3_configs,
            "debug": args.debug,
        }
        kwargs_worker.append(deepcopy(kwarg_dict))

    # Run the job queue
    if not args.webapp:  # Normal mode - i.e. multithreaded
        n_jobs = len(kwargs_worker)
        with mp.Pool(processes=args.n_tasks, maxtasksperchild=1) as p:
            if args.debug:
                _ = list(tqdm.tqdm(
                    p.imap_unordered(cprofile_worker, kwargs_worker),
                    total=n_jobs,
                ))
            else:
                _ = list(tqdm.tqdm(
                    p.imap_unordered(_lib_kasp.main, kwargs_worker),
                    total=n_jobs,
                    disable=args.silent,
                ))
    else:  # Webapp mode - single jobs
        logger.info("nanobar - {:d}/{:d}".format(0, len(kwargs_worker)))
        for i, kwargs_job in enumerate(kwargs_worker):
            _lib_kasp.main(kwargs_job)
            logger.info("nanobar - {:d}/{:d}".format(i, len(kwargs_worker)))

    # Concatenate all primers for all markers into a single report
    logger.info("Preparing report ...")
    fp_out = join(out_path, args.output + ".txt")
    markers.write_report(fp_out)

    if not args.debug:
        shutil.rmtree(temp_path)

    if not args.webapp:
        logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
        logger.info("Report written to -> {}".format(fp_out))
    else:
        logger.info("Report ready !")


if __name__ == "__main__":
    main()
