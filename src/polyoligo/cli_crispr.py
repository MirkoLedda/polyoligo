from __future__ import print_function, division
import sys
import argparse
import logging
import os
from os.path import join, exists
import shutil
import multiprocessing as mp
import cProfile

from . import lib_blast, _lib_crispr, _logger_config, lib_utils, _version

__version__ = _version.__version__

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "bin/win_x64"),
}


def cprofile_worker(kwargs):
    """Wrapper to cProfile subprocesses."""
    cProfile.runctx('_getcrispr.main(kwargs)', globals(), locals(), 'profile_{}.out'.format(kwargs["marker"].name))


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="polyoligo-crispr",
                                     description="Design guide RNAs for CRISPR/Cas9 assays.",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "roi",
        metavar="ROI",
        type=str,
        help="Target region as CHR:START-END.",
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
        "--pam",
        metavar="<TXT>",
        type=str,
        default="NGG",
        # help="PAM site.",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Silent mode, will not print to STDOUT.",
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

    # # Make a FASTA reference genome by default
    # if not blast_db.has_fasta:
    #     logger.info("Converting the input reference genome to FASTA ...")
    #     blast_db.db2fasta()
    #
    # logger.info("Loading the input reference genome in memory ...")
    # blast_db.load_fasta()

    # Make a FASTA reference genome if fast mode is activated
    blast_db.has_fasta = False
    blast_db.purge()

    # Check that the PAM site is valid
    if not lib_utils.is_dna(args.pam):
        logger.error("PAM site invalid: {}".format(args.pam))
        sys.exit(1)

    # Init Crispr object and parse ROI
    crispr = _lib_crispr.Crispr(
        roi=args.roi,
        pam=args.pam,
        blast_db=blast_db)
    crispr.fetch_roi()

    if args.webapp:
        logger.info("nanobar - {:d}/{:d}".format(0, 1))

    # Find putative guide RNAs
    logger.info("Searching guide RNAs - PAM: {} - Region size: {} nts ...".format(crispr.pam, len(crispr.seq)))
    crispr.find_gRNAs()
    logger.info("Found {} possible guide RNAs".format(len(crispr.gRNAs)))

    if args.webapp:
        logger.info("nanobar - {:d}/{:d}".format(int(len(crispr.gRNAs)/3), len(crispr.gRNAs)))

    # Check for offtargeting
    logger.info("Determining offtargets - this may take a while ...")
    crispr.check_offtargeting(logger=logger, webapp=args.webapp)

    # Check seed regions for homologs
    logger.info("Finding seed homologs - this may also take a while ...")
    crispr.check_seeds(logger=logger,  webapp=args.webapp)

    # Compute some gRNA statistics
    logger.info("Computing gRNA statistics ...")
    crispr.calc_TM()  # Get melting temperatures of gRNAs
    crispr.get_run_Ts()  # Identify runs of TTTT

    # Print report
    logger.info("Preparing report ...")
    fp_report = join(out_path, args.output + ".txt")
    fp_report_bed = join(out_path, args.output + ".bed")
    crispr.write_report_header(fp_report, fp_report_bed)
    crispr.write_report(fp_report, fp_report_bed)

    if not args.debug:
        shutil.rmtree(temp_path)

    if not args.webapp:
        logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
        logger.info("Report written to -> {}".format(fp_report))
    else:
        logger.info("Report ready !")


if __name__ == "__main__":
    main()
