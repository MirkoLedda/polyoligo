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

from . import lib_blast, _lib_pcr, _logger_config, lib_utils, lib_vcf, _version, _lib_markers

__version__ = _version.__version__

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "bin/win_x64"),
}

PRIMER3_DEFAULTS = join(os.path.dirname(__file__), "data/PRIMER3_PCR.yaml")
MARKER_FLANKING_N = 50  # Number of nucleotides on each sides when retrieving the sequence flanking the marker
MIN_ALIGN_LEN = 50  # Minimum alignment to declare homologs
MIN_ALIGN_ID = 88  # Minimum alignment identity to declare homologs
# HOMOLOG_FLANKING_N = 250  # Number of nucleotides on each sides of the SNP when retrieving homolog sequences
WEBAPP_MAX_N = 100  # Limit for the number of markers to be designed when using the webapp mode


def cprofile_worker(kwargs):
    """Wrapper to cProfile subprocesses."""
    cProfile.runctx('_getkasp.main(kwargs)', globals(), locals(), 'profile_{}.out'.format(kwargs["marker"].name))


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="polyoligo-pcr",
                                     description="Design primers for PCR assays",
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
        help="Target region around which to design primers. Declared as CHR:START-END.",
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
        "--fast",
        action="store_true",
        # help="Run in fast mode. This mode loads the entire reference genome in memory making the design "
        #      "of large number of probes (> 1000) faster at the expense of heavy RAM consumption.",
        help=argparse.SUPPRESS,
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
        help="Maximum number of primer pairs to return.",
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
        # help="The minimum difference in melting temperature between the primer tm and the various "
        #      "structures that could form (homo/heterodimer and hairpins).",
        help=argparse.SUPPRESS,
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

    # Set primer3 default values for unset values
    with open(PRIMER3_DEFAULTS, "r") as f:
        primer3_defaults = yaml.safe_load(f)
        for k, v in primer3_defaults.items():
            if k not in primer3_configs.keys():  # overwrite only if not set
                primer3_configs[k] = v

    # Init a VCF object if a file is provided
    if args.vcf == "":
        vcf_obj = None
    else:
        logger.info("Loading VCF information ...")
        vcf_obj = lib_vcf.VCF(fp=args.vcf, fp_inc_samples=args.vcf_include, fp_exc_samples=args.vcf_exclude)

    if args.webapp:
        logger.info("nanobar - {:d}/{:d}".format(0, 100))

    # Get target sequence
    logger.info("Retrieving target sequence ...")
    roi = _lib_pcr.ROI(
        roi=args.roi,
        blast_db=blast_db,
    )
    roi.fetch_roi()

    # Create a mock marker object (center of the considered region)
    markers = _lib_markers.Markers(
        blast_db=blast_db,
        MARKER_FLANKING_N=MARKER_FLANKING_N,
        MIN_ALIGN_LEN=MIN_ALIGN_LEN,
        MIN_ALIGN_ID=MIN_ALIGN_ID,
        HOMOLOG_FLANKING_N=int(len(roi.seq)/2),
    )

    seq_center = int((roi.end - roi.start) / 2)
    markers.markers = [
        _lib_markers.Marker(
            chrom=roi.chrom,
            pos=roi.start + seq_center - 1,  # to 0-based indexing
            ref_allele=roi.seq[seq_center],
            alt_allele=roi.seq[seq_center],
            name="target",
            HOMOLOG_FLANKING_N=markers.HOMOLOG_FLANKING_N,
        )
    ]

    if len(roi.seq) > 2001:
        # Skip homolog search
        logger.info("Region is > 2001 nt, skipping the homology search ...")
        roi.fasta_name = "{}_{}_{}".format(roi.chrom, roi.start, roi.end)
        seqs = {roi.fasta_name: roi.seq}
        lib_blast.write_fasta(seqs, fp_out=join(blast_db.temporary, "target.fa"))

    else:
        # Find homologs
        logger.info("Finding homeologs/duplications by sequence homology ...")
        seqs = markers.get_marker_flanks()
        markers.find_homologs(seqs)

        if args.webapp:
            logger.info("nanobar - {:d}/{:d}".format(33, 100))

        # Merge markers with the roi object
        roi.fasta_name = markers.markers[0].fasta_name

    # Upload VCF information
    if vcf_obj is not None:
        logger.info("Uploading the VCF file ...")
        roi.upload_mutations(vcf_obj)  # Upload mutations for each markers from a VCF file (if provided)

        # Write subjects containing alternative alleles for each mutations
        if args.report_alts:
            args.report_alts = join(out_path, args.output + "_altlist.txt")
            logger.info("Writing subjects with alternative alleles -> {}".format(args.report_alts))

            if os.path.exists(args.report_alts):
                os.remove(args.report_alts)

            roi.print_alt_subjects(
                vcf_obj=vcf_obj,
                fp=args.report_alts,
            )

    # Start the search for candidate KASP primers
    logger.info("Designing primers - Region size: {} nts".format(len(roi.seq)))
    if not args.webapp:
        logger.info("Using {} parallel processes ...".format(args.n_tasks))
    else:
        logger.info("nanobar - {:d}/{:d}".format(50, 100))

    # Purge sequences from BlastDB
    blast_db.purge()

    # Kwargs
    fp_out = join(out_path, args.output + ".txt")
    kwarg_dict = {
        "fp_fasta": join(blast_db.temporary, "target.fa"),
        "roi": roi,
        "blast_db": blast_db,
        "muscle": muscle,
        "n_primers": args.n_primers,
        "p3_search_depth": args.depth,
        "debug": args.debug,
        "fp_out": fp_out,
        "tm_delta": args.tm_delta,
        "offtarget_size": [args.offtarget_min_size, args.offtarget_max_size],
        "primer_seed": args.seed,
        "primer3_configs": primer3_configs,
    }

    _lib_pcr.main(kwarg_dict)

    if not args.debug:
        shutil.rmtree(temp_path)

    if args.webapp:
        logger.info("nanobar - {:d}/{:d}".format(100, 100))

    if not args.webapp:
        logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
        logger.info("Report written to -> {}".format(fp_out))
    else:
        logger.info("Report ready !")


if __name__ == "__main__":
    main()
