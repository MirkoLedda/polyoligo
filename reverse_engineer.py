import os
from os.path import join, exists
import argparse
import sys
import multiprocessing as mp
import logging
import pandas as pd
from copy import deepcopy
import tqdm
import shutil
import numpy as np

# noinspection PyProtectedMember
from src.polyoligo import cli_kasp, blast_lib, _getkasp, _markers_lib, utils, vcf_lib, _logger_config

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}

HOMOLOG_FLANKING_N = cli_kasp.HOMOLOG_FLANKING_N


class Marker(_markers_lib.Marker):
    def __init__(self, chrom, pos, ref_allele, alt_allele, name=None):
        super().__init__(chrom, pos, ref_allele, alt_allele, name)
        self.pp = None

    def add_primers(self, pp_dir, starts, stops, seqs):
        if pp_dir == "F":
            self.start = starts["REF"]
            self.stop = stops["COM"]
        else:
            self.start = starts["COM"]
            self.stop = stops["REF"]

        self.pp = _getkasp.PrimerPair()

        for d in ["F", "R", "A"]:
            if d == "A":
                ptype = "ALT"
            elif pp_dir == d:
                ptype = "REF"
            else:
                ptype = "COM"

            self.pp.primers[d].start = starts[ptype]
            self.pp.primers[d].stop = stops[ptype]
            self.pp.primers[d].sequence = seqs[ptype]

        self.pp.chrom = self.chrom
        self.pp.dir = pp_dir
        self.pp.product_size = self.stop - self.start + 1

    def score_primers(self, fp_reporters):
        pcr = _getkasp.PCR(self.name, self.chrom, self.pos, self.ref, self.alt)
        pcr.pps = [self.pp]

        pcr.add_mutations(self.mutations)
        pcr.add_reporter_dyes(fp_reporters)
        pcr.check_heterodimerization()
        pcr.classify()
        pcr.prune(n=np.inf)

        return pcr


class Markers(_markers_lib.Markers):
    def __init__(self, blast_db, **kwargs):
        super().__init__(blast_db, **kwargs)
        self.markers = []

    def read_markers(self, fp):
        df = pd.read_csv(fp, delim_whitespace=True, header=None)
        df.columns = ["name", "chr", "pos", "ref", "alt", "start", "stop", "dir", "type", "id", "seq"]

        marker_names = []
        starts = {}
        stops = {}
        seqs = {}
        marker = None

        # Each markers spans 3 rows for REF/ALT/COM primers
        collected_primers = []
        for i in range(len(df)):

            if pd.isna(df["type"][i]):
                continue

            if df["type"][i] == "REF":
                collected_primers.append("REF")
                starts = {}
                stops = {}
                seqs = {}

                marker = Marker(
                    chrom=df["chr"][i],
                    pos=df["pos"][i],
                    ref_allele=df["ref"][i],
                    alt_allele=df["alt"][i],
                    name=df["name"][i],
                )

            if df["type"][i] == "ALT":
                collected_primers.append("ALT")
                pp_dir = df["dir"][i]

            if df["type"][i] == "COM":
                collected_primers.append("COM")

            starts[df["type"][i]] = df["start"][i]
            stops[df["type"][i]] = df["stop"][i]
            seqs[df["type"][i]] = df["seq"][i].upper()

            if len(collected_primers) == 3:
                # Check for duplicate names in the markers
                j = 0
                n = marker.name
                while True:
                    if n in marker_names:
                        n = marker.name + "_{:d}".format(j)
                        j += 1
                    else:
                        marker.name = n
                        marker_names.append(marker.name)
                        break
                marker.add_primers(pp_dir, starts, stops, seqs)
                self.markers.append(deepcopy(marker))

                collected_primers = []


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
        # help="Maximum number of KASP primers to return for each marker.",
        help=argparse.SUPPRESS,
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
        "--reporters",
        metavar="<FASTA>",
        type=str,
        default="",
        help="Sequence of two reporter dyes named 'dye1' and 'dye2' in FASTA format and in a 5'-3' orientation",
    )
    parser.add_argument(
        "--multiplier",
        metavar="<FLOAT>",
        type=float,
        default=2.0,
        # help="This parameter controls the exhaustiveness of the primer pair search, which is given by "
        #      "'n' * 'multiplier'. By increasing this value, more primer pairs will be considered but the process
        #      will "
        #      "be computationally heavier.",
        help=argparse.SUPPRESS,
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
        "--debug",
        action="store_true",
        # help="Debug mode. FOR DEVS ONLY.",
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args(inputargs)  # Parse input args

    return args


def score_marker(kwargs_dict):
    marker = kwargs_dict["marker"]
    pcr = marker.score_primers(kwargs_dict["fp_reporters"])
    _getkasp.print_report(pcr, kwargs_dict["fp_out"])


def main():
    main_time = utils.timer_start()  # Set main timer
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
    curr_os = utils.get_os()
    if curr_os is None:
        logger.error("OS not supported or not detected properly.")
        sys.exit()
    bin_path = BINARIES[curr_os]

    # Init the BlastDB
    blast_db = blast_lib.BlastDB(
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

    # Init Markers object
    markers = Markers(
        blast_db=blast_db
    )

    # Read markers
    markers.read_markers(args.markers)
    logger.info("Number of target markers = {}".format(len(markers)))

    # Init a VCF object if a file is provided
    if args.vcf == "":
        vcf_obj = None
    else:
        logger.info("Loading VCF information ...")
        vcf_obj = vcf_lib.VCF(fp=args.vcf, fp_inc_samples=args.vcf_include, fp_exc_samples=args.vcf_exclude)

    markers.upload_mutations(vcf_obj)  # Upload mutations for each markers from a VCF file (if provided)

    # Start the search for candidate KASP primers
    logger.info("Scoring KASP candidates using {} parallel processes ...".format(args.n_tasks))

    # Initialize a list of kwargs for the worker and a job queue
    kwargs_worker = []

    for marker in markers.markers:
        blast_db.job_id = marker.name  # To make sure files are not overwritten during multiprocessing
        blast_db.n_cpu = 1  # Blast should now only use 1 thread as we will be spawning jobs

        kwarg_dict = {
            "marker": marker,
            "fp_out": join(temp_path, marker.name + ".txt"),
            "fp_reporters": args.reporters,
        }
        kwargs_worker.append(deepcopy(kwarg_dict))

    # Run the job queue
    n_jobs = len(kwargs_worker)
    with mp.Pool(processes=args.n_tasks, maxtasksperchild=1) as p:
        _ = list(tqdm.tqdm(p.imap_unordered(score_marker, kwargs_worker), total=n_jobs, disable=args.silent))

    # Concatenate all primers for all markers into a single report
    logger.info("Preparing report ...")
    fp_out = join(out_path, args.output + ".txt")
    markers.write_report(fp_out)

    if not args.debug:
        shutil.rmtree(temp_path)

    logger.info("Total time elapsed: {}".format(utils.timer_stop(main_time)))
    logger.info("Report written to -> {}".format(fp_out))


if __name__ == "__main__":
    main()
