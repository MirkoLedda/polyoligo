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
import yaml

# noinspection PyProtectedMember
from src.polyoligo import cli_kasp, lib_blast, _lib_kasp, _lib_markers, lib_utils, lib_vcf, _logger_config, lib_primer3

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}

HOMOLOG_FLANKING_N = cli_kasp.HOMOLOG_FLANKING_N
PRIMER3_DEFAULTS = join(os.path.dirname(__file__), "src/polyoligo/data/PRIMER3_KASP.yaml")


class Marker(_lib_markers.Marker):
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

        self.pp = _lib_kasp.PrimerPair()

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

    def score_primers(self, fp_reporters, blast_db):
        pcr = _lib_kasp.PCR(self.name, self.chrom, self.pos, self.ref, self.alt)
        pcr.pps = [self.pp]

        pcr.check_offtargeting(blast_db)
        pcr.add_mutations(self.mutations)
        pcr.add_reporter_dyes(fp_reporters)
        pcr.check_heterodimerization()
        pcr.classify()
        pcr.prune(n=np.inf)

        return pcr


class Markers(_lib_markers.Markers):
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

            seq = df["seq"][i]
            seqp = ""
            flag_get = False
            for c in seq:
                if c.islower():
                    flag_get = True

                if flag_get:
                    seqp += c

            seqs[df["type"][i]] = seqp.upper()

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
    if not blast_db.has_fasta:
        logger.info("Converting the input reference genome to FASTA ...")
        blast_db.db2fasta()

    blast_db.load_fasta()

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
        vcf_obj = lib_vcf.VCF(fp=args.vcf, fp_inc_samples=args.vcf_include, fp_exc_samples=args.vcf_exclude)

    markers.upload_mutations(vcf_obj)  # Upload mutations for each markers from a VCF file (if provided)

    # Purge sequences from BlastDB
    blast_db.purge()

    # Start the search for candidate KASP primers
    logger.info("Scoring KASP candidates using {} parallel processes ...".format(args.n_tasks))

    # Initialize a list of kwargs for the worker and a job queue
    kwargs_worker = []

    for marker in markers.markers:
        blast_db.job_id = marker.name  # To make sure files are not overwritten during multiprocessing
        blast_db.n_cpu = 1  # Blast should now only use 1 thread as we will be spawning jobs

        kwarg_dict = {
            "marker": marker,
            "fp_out": join(temp_path, marker.name),
            "fp_reporters": reporters,
            "blast_db": blast_db,
        }
        kwargs_worker.append(deepcopy(kwarg_dict))

    # Run the job queue
    n_jobs = len(kwargs_worker)
    with mp.Pool(processes=args.n_tasks, maxtasksperchild=1) as p:
        _ = list(tqdm.tqdm(p.imap_unordered(score_marker, kwargs_worker), total=n_jobs, disable=args.silent))

    # Concatenate all primers for all markers into a single report
    logger.info("Preparing report ...")
    markers.write_kasp_reports(join(out_path, args.output))

    if not args.debug:
        shutil.rmtree(temp_path)

    logger.info("Total time elapsed: {}".format(lib_utils.timer_stop(main_time)))
    logger.info("Report written to -> {} [{}, {}]".format(join(out_path, args.output) + ".txt", ".bed", ".log"))


if __name__ == "__main__":
    main()
