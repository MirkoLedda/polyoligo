from __future__ import print_function, division
import sys
import argparse
import logging
import os
from os.path import join, exists
# noinspection PyPackageRequirements
from Bio.Seq import Seq
import shutil
import pandas as pd

from src.polyoligo import lib_blast, lib_markers, _logger_config, lib_utils

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}
NUC = ["A", "T", "G", "C"]


class Marker(lib_markers.Marker):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def assert_marker(self, blast_hook):
        if self.ref != "":
            query = [{
                "chr": self.chrom,
                "start": int(self.pos1),
                "stop": int(self.pos1 + self.n - 1),
            }]

            seq = "N"
            for _, seq in blast_hook.fetch(query).items():
                seq = seq.upper()
                break

            if self.ref == seq:
                pass
            else:
                # Nucleotides do not match so look for a solution
                n1 = self.ref
                n2 = self.alt
                n1c = str(Seq(n1).complement())
                n2c = str(Seq(n2).complement())

                if self.alt == seq:  # Check for an inversion
                    self.ref = n2
                    self.alt = n1
                elif n1c == seq:  # Check for complementarity
                    self.ref = n1c
                    self.alt = n2c
                elif n2c == seq:  # Check for reverse complementarity
                    self.ref = n2c
                    self.alt = n1c
                else:
                    print("{} {}/{} -> ERROR skipped as reference is {}".format(
                        self.name,
                        n1,
                        n2,
                        seq,
                    ))
                    return None

                print("{} {}/{} -> {}/{}".format(
                    self.name,
                    n1,
                    n2,
                    self.ref,
                    self.alt,
                ))


def read_markers(fp):
    """Reads a text file of input markers and returns a list of Marker objects.

    File should be tab delimited and no headers and contain the information columns:
    CHR POS NAME REF ALT

    """

    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    df.columns = ["chr", "pos", "name", "ref", "alt"]

    markers = []
    for _, row in df.iterrows():
        markers.append(Marker(
            chrom=row["chr"],
            pos=row["pos"] - 1,  # to 0-based indexing
            ref_allele=row["ref"],
            alt_allele=row["alt"],
            name=row["name"],
        ))

    return markers


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="check_input_markers",
                                     description="",
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
        help="File containing a list of input markers in tab-separated format and with no header as: "
             "CHR POS NAME REF ALT",
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        type=str,
        help="Output marker filename.",
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

    args = parser.parse_args(inputargs)  # Parse input args

    return args


def main(strcmd=None):
    # Input arguments handling
    if strcmd:  # Means we are running the script using a string of arguments (e.g. for testing)
        testcmd = lib_utils.absolute_paths(strcmd)  # Make paths absolute
        args = parse_args(testcmd.split()[1:])
    else:
        args = parse_args(sys.argv[1:])

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

    # Initialize hooks
    blast_hook = lib_blast.BlastDB(
        path_db=args.refgenome,
        path_temporary=temp_path,
        path_bin=bin_path,
        job_id="main",
        n_cpu=1,
    )

    # Build the BlastDB if a fasta was provided
    if not blast_hook.has_db:
        logger.info("Converting the input reference genome to BlastDB ...")
        blast_hook.fasta2db()

    # Read markers
    try:
        markers = read_markers(args.markers)
    except:
        logger.error("Failed to read input markers. Please check the input file path or format.")
        sys.exit(1)

    fp_out = join(out_path, args.output)  # Write the cleaned markers to a new file
    with open(fp_out, "w") as f:
        for marker in markers:
            marker.assert_marker(blast_hook)
            f.write("{}\t{}\t{}\t{}\t{}\n".format(
                marker.chrom,
                marker.pos1,
                marker.name,
                marker.ref,
                marker.alt,
            ))

    shutil.rmtree(temp_path)
    os.remove(fp_out + ".log")  # Remove the log, this is redundant here


if __name__ == "__main__":
    main()
