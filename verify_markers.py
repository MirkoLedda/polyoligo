from __future__ import print_function, division
import sys
import argparse
import logging
import os
from os.path import join, exists
import multiprocessing as mp
from copy import deepcopy
# noinspection PyPackageRequirements
from Bio.Seq import Seq
import shutil
import pandas as pd

from src.polyoligo import lib_blast, _lib_kasp, _lib_markers, _logger_config, lib_utils, lib_vcf

BINARIES = {
    "macosx": join(os.path.dirname(__file__), "src/polyoligo", "bin/macosx_x64"),
    "linux": join(os.path.dirname(__file__), "src/polyoligo", "bin/linux_x64"),
    "win": join(os.path.dirname(__file__), "src/polyoligo", "bin/win_x64"),
}
NUC = ["A", "T", "G", "C"]


class Markers(_lib_markers.Markers):
    def __init__(self, blast_db=None, **kwargs):
        super().__init__(blast_db, **kwargs)

    def read_markers(self, fp, logger):
        """Reads a text file of input markers and returns a list of Marker objects.

        File should be tab delimited and no headers and contain the information columns:
        CHR POS NAME REF ALT

        """

        df = pd.read_csv(fp, delim_whitespace=True, header=None)
        df.columns = ["chr", "pos", "name", "ref", "alt"]

        self.markers = []
        for _, row in df.iterrows():

            # Make sure the SNP is biallelic
            if (row["ref"] not in NUC) or (row["alt"] not in NUC):
                logger.warning("Marker {} is not a biallelic SNP [{}]/[{}]".format(
                    row["name"],
                    row["ref"],
                    row["alt"],
                ))
                continue

            self.markers.append(_lib_markers.Marker(
                chrom=row["chr"],
                pos=row["pos"] - 1,  # to 0-based indexing
                ref_allele=row["ref"],
                alt_allele=row["alt"],
                name=row["name"],
            ))


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="check_input_markers",
                                     description="",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

    args = parser.parse_args(inputargs)  # Parse input args

    return args


def main(strcmd=None):
    # Input arguments handling
    if strcmd:  # Means we are running the script using a string of arguments (e.g. for testing)
        testcmd = lib_utils.absolute_paths(strcmd)  # Make paths absolute
        args = parse_args(testcmd.split()[1:])
    else:
        args = parse_args(sys.argv[1:])

    # Set the number of CPUs
    n_tasks = mp.cpu_count() - 1

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

    # Init BlastDB
    blast_db = lib_blast.BlastDB(
        path_db=args.refgenome,
        path_temporary=temp_path,
        path_bin=bin_path,
        job_id="check_markers",
        n_cpu=n_tasks,
    )

    # Make a FASTA reference genome if fast mode is activated
    if args.fast:
        if not blast_db.has_fasta:
            logger.info("Converting the input reference genome to FASTA ...")
            blast_db.db2fasta()

        logger.info("Fast mode On: Loading the input reference genome in memory ...")
        blast_db.load_fasta()

    else:
        if not blast_db.has_db:
            logger.info("Converting the input reference genome to BlastDB ...")
            blast_db.fasta2db()

        blast_db.has_fasta = False
        blast_db.purge()

    # Init Markers object
    markers = Markers(
        blast_db=blast_db,
    )

    # Read markers
    markers.read_markers(args.markers, logger)
    logger.info("Number of target markers = {}".format(len(markers)))

    # Build BLASTDBCMD queries
    queries = []
    for marker in markers.markers:
        q = {"chr": marker.chrom,
             "start": marker.pos1,
             "stop": marker.pos1,
             "strand": "plus",
             }
        queries.append(deepcopy(q))

    # Retrieve the expected nucleotide locations of the markers in the BLASTDB
    blast_nucs = blast_db.fetch(queries)

    # Loop across markers and check if nucleotides match and do the necessary change if not
    fp_out = join(out_path, args.output)  # Write the cleaned markers to a new file
    with open(fp_out, "w") as f:
        for marker in markers.markers:
            q_key = "{}:{}-{}".format(marker.chrom, marker.pos1, marker.pos1)
            bn = blast_nucs[q_key]

            if marker.ref == bn:
                pass
            else:
                # Nucleotides do not match so look for a solution
                n1 = marker.ref
                n2 = marker.alt
                n1c = str(Seq(n1).complement())
                n2c = str(Seq(n2).complement())

                if marker.alt == bn:  # Check for an inversion
                    marker.ref = n2
                    marker.alt = n1
                elif n1c == bn:  # Check for complementarity
                    marker.ref = n1c
                    marker.alt = n2c
                elif n2c == bn:  # Check for reverse complementarity
                    marker.ref = n2c
                    marker.alt = n1c
                else:
                    logger.warning("Marker {}: {}/{} | The reference is {}".format(
                        marker.name,
                        n1,
                        n2,
                        bn,
                    ))
                    continue

                logger.info("Marker {}: {}/{} -> {}/{}".format(
                    marker.name,
                    n1,
                    n2,
                    marker.ref,
                    marker.alt,
                ))

            f.write("{}\t{}\t{}\t{}\t{}\n".format(
                marker.chrom,
                marker.pos1,
                marker.name,
                marker.ref,
                marker.alt,
            ))

    logger.info("Output written to -> {}".format(fp_out))

    shutil.rmtree(temp_path)


if __name__ == "__main__":
    main()
