from __future__ import print_function, division
import sys
import argparse
import logging
import os
from os.path import join, exists
import pandas as pd

# noinspection PyProtectedMember
from src.polyoligo import _logger_config, lib_utils, lib_vcf


def parse_args(inputargs):
    # Define the args parser
    parser = argparse.ArgumentParser(prog="vcf2markers",
                                     description="",
                                     epilog="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "vcf",
        metavar="VCF",
        type=str,
        help="Input VCF file.",
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        type=str,
        help="Output filename for returning markers.",
    )
    parser.add_argument(
        "roi",
        metavar="ROI",
        type=str,
        help="File listing regions of interest with with no header as: CHR START END. START and END are optional, "
             "if not provided the entire chromosome will be considered",
    )
    args = parser.parse_args(inputargs)  # Parse input args

    return args


def read_roi(fp):
    df = pd.read_csv(fp, delim_whitespace=True, header=None)
    return df


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
    args.output = join(out_path, os.path.basename(args.output))

    if not exists(out_path):
        os.makedirs(out_path)

    if exists(args.output):
        os.remove(args.output)

    # Init the logger
    _logger_config.setup_console_logging(verbose=True)
    logger = logging.getLogger(__name__)

    # Init the VCF object
    vcf_obj = lib_vcf.VCF(fp=args.vcf)

    # Read ROIs
    rois = read_roi(args.roi)

    # Loop across ROIs
    for _, roi in rois.iterrows():
        if len(roi) == 1:
            # fetch VCF based on chromosomes only and write to file
            vcf_obj.fetch_kasp_markers(fp=args.output, chrom=roi[0])
        else:
            # fetch VCF and write to file
            vcf_obj.fetch_kasp_markers(fp=args.output, chrom=roi[0], start=roi[1], stop=roi[2])

    logger.info("Output written to -> {}".format(args.output))


if __name__ == "__main__":
    main()