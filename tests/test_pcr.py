import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_pcr

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

cli_pcr.main(strcmd=" ".join([
    "polyoligo-pcr",
    KWARGS["roi_pcr"],
    KWARGS["out"],
    KWARGS["reference"],
    "--depth 100",
    "--seed 14",
    "--webapp",
    "--debug",
]))

cli_pcr.main(strcmd=" ".join([
    "polyoligo-pcr",
    KWARGS["roi_pcr"],
    KWARGS["out"],
    KWARGS["ref_fasta"],
]))

cli_pcr.main(strcmd=" ".join([
    "polyoligo-pcr",
    KWARGS["roi_pcr_tomato"],
    KWARGS["out"],
    KWARGS["ref_tomato"],
    "--webapp",
    "--debug",
]))

cli_pcr.main(strcmd=" ".join([
    "polyoligo-pcr",
    KWARGS["roi_pcr_lim"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n {}".format(KWARGS["n"]),
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_include {}".format(KWARGS["vcf_include"]),
    "--depth {}".format(KWARGS["depth"]),
    "--tm_delta {}".format(KWARGS["tm_delta"]),
    "--seed {}".format(KWARGS["seed"]),
    "--offtarget_min_size {}".format(KWARGS["offtarget_min_size"]),
    "--offtarget_max_size {}".format(KWARGS["offtarget_max_size"]),
    "--primer3 {}".format(KWARGS["primer3_pcr"]),
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
]))
