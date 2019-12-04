import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_hrm

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

cli_hrm.main(strcmd=" ".join([
    "polyoligo-hrm",
    KWARGS["marker"],
    KWARGS["out"],
    KWARGS["reference"],
    "--webapp",
    "--debug",
]))

cli_hrm.main(strcmd=" ".join([
    "polyoligo-hrm",
    KWARGS["marker"],
    KWARGS["out"],
    KWARGS["ref_fasta_sample"],
]))

cli_hrm.main(strcmd=" ".join([
    "polyoligo-hrm",
    KWARGS["marker_indels"],
    KWARGS["out"],
    KWARGS["reference"],
    "--webapp",
    "--debug",
]))

cli_hrm.main(strcmd=" ".join([
    "polyoligo-hrm",
    KWARGS["marker_tomato"],
    KWARGS["out"],
    KWARGS["ref_tomato"],
    "--webapp",
    "--debug"
]))

cli_hrm.main(strcmd=" ".join([
    "polyoligo-hrm",
    KWARGS["marker_lim"],
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
    "--primer3 {}".format(KWARGS["primer3_hrm"]),
    "-nt {}".format(KWARGS["nt"]),
]))
