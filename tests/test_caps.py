import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_caps

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

cli_caps.main(strcmd=" ".join([
    "polyoligo-caps",
    KWARGS["marker"],
    KWARGS["out"],
    KWARGS["ref_sample"],
    "--webapp",
    "--debug",
]))

cli_caps.main(strcmd=" ".join([
    "polyoligo-caps",
    KWARGS["marker_tomato"],
    KWARGS["out"],
    KWARGS["ref_tomato"],
    "--webapp",
    "--debug"
]))

cli_caps.main(strcmd=" ".join([
    "polyoligo-caps",
    KWARGS["marker_lim"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n {}".format(KWARGS["n"]),
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_include {}".format(KWARGS["vcf_include"]),
    "--report_alts",
    "--depth {}".format(KWARGS["depth"]),
    "--enzymes {}".format(KWARGS["enzymes"]),
    "--fragment_min_size {}".format(KWARGS["fragment_min_size"]),
    "--tm_delta {}".format(KWARGS["tm_delta"]),
    "--seed {}".format(KWARGS["seed"]),
    "--offtarget_min_size {}".format(KWARGS["offtarget_min_size"]),
    "--offtarget_max_size {}".format(KWARGS["offtarget_max_size"]),
    "--primer3 {}".format(KWARGS["primer3_caps"]),
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
]))
