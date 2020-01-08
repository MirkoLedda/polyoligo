import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_kasp

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker_pistachio"],
    KWARGS["out"],
    KWARGS["ref_pistachio"],
    "--vcf {}".format(KWARGS["vcf_pistachio"]),
    "--report_alts",
    "--debug",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker"],
    KWARGS["out"],
    KWARGS["reference"],
    "--report_alts",
    "--webapp",
    "--debug",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker_indels"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n 50",
    "--depth 1",
    "--webapp",
    "--debug",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker"],
    KWARGS["out"],
    KWARGS["ref_fasta"],
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker_tomato"],
    KWARGS["out"],
    KWARGS["ref_tomato"],
    "--debug",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker_cotton"],
    KWARGS["out"],
    KWARGS["ref_cotton"],
    "--debug",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo-kasp",
    KWARGS["marker_lim"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n {}".format(KWARGS["n"]),
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_include {}".format(KWARGS["vcf_include"]),
    "--dye1 {}".format(KWARGS["dye1"]),
    "--dye2 {}".format(KWARGS["dye2"]),
    "--depth {}".format(KWARGS["depth"]),
    "--tm_delta {}".format(KWARGS["tm_delta"]),
    "--seed {}".format(KWARGS["seed"]),
    "--offtarget_min_size {}".format(KWARGS["offtarget_min_size"]),
    "--offtarget_max_size {}".format(KWARGS["offtarget_max_size"]),
    "--primer3 {}".format(KWARGS["primer3_kasp"]),
    "-nt {}".format(KWARGS["nt"]),
]))
