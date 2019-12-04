import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_crispr

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

cli_crispr.main(strcmd=" ".join([
    "polyoligo-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
    "--debug",
]))

cli_crispr.main(strcmd=" ".join([
    "polyoligo-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["ref_fasta"],
]))


cli_crispr.main(strcmd=" ".join([
    "polyoligo-caps",
    KWARGS["roi_tomato"],
    KWARGS["out"],
    KWARGS["ref_tomato"],
    "--webapp",
    "--debug",
]))

cli_crispr.main(strcmd=" ".join([
    "polyoligo-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
]))