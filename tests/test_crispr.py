import sys
import os
import yaml

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_crispr

with open("tests/KWARGS.yaml", "r") as f:
    KWARGS = yaml.safe_load(f)

# STD
cli_crispr.main(strcmd=" ".join([
    "polyoligo-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
]))

# WEBAPP
cli_crispr.main(strcmd=" ".join([
    "polyoligo-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
    "--webapp",
]))