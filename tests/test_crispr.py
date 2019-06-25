import sys
import os
import shutil
from os.path import join

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_crispr

KWARGS = {
    "roi": "Fvb2-4 10000 11000",
    "out": "out",
    "reference": "sample_data/blastdb",
    "nt": 4,
}
TEMP_DIR = join(os.getcwd(), "temporary")

cli_crispr.main(strcmd=" ".join([
    "cg-crispr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
    # "--pam {}".format(KWARGS["pam"]),
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
]))

# Cleanup
# os.remove(KWARGS["out"] + ".txt")
# os.remove(KWARGS["out"] + ".log")
# os.remove(KWARGS["out"] + "_altlist.txt")
# shutil.rmtree(TEMP_DIR)
