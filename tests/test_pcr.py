import sys
import os
import shutil
from os.path import join

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_pcr

KWARGS = {
    "roi": "Fvb2-4:100000-101000",
    "out": "out",
    "reference": "sample_data/blastdb",
    "vcf": "sample_data/vcf.txt.gz",
    "vcf_include": "sample_data/vcf_include.txt",
    "vcf_exclude": "sample_data/vcf_include.txt",
    "nt": 4,
    "tm": 50,
    "seed": 10,
}
TEMP_DIR = join(os.getcwd(), "temporary")

cli_pcr.main(strcmd=" ".join([
    "cg-pcr",
    KWARGS["roi"],
    KWARGS["out"],
    KWARGS["reference"],
    # "--pam {}".format(KWARGS["pam"]),
    "-nt {}".format(KWARGS["nt"]),
    "--debug",
    "--primer3 sample_data/primer3_example.yaml",
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_include {}".format(KWARGS["vcf_include"]),
    "--seed {}".format(KWARGS["seed"])
]))

# Cleanup
# os.remove(KWARGS["out"] + ".txt")
# os.remove(KWARGS["out"] + ".log")
# os.remove(KWARGS["out"] + "_altlist.txt")
# shutil.rmtree(TEMP_DIR)
