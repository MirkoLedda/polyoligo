import sys
import os
import cProfile

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_kasp

KWARGS = {
    "markers": "sample_data/markers.txt",
    "out": "out",
    "reference": "sample_data/blastdb",
    # "reference": "sample_data/blastdb.fa.gz",
    "n": 10,
    "vcf": "sample_data/vcf.txt.gz",
    "vcf_include": "sample_data/vcf_include.txt",
    "vcf_exclude": None,
    "reporters": "sample_data/VIC_FAM_reporters.txt",
    "multiplier": 2,
    "tm_delta": 10,
    "nt": 4,
    "primer3": "sample_data/primer3_configs.yaml",
}


def cprofile_job():
    cProfile.run('main()', 'profile.out')


def main():
    cli_kasp.main(strcmd=" ".join([
        "polyoligo",
        KWARGS["markers"],
        KWARGS["out"],
        KWARGS["reference"],
        "--vcf {}".format(KWARGS["vcf"]),
        "--vcf_include {}".format(KWARGS["vcf_include"]),
        "--reporters {}".format(KWARGS["reporters"]),
        "--multiplier {}".format(KWARGS["multiplier"]),
        "--tm_delta {}".format(KWARGS["tm_delta"]),
        "-nt {}".format(KWARGS["nt"]),
        "-n {}".format(KWARGS["n"]),
        "--debug",
        "--report_alts",
        "--primer3 {}".format(KWARGS["primer3"]),
    ]))


if __name__ == "__main__":
    cprofile_job()
