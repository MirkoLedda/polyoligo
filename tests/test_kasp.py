import sys
import os
import shutil
from os.path import join

sys.path.insert(0, os.path.abspath("."))

from src.polyoligo import cli_kasp, lib_utils, lib_vcf

KWARGS = {
    "markers": "sample_data/markers.txt",
    "out": "out",
    "reference": "sample_data/blastdb",
    "n": 100,
    "vcf": "sample_data/vcf.txt.gz",
    "vcf_include": "sample_data/vcf_include.txt",
    "vcf_exclude": "sample_data/vcf_include.txt",
    "dye1": "GAAGGTCGGAGTCAACGGATT",
    "dye2": "GAAGGTGACCAAGTTCATGCT",
    "multiplier": 1,
    "tm_delta": 10,
    "primer3": "sample_data/primer3_configs.yaml",
    "nt": 4,
    "seed": 14,
}
TEMP_DIR = join(os.getcwd(), "temporary")

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n {}".format(KWARGS["n"]),
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_include {}".format(KWARGS["vcf_include"]),
    "--dye1 {}".format(KWARGS["dye1"]),
    "--dye2 {}".format(KWARGS["dye2"]),
    "--multiplier {}".format(KWARGS["multiplier"]),
    "--tm_delta {}".format(KWARGS["tm_delta"]),
    "-nt {}".format(KWARGS["nt"]),
    "--primer3 {}".format(KWARGS["primer3"]),
    "--report_alts",
    "--debug",
    "--seed {}".format(KWARGS["seed"]),
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n {}".format(KWARGS["n"]),
    "--dye1 {}".format(KWARGS["dye1"]),
    "--dye2 {}".format(KWARGS["dye2"]),
    "--webapp",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "-n 1",
    "--vcf {}".format(KWARGS["vcf"]),
    "--vcf_exclude {}".format(KWARGS["vcf_exclude"]),
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "-nt -10",
]))

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "--fast",
    "-n 1",
]))

KWARGS["reference"] = "sample_data/blastdb.fa.gz"

cli_kasp.main(strcmd=" ".join([
    "polyoligo",
    KWARGS["markers"],
    KWARGS["out"],
    KWARGS["reference"],
    "--fast",
    "--debug",
    "-n 1",
]))

vcf_obj = lib_vcf.VCF(
    fp=KWARGS["vcf"],
)
vcf_obj.fetch_kasp_markers(
    fp=join(TEMP_DIR, "vcf.out"),
    chrom="Fvb2-4",
)

lib_utils.read_json(fp=join(TEMP_DIR, "AX-184205870_blast_F.json"))

lib_utils.write_fasta(
    seqs={
        "A": "A",
        "B": "A",
    },
    fp_out=join(TEMP_DIR, "test.fa"),
)

lib_utils.read_fasta(
    fp=join(TEMP_DIR, "test.fa"),
    filters=["A"],
)

_ = lib_utils.padding('a')
_ = lib_utils.padding_left('a')
_ = lib_utils.padding_right('a')

# Cleanup
os.remove(KWARGS["out"] + ".txt")
os.remove(KWARGS["out"] + ".log")
os.remove(KWARGS["out"] + "_altlist.txt")
shutil.rmtree(TEMP_DIR)
