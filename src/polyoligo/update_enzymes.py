import json
import re
# noinspection PyPackageRequirements
from Bio import Seq
from itertools import product

REBASE = "src/polyoligo/data/type2.txt"
OUT_JSON = "src/polyoligo/data/type2.json"


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input
    (Credit: https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence)"""
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def rebase2json(fp):
    skipped_header = False
    enzymes = []

    with open(fp, "r") as f:

        for line in f:
            line = line.strip()
            if line.startswith("<1>"):
                if skipped_header:
                    enzyme = {
                        "name": line.split("<1>")[1],
                        "prototype": next(f).strip().split("<2>")[1],
                        "recognition": next(f).strip().split("<3>")[1],
                        "methylation": next(f).strip().split("<4>")[1],
                        "supplier": next(f).strip().split("<5>")[1],
                    }
                    enzymes.append(enzyme)
                else:
                    skipped_header = True

    # Filter enzymes
    enzymes_clean = []
    for enzyme in enzymes:
        if (enzyme["supplier"] != "") and (enzyme["prototype"] == ""):
            enzymes_clean.append(enzyme)

    enzymes = []
    max_n = 0
    for enzyme in enzymes_clean:
        re_seq = enzyme["recognition"]
        re_seq = re_seq.replace("^", "")
        re_seq = re.findall(re.compile("[a-zA-Z]+"), re_seq)[0]
        re_seq = re_seq.replace("^", "")
        re_seq_rc = Seq.reverse_complement(re_seq)
        enzyme["n"] = len(re_seq)
        max_n = max([max_n, enzyme["n"]])

        # build the REGEX for each recognition site (both forward and reverse)
        enzyme["regex"] = {}

        # Forward regex
        regex = ""
        n_non_amb = 0
        for nuc in re_seq:
            if nuc in ["A", "T", "G", "C"]:
                regex += nuc
                n_non_amb += 1
            else:
                nuc = Seq.IUPAC.IUPACData.ambiguous_dna_values[nuc]
                regex += "[{}]".format(nuc)
        enzyme["regex"]["F"] = regex
        enzyme["n_certain_nuc"] = n_non_amb

        # Reverse regex
        regex = ""
        for nuc in re_seq_rc:
            if nuc in ["A", "T", "G", "C"]:
                regex += nuc
            else:
                nuc = Seq.IUPAC.IUPACData.ambiguous_dna_values[nuc]
                regex += "[{}]".format(nuc)
        enzyme["regex"]["R"] = regex

        enzymes.append(enzyme)

        with open(OUT_JSON, "w") as f:
            json.dump(enzymes, f)

    print("Largest recognition site: {}".format(max_n))


def main():
    rebase2json(fp=REBASE)


if __name__ == "__main__":
    main()
