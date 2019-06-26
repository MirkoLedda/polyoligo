import json
import re
from Bio import Seq
from itertools import product


REBASE = "src/polyoligo/data/type2.txt"
OUT_JSON = "src/polyoligo/data/type2.json"



def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input (Credit: https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence)"""
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

    seq2enzyme = {}
    seq_list = []
    pattern = re.compile("[a-zA-Z]+")

    for enzyme in enzymes_clean:
        re_seq = enzyme["recognition"]
        re_seq = re_seq.replace("^", "")
        re_seq = re.findall(pattern, re_seq)[0]
        if re_seq.count("N") <= 4:
            re_seq = re_seq.replace("^", "")
            re_seq = re.findall(pattern, re_seq)[0]
            re_seqs = extend_ambiguous_dna(re_seq)

            seq_list += re_seqs

            for re_seq in re_seqs:
                if re_seq in seq2enzyme.keys():
                    seq2enzyme[re_seq].append(enzyme["name"])
                else:
                    seq2enzyme[re_seq] = [enzyme["name"]]

    print(len(seq_list))


def main():
    rebase2json(fp=REBASE)


if __name__ == "__main__":
    main()
