import numpy as np
import sys
import time
from datetime import timedelta
import os
import json
from copy import deepcopy
import re
# noinspection PyPackageRequirements
from Bio import Seq

IUPAC_DNA_REGEX = "^[atgcrynmkATGCRYNMK]+$"
IUPAC_RNA_REGEX = "^[augcrynmkAUGCRYNMK]+$"

IUPAC_DNA = Seq.IUPAC.IUPACData.ambiguous_dna_values
IUPAC_DNA_I = {}
for k, v in IUPAC_DNA.items():
    IUPAC_DNA_I[v] = k


def round_tidy(x, n):
    """Return 'x' rounded to 'n' significant digits."""
    y = np.abs(x)
    if y <= sys.float_info.min:
        return 0.0
    else:
        return np.round(x, int(n - np.ceil(np.log10(y))))


def seconds_to_hms(t):
    """Formats a datetime.timedelta object to a HH:MM:SS string."""

    hours, remainder = divmod(t.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    hours, minutes, seconds = int(hours), int(minutes), int(seconds)
    if hours < 10:
        hours = "0%s" % int(hours)
    if minutes < 10:
        minutes = "0%s" % minutes
    if seconds < 10:
        seconds = "0%s" % seconds
    return "%s:%s:%s" % (hours, minutes, seconds)


def timer_start():
    """Start a timer."""

    return time.time()


def timer_stop(t0):
    """Stop a timer and return time as a formatted string."""

    t = timedelta(seconds=time.time() - t0)
    return seconds_to_hms(t)


def absolute_paths(txt):
    """Replace paths within a string with absolute paths"""

    txt = txt.replace("~", os.path.expanduser("~"))
    txt = txt.replace(" .", " " + os.path.abspath("."))
    txt = txt.replace(" ..", " " + os.path.abspath(".."))

    return txt


def get_os():
    platform = sys.platform
    cwos = None

    if platform == "darwin":  # Mac OSX
        cwos = "macosx"
    elif platform.startswith("linux"):  # Linux
        cwos = "linux"
    elif (platform == "win32") or (platform == "cygwin"):  # Windows
        cwos = "win"

    return cwos


def read_json(fp):
    with open(fp, "r") as f:
        json_dict = json.load(f)
    return json_dict


def write_fasta(seqs, fp_out):
    with open(fp_out, "w") as f:
        for name, seq in seqs.items():
            f.write(">{}\n{}\n".format(name, seq))


def read_fasta(fp, filters=None):
    seqs = {}

    with open(fp, "r") as f:

        name = None
        seq = ""

        for line in f:
            if line.startswith(">"):
                curr_name = line.split(">")[1].strip()

                if name is not None:
                    if filters is not None:
                        if name in filters:
                            seqs[name] = seq
                    else:
                        seqs[name] = seq

                name = curr_name
                seq = ""
            else:
                seq += line.strip()

        seqs[name] = seq  # Add the last entry

    return seqs


def list_2_intervals(ixs):
    """Take a list of indexes and returns a list of intervals in the form of (start, stop) tuples."""
    ivs = []

    if len(ixs) != 0:
        ixs = np.sort(np.array(ixs))
        delta = np.diff(ixs)
        anchors = np.array(ixs[1:])
        anchors = anchors[delta != 1]
        anchors_loc = np.where(delta != 1)[0]

        i = ixs[0]

        for ia, a in enumerate(anchors):
            ivs.append((i, ixs[anchors_loc[ia]]))
            i = ixs[anchors_loc[ia] + 1]
        ivs.append((i, ixs[-1]))

    return ivs


def list_2_ranges(ixs):
    """Take a list of indexes and returns a list of intervals in the form of [start, length] list."""
    ivs = []

    if len(ixs) != 0:
        ixs = np.sort(np.array(ixs))
        delta = np.diff(ixs)
        anchors = np.array(ixs[1:])

        i = ixs[0]
        n = 1
        for d, a in zip(delta, anchors):
            if d == 1:
                n += 1
            else:
                ivs.append([int(i), int(n)])
                i = a
                n = 1
        ivs.append([int(i), int(n)])

    return ivs


def kwargs2attr_deep(obj, kwargs):
    """Deep copy attributes values based on keyword arguments. Assign only if the attribute already exists in obj."""
    if kwargs is None:
        return

    if not isinstance(kwargs, dict):
        raise TypeError("Attempted to set arguments not using a dictionary.")

    attr = obj.__dict__.keys()
    for key in kwargs.keys():
        if key in attr:
            setattr(obj, key, deepcopy(kwargs[key]))


def padding(x, left=0, right=0, char="N"):
    x = char * left + x + char * right
    return x


def padding_left(x, n=0, char="N"):
    return (n - len(x)) * char + x


def padding_right(x, n=0, char="N"):
    return x + (n - len(x)) * char


def get_n_padding(x, n):
    lp = 0
    rp = 0

    if x < n:
        d = n - x
        if (d % 2) != 0:
            lp = int(np.floor(d / 2))
            rp = int(np.ceil(d / 2))
        else:
            lp = int(d / 2)
            rp = lp
    return lp, rp


def is_dna(x):
    if x == "":
        return True
    else:
        return bool(re.match(IUPAC_DNA_REGEX, x))


def seqs2ambiguous_dna(seqs):
    """Takes a list of same length sequences and returns their ambiguous form"""
    nucs = list(map(set, zip(*seqs)))
    nucs = ["".join(sorted(list(nuc))) for nuc in nucs]

    seq_amb = []
    for nuc in nucs:
        try:
            seq_amb.append(IUPAC_DNA_I[nuc])
        except KeyError:
            seq_amb.append("N")

    seq_amb = "".join(seq_amb)

    return seq_amb


def reduce_ivs(ivs, n):
    """Take a list of intervals and merges them until n survive."""

    if len(ivs) > 0:
        ivs = list(ivs)  # Make a deep copy

        # Compute spaces between intervals
        past_iv = ivs[0]
        spaces = []
        for iv in ivs[1:]:
            i = past_iv[0] + (past_iv[1] - 1)
            j = iv[0]
            spaces.append(j - i - 1)
            past_iv = iv

        while len(ivs) > n:
            # Reduce the interval list by merging the two closest intervals
            i = int(np.argmin(spaces))
            _ = spaces.pop(i)

            iv1 = ivs.pop(i)
            iv2 = ivs.pop(i)

            nstop = iv2[0] + (iv2[1] - 1)
            niv = [iv1[0], nstop - iv1[0] + 1]

            ivs.insert(i, niv)

    return ivs
