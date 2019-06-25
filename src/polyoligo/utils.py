import numpy as np
import sys
import time
from datetime import timedelta
import os
import json
from copy import deepcopy
import re

IUPAC_DNA_REGEX = "^[atgcrynmkATGCRYNMK]+$"
IUPAC_RNA_REGEX = "^[augcrynmkAUGCRYNMK]+$"

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
    k = None

    if platform == "darwin":  # Mac OSX
        k = "macosx"
    elif platform.startswith("linux"):  # Linux
        k = "linux"
    elif (platform == "win32") or (platform == "cygwin"):  # Windows
        k = "win"

    return k


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
    return (n-len(x)) * char + x


def padding_right(x, n=0, char="N"):
    return x + (n-len(x)) * char


def is_dna(x):
    if x == "":
        return True
    else:
        return bool(re.match(IUPAC_DNA_REGEX, x))
