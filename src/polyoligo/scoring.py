import numpy as np


def score_base(pp):
    score = 0
    qcode = ""

    tms = np.array([pp.primers[d].tm for d in pp.primers.keys()])
    tms_l1_norm = np.sum(np.abs(tms - np.mean(tms)))
    if tms_l1_norm <= 5:
        score += 1
    else:
        qcode += "t"

    if len(pp.offtargets) == 0:
        score += 3
    else:
        qcode += "O"

    if not pp.dimer:
        score += 1
    else:
        qcode += "d"

    max_aafs = np.array([pp.primers[d].max_aaf for d in pp.primers.keys()])
    if np.all(max_aafs < 0.1):
        score += 2

        if np.all(max_aafs == 0):
            score += 1
        else:
            qcode += "m"

    else:
        qcode += "M"

    if pp.max_indel_size < 50:
        score += 1

        if pp.max_indel_size == 0:
            score += 1
        else:
            qcode += "i"

    else:
        qcode += "I"

    return score, qcode


def score_hrm(pp):
    score = 0
    qcode = ""

    if pp.hrm["delta"] >= 0.5:
        score += 2
        if pp.hrm["delta"] >= 1:
            score += 1
        else:
            qcode += "h"
    else:
        qcode += "H"

    tms = np.array([pp.primers[d].tm for d in pp.primers.keys()])
    tms_l1_norm = np.sum(np.abs(tms - np.mean(tms)))
    if tms_l1_norm <= 5:
        score += 1
    else:
        qcode += "t"

    if len(pp.offtargets) == 0:
        score += 1
    else:
        qcode += "O"

    if not pp.dimer:
        score += 1
    else:
        qcode += "d"

    max_aafs = np.array([pp.primers[d].max_aaf for d in pp.primers.keys()])
    if np.all(max_aafs < 0.1):
        score += 1

        if np.all(max_aafs == 0):
            score += 1
        else:
            qcode += "m"

    else:
        qcode += "M"

    if pp.max_indel_size < 50:
        score += 1

        if pp.max_indel_size == 0:
            score += 1
        else:
            qcode += "i"

    else:
        qcode += "I"

    return score, qcode


def score_kasp(pp):
    score = 0
    qcode = ""

    tms = np.array([pp.primers[d].tm for d in pp.primers.keys()])
    tms_l1_norm = np.sum(np.abs(tms - np.mean(tms)))
    if tms_l1_norm <= 5:
        score += 1
    else:
        qcode += "t"

    if len(pp.offtargets) == 0:
        score += 3
    else:
        qcode += "O"

    if (not pp.ref_dimer) and (not pp.alt_dimer):
        score += 1
    else:
        qcode += "d"

    max_aafs = np.array([pp.primers[d].max_aaf for d in pp.primers.keys()])
    if np.all(max_aafs < 0.1):
        score += 2

        if np.all(max_aafs == 0):
            score += 1
        else:
            qcode += "m"

    else:
        qcode += "M"

    if pp.max_indel_size < 50:
        score += 1

        if pp.max_indel_size == 0:
            score += 1
        else:
            qcode += "i"

    else:
        qcode += "I"

    return score, qcode


def get_score_fnc(assay_name):
    if assay_name == "KASP":
        return score_kasp
    if assay_name == "PCR":
        return score_base
    if assay_name == "CAPS":
        return score_base
    if assay_name == "HRM":
        return score_hrm
