from __future__ import print_function, division
from subprocess import call
import pandas as pd
import json
import numpy as np
import time
from os.path import join, exists
import gzip
# noinspection PyPackageRequirements
from Bio.Seq import Seq
from uuid import uuid4


class BlastDB:
    def __init__(self, path_db, path_temporary, path_bin, job_id="", n_cpu=1):
        self.db = path_db
        self.has_db = True
        self.has_fasta = False
        self.fasta = None
        self.temporary = path_temporary
        self.exe = path_bin
        self.job_id = job_id
        self.n_cpu = n_cpu
        self.seqs = None

        if exists(self.db) and self.db.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
            self.fasta = self.db
            self.has_fasta = True
            self.has_db = False

    def locked_call(self, cmd):
        task_id = str(time.time()) + str(uuid4())
        lock = join(self.temporary, "{}.lock".format(task_id))
        lock_cmd = "touch {}".format(lock)
        unlock_cmd = "rm {}".format(lock)

        lcmd = "{};{};{};".format(lock_cmd, cmd, unlock_cmd)
        call(lcmd, shell=True)

        while exists(lock):
            time.sleep(0.1)

    def load_fasta(self):
        if self.has_fasta:
            self.seqs = read_fasta(self.fasta)

    def fasta2db(self):
        self.db = join(self.temporary, "blastdb")

        if self.fasta.endswith(".gz"):
            cmd = [
                "gunzip",
                "-c",
                self.fasta,
                "|",
                join(self.exe, "makeblastdb"),
                "-in -",
                "-title blastdb",
                "-input_type fasta",
                "-dbtype nucl",
                "-out {}".format(self.db),
                "-parse_seqids",
            ]
        else:
            cmd = [
                join(self.exe, "makeblastdb"),
                "-in {}".format(self.fasta),
                "-input_type fasta",
                "-dbtype nucl",
                "-out {}".format(self.db),
                "-parse_seqids",
            ]

        cmd.append("> {}".format(join(self.temporary, "devnull")))
        cmd = " ".join(cmd)

        self.locked_call(cmd)
        self.has_db = True

    def db2fasta(self):
        if not self.has_fasta:
            self.fasta = join(self.temporary, "refG.fa")
            cmd = [
                join(self.exe, "blastdbcmd"),
                "-entry all",
                "-db {}".format(self.db),
                "-out {}".format(self.fasta),
            ]

            cmd.append("> {}".format(join(self.temporary, "devnull")))
            cmd = " ".join(cmd)

            self.locked_call(cmd)
            self.has_fasta = True
            self.load_fasta()

    def purge(self):
        self.seqs = None

    def fetch(self, *args, **kwargs):
        if self.has_fasta:
            return self.fetch_fasta(*args, **kwargs)
        else:
            return self.fetch_blastdb(*args, **kwargs)

    def fetch_blastdb(self, queries, fp_out=None, **kwargs):
        """Run a blastdbcmd command based on a dictionary of queries. Accepts additional blastdbcmd arguments.

        Args:
            queries (list): List of dictionary with keys=["chr", "start", "stop", "strand" (optional)].
            fp_out (str): Output filename.

        Returns:
            seqs (dict): Dictionary of sequences keyed by fasta id.

        """

        if not self.has_db:
            self.fasta2db()

        fp_blast_in = join(self.temporary, "{}_blastdbcmd_in.txt".format(self.job_id))

        if fp_out is None:
            fp_blast_out = join(self.temporary, "{}_blastdbcmd_out.fa".format(self.job_id))
        else:
            fp_blast_out = fp_out

        with open(fp_blast_in, "w") as f:
            for q in queries:
                if "strand" in q.keys():
                    if q["strand"] == "plus":
                        f.write("{}\t{:d}-{:d}\t{}\n".format(q["chr"], q["start"], q["stop"], q["strand"]))
                    else:
                        f.write("{}\t{:d}-{:d}\t{}\n".format(q["chr"], q["start"], q["stop"], q["strand"]))
                else:
                    f.write("{}\t{:d}-{:d}\n".format(q["chr"], q["start"], q["stop"]))

        blast_kwargs = []
        for k, v in kwargs.items():
            blast_kwargs.append("-{} {}".format(k, v))

        cmd = [
            join(self.exe, "blastdbcmd"),
            "-entry_batch {}".format(fp_blast_in),
            "-db {}".format(self.db),
            "-dbtype nucl",
        ]

        cmd += blast_kwargs  # Add blast kwargs
        cmd += ["> {}".format(fp_blast_out)]  # Add output pipe
        cmd = " ".join(cmd)

        self.locked_call(cmd)

        return read_fasta(fp_blast_out)

    def fetch_fasta(self, queries, fp_out=None, **kwargs):
        """Fetch sequences based on a dictionary of queries.

        Args:
            queries (list): List of dictionary with keys=["chr", "start", "stop", "strand" (optional)].
            fp_out (str): Output filename.

        Returns:
            seqs (dict): Dictionary of sequences keyed by fasta id.

        """

        fetched_seqs = {}

        for q in queries:
            start0 = q["start"] - 1
            stop0 = q["stop"]
            qkey = "{}:{}-{}".format(q["chr"], q["start"], q["stop"])

            if "strand" in q.keys():
                if q["strand"] == "plus":
                    fetched_seqs[qkey] = self.seqs[q["chr"]][start0:stop0]
                else:
                    qkey = "{}:c{}-{}".format(q["chr"], q["stop"], q["start"])
                    fetched_seqs[qkey] = self.seqs[q["chr"]][start0:stop0]
                    fetched_seqs[qkey] = str(Seq(fetched_seqs[qkey]).reverse_complement())
            else:
                fetched_seqs[qkey] = self.seqs[q["chr"]][start0:stop0]

        if fp_out is not None:
            write_fasta(fetched_seqs, fp_out)

        return fetched_seqs

    def blastn(self, fp_query, fp_out=None, **kwargs):
        """Run a blastn command with specified blastn arguments.."""

        if not self.has_db:
            self.fasta2db()

        # Parse blast kwargs
        blast_kwargs = []
        for k, v in kwargs.items():
            blast_kwargs.append("-{} {}".format(k, v))

        if fp_out is None:
            fp_out = join(self.temporary, "{}_blastn_out.txt".format(self.job_id))

        cmd = [
            join(self.exe, "blastn"),
            "-task blastn",
            "-db {}".format(self.db),
            "-query {}".format(fp_query),
            "-num_threads {:d}".format(self.n_cpu),
            "-out {}".format(fp_out),
        ]

        cmd += blast_kwargs  # Add blast kwargs
        cmd = " ".join(cmd)

        self.locked_call(cmd)

    @staticmethod
    def prune_blastn_output(fp, fp_out, min_align_len=50, max_n_seqs=np.inf):
        cnt = 0
        with open(fp, "r") as f, open(fp_out, "w") as fout:
            for line in f:
                fields = line.strip().split("\t")
                alen = int(fields[3])

                if alen >= min_align_len:
                    fout.write(line)
                    cnt += 1

                if cnt > max_n_seqs:
                    break

    @staticmethod
    def parse_blastn_output(fp):

        df = pd.read_csv(fp, delim_whitespace=True, header=None)
        df.columns = ["qname", "schr", "pct_id", "alen", "amis", "agap", "qstart", "qstop", "sstart", "sstop",
                      "evalue", "bitscore", "qseq", "sseq", "slen"]

        return df

    @staticmethod
    def parse_blastn_json(fp, min_identity=-np.inf):
        while True:
            try:
                json_dict = read_json(fp)
                break
            except json.decoder.JSONDecodeError:
                pass

        hits = {}
        for qhit in json_dict["BlastOutput2"]:
            query_title = qhit["report"]["results"]["search"]["query_title"]
            query_len = qhit["report"]["results"]["search"]["query_len"]
            hits[query_title] = {}

            for hit in qhit["report"]["results"]["search"]["hits"]:
                chrom = hit["description"][0]["id"]

                for i, hsp in enumerate(hit["hsps"]):
                    qstart = hsp["query_from"]
                    qstop = hsp["query_to"]
                    midline = hsp["midline"].replace(" ", "X")
                    alen = hsp["identity"]

                    if alen >= min_identity:
                        # Pad the midline with mismatches on both sides
                        ngaps = hsp["qseq"].count("-")
                        midline = (qstart - 1) * "X" + midline + (query_len + ngaps - qstop) * "X"

                        try:
                            hits[query_title][chrom][i] = {"midline": midline,
                                                           "sstart": hsp["hit_from"],
                                                           "sstop": hsp["hit_to"]}
                        except KeyError:
                            hits[query_title][chrom] = {}
                            hits[query_title][chrom][i] = {"midline": midline,
                                                           "sstart": hsp["hit_from"],
                                                           "sstop": hsp["hit_to"]}

        return hits


class Muscle:
    def __init__(self, path_temporary, exe):
        self.exe = exe
        self.temporary = path_temporary

    def locked_call(self, cmd):
        lock = join(self.temporary, "{}.lock".format(time.time()))
        lock_cmd = "touch {}".format(lock)
        unlock_cmd = "rm {}".format(lock)

        lcmd = "{};{};{};".format(lock_cmd, cmd, unlock_cmd)
        call(lcmd, shell=True)

        while exists(lock):
            time.sleep(0.1)

    def align_fasta(self, fp, fp_out=None):
        if fp_out is None:
            fp_out = join(self.temporary, "muscle_alignment.fa")

        cmd = " ".join([
            join(self.exe, "muscle"),
            "-in {}".format(fp),
            "-out {}".format(fp_out),
            "-quiet",
        ])
        self.locked_call(cmd)


def write_fasta(seqs, fp_out):
    with open(fp_out, "w") as f:
        if isinstance(seqs, list) or isinstance(seqs, np.ndarray):
            for i, seq in enumerate(seqs):
                f.write(">{}\n{}\n".format(i, seq))
        elif isinstance(seqs, dict):
            for k, seq in seqs.items():
                f.write(">{}\n{}\n".format(k, seq))


def read_fasta(fp, filters=None, return_ordering=False):

    if fp.endswith(".gz"):
        opener = gzip.open(fp, "rb")
    else:
        opener = open(fp, "r")

    seqs = {}
    key_ordered = []

    with opener as f:

        name = None
        seq = ""

        for line in f:
            if isinstance(line, bytes):
                line = line.decode("utf-8")

            if line.startswith(">"):
                curr_name = line.split(">")[1].strip()

                if name is not None:
                    if filters is not None:
                        if name in filters:
                            seqs[name] = seq
                            key_ordered.append(name)
                    else:
                        seqs[name] = seq
                        key_ordered.append(name)

                name = curr_name
                seq = ""
            else:
                seq += line.strip()

        seqs[name] = seq  # Add the last entry
        key_ordered.append(name)

    if return_ordering:
        return seqs, key_ordered
    else:
        return seqs


def read_json(fp):
    with open(fp, "r") as f:
        json_dict = json.load(f)
    return json_dict


if __name__ == "__main__":
    pass