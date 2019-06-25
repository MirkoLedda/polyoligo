import re
# noinspection PyPackageRequirements
from Bio.Seq import Seq
from os.path import join
from primer3.thermoanalysis import ThermoAnalysis

from . import lib_blast


# noinspection PyPep8Naming
class gRNA:
    def __init__(self, chrom, start, end, seq, strand, pam_len):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq[:-pam_len]
        self.seq_pam = seq
        self.strand = strand
        if self.strand == "plus":
            self.name = "{}_{}_{}".format(self.chrom, self.start, self.end + pam_len)
        else:
            self.name = "{}_{}_{}".format(self.chrom, self.start, self.end - pam_len)

        # self.offtargets = None
        # self.n_offtargets = None
        self.offtargets_8mer = None
        self.n_offtargets_8mer = None
        self.offtargets_12mer = None
        self.n_offtargets_12mer = None
        self.perfect_offtargets = None
        self.n_perfect_offtargets = None
        self.do_seed_analysis = False
        self.n_8mers = "NA"
        self.n_12mers = "NA"
        self.tm = None
        self.T_runs = None

    def calc_TM(self):
        thals_obj = ThermoAnalysis()  # Initialize a primer3 thal object
        self.tm = thals_obj.calcTm(self.seq)

    def get_run_Ts(self):
        if "TTTT" in self.seq:
            self.T_runs = 1
        else:
            self.T_runs = 0


# noinspection PyPep8Naming
class Crispr:
    def __init__(self, chrom, start, end, pam, blast_db):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.pam = pam
        self.pam_len = len(pam)
        self.pam_regex = pam
        self.pam_regex = re.compile(self.pam_regex.replace("N", "[ATGCN]"))
        self.blast_db = blast_db

        self.seq = None
        self.gRNAs = None

    def is_pam(self, seq):
        return self.pam_regex.search(seq)

    def fetch_roi(self):
        query = [{
            "chr": self.chrom,
            "start": self.start,
            "stop": self.end,
        }]

        seqs = self.blast_db.fetch(query)
        self.seq = list(seqs.values())[0]
        self.seq = self.seq.upper()

    def find_gRNAs(self):
        n_nuc = len(self.seq)
        grna_len = 20

        grnas = []
        # Find forward matches
        for match in re.finditer(self.pam_regex, self.seq):
            j = match.span()[0]
            if (j - grna_len) >= 0:
                grna = gRNA(
                    chrom=self.chrom,
                    start=self.start + (j - 1) - (grna_len - 1),
                    end=self.start + (j - 1),
                    seq=self.seq[(j - grna_len):(j + self.pam_len)],
                    strand="plus",
                    pam_len=self.pam_len,
                )
                grnas.append(grna)

        # Find reverse matches
        iseq = str(Seq(self.seq).reverse_complement())

        for match in re.finditer(self.pam_regex, iseq):
            j = match.span()[0]
            if (j - grna_len) >= 0:
                grna = gRNA(
                    chrom=self.chrom,
                    start=self.start + (n_nuc - j) + (grna_len - 1),
                    end=self.start + (n_nuc - j),
                    seq=iseq[(j - grna_len):(j + self.pam_len)],
                    strand="minus",
                    pam_len=self.pam_len,
                )
                grnas.append(grna)

        # Remove guides containing Ns
        grnas_parsed = []
        for grna in grnas:
            if "N" not in grna.seq:
                grnas_parsed.append(grna)
        self.gRNAs = grnas_parsed

    def check_offtargeting(self, batch_size=100, logger=None):

        n_total = len(self.gRNAs)
        batch_id = 0
        gRNAs_batches = {batch_id: []}
        for grna in self.gRNAs:
            gRNAs_batches[batch_id].append(grna)

            if len(gRNAs_batches[batch_id]) > batch_size:
                batch_id += 1
                gRNAs_batches[batch_id] = []

        # Run each batch iteratively
        self.gRNAs = []

        for k, gRNAs_batch in gRNAs_batches.items():
            gRNAs_batch = self.check_offtargeting_minibatch(gRNAs_batch=gRNAs_batch)
            self.gRNAs += gRNAs_batch
            if logger is not None:
                logger.info("{}/{}".format(min((k + 1) * batch_size, n_total), n_total))

    def check_offtargeting_minibatch(self, gRNAs_batch):
        # Build query dictionaries
        queries = {}

        for grna in gRNAs_batch:
            queries[grna.name] = grna.seq + self.pam

        # Forward guide lookup
        fp_query = join(self.blast_db.temporary, "{}_blast.fa".format(self.blast_db.job_id))
        fp_out = join(self.blast_db.temporary, "{}_blast.json".format(self.blast_db.job_id))
        lib_blast.write_fasta(queries, fp_query)

        self.blast_db.blastn(
            fp_query=fp_query,
            fp_out=fp_out,
            word_size=8,
            evalue=1000,
            # max_hsps=5,
            perc_identity=86,
            dust="no",
            soft_masking="false",
            outfmt='"15"',
        )

        hits = self.blast_db.parse_blastn_json(fp_out, min_identity=20)

        # Remove target site
        offtargets_8mer = {}
        offtargets_12mer = {}
        perfect_offtargets = {}
        for target_name, chrhits in hits.items():
            found_target = False
            offtargets_8mer[target_name] = []
            offtargets_12mer[target_name] = []
            perfect_offtargets[target_name] = []
            for chrom, shits in chrhits.items():
                for shit in shits.values():
                    if not target_name == "{}_{}_{}".format(chrom, shit["sstart"], shit["sstop"]):
                        offtarget = {
                            "chrom": chrom,
                            "start": shit["sstart"],
                            "end": shit["sstop"],
                        }

                        offt8, offt12 = self.score_offtargeting(shit["midline"])

                        # will_offtarget = self.score_offtargeting(shit["midline"])
                        # if will_offtarget:
                        #     offtargets[target_name].append(offtarget)
                        if offt8:
                            offtargets_8mer[target_name].append(offtarget)
                        if offt12:
                            offtargets_12mer[target_name].append(offtarget)

                        if offt8 or offt12:
                            if self.is_perfect_offtarget(shit["midline"]):
                                perfect_offtargets[target_name].append(offtarget)
                    else:
                        found_target = True
            if not found_target:
                print("WARNING - Target not found: {}".format(target_name))

        # Add offtargets to the gRNA objects
        for grna in gRNAs_batch:
            # grna.offtargets = offtargets[grna.name]
            # grna.n_offtargets = len(grna.offtargets)
            grna.offtargets_8mer = offtargets_8mer[grna.name]
            grna.n_offtargets_8mer = len(grna.offtargets_8mer)
            grna.offtargets_12mer = offtargets_12mer[grna.name]
            grna.n_offtargets_12mer = len(grna.offtargets_12mer)
            grna.perfect_offtargets = perfect_offtargets[grna.name]
            grna.n_perfect_offtargets = len(grna.perfect_offtargets)

            if grna.n_perfect_offtargets == 0:
                grna.do_seed_analysis = True

        return gRNAs_batch

    def check_seeds(self, batch_size=100, logger=None):

        n_total = len(self.gRNAs)
        batch_id = 0
        gRNAs_batches = {batch_id: []}
        for grna in self.gRNAs:
            gRNAs_batches[batch_id].append(grna)

            if len(gRNAs_batches[batch_id]) > batch_size:
                batch_id += 1
                gRNAs_batches[batch_id] = []

        # Run each batch iteratively
        self.gRNAs = []

        for k, gRNAs_batch in gRNAs_batches.items():
            gRNAs_batch = self.check_seeds_minibatch(gRNAs_batch=gRNAs_batch)
            self.gRNAs += gRNAs_batch
            if logger is not None:
                logger.info("{}/{}".format(min((k + 1) * batch_size, n_total), n_total))

    def check_seeds_minibatch(self, gRNAs_batch):

        # Check first it we need to do this
        do_seeds = False
        for grna in gRNAs_batch:
            if grna.do_seed_analysis:
                do_seeds = True
                break

        if do_seeds:
            # 8-MERS
            # Build query dictionaries
            queries = {}

            for grna in gRNAs_batch:
                if grna.do_seed_analysis:
                    queries[grna.name] = grna.seq[12:] + self.pam

            # Forward strand lookup
            fp_query = join(self.blast_db.temporary, "{}_blast.fa".format(self.blast_db.job_id))
            fp_out = join(self.blast_db.temporary, "{}_blast.json".format(self.blast_db.job_id))
            lib_blast.write_fasta(queries, fp_query)

            self.blast_db.blastn(
                fp_query=fp_query,
                fp_out=fp_out,
                word_size=8,
                evalue=1000000,
                max_hsps=1000,
                perc_identity=90,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
            )

            hits = self.blast_db.parse_blastn_json(fp_out, min_identity=10)

            # Remove target site
            n_8mers = {}
            for target_name, chrhits in hits.items():

                n_8mers[target_name] = -1  # One is going to be the target
                for chrom, shits in chrhits.items():
                    for _ in shits.values():
                        n_8mers[target_name] += 1

            # Add offtargets to the gRNA objects
            for grna in gRNAs_batch:
                if grna.do_seed_analysis:
                    grna.n_8mers = n_8mers[grna.name]

            # 12-MERS
            # Build query dictionaries
            queries = {}

            for grna in gRNAs_batch:
                if grna.do_seed_analysis:
                    queries[grna.name] = grna.seq[8:] + self.pam

            # Forward strand lookup
            fp_query = join(self.blast_db.temporary, "{}_blast.fa".format(self.blast_db.job_id))
            fp_out = join(self.blast_db.temporary, "{}_blast.json".format(self.blast_db.job_id))
            lib_blast.write_fasta(queries, fp_query)

            self.blast_db.blastn(
                fp_query=fp_query,
                fp_out=fp_out,
                word_size=12,
                evalue=1000000,
                max_hsps=1000,
                perc_identity=93,
                dust="no",
                soft_masking="false",
                outfmt='"15"',
            )

            hits = self.blast_db.parse_blastn_json(fp_out, min_identity=14)

            # Remove target site
            n_12mers = {}
            for target_name, chrhits in hits.items():

                n_12mers[target_name] = -1  # One is going to be the target
                for chrom, shits in chrhits.items():
                    for _ in shits.values():
                        n_12mers[target_name] += 1

            # Add offtargets to the gRNA objects
            for grna in gRNAs_batch:
                if grna.do_seed_analysis:
                    grna.n_12mers = n_12mers[grna.name]

        return gRNAs_batch

    @staticmethod
    def score_offtargeting(midline):
        n = len(midline)
        midline_pam = midline[(n - 3):n]
        midline = midline[:-3]

        mer8 = midline[(len(midline) - 8):len(midline)].count("X") == 0  # Check the 8mer seed region
        mer12 = midline[(len(midline) - 12):len(midline)].count("X") == 0  # Check the 12mer seed region
        pam = midline_pam == "X||"  # Check the PAM site

        if pam:
            fmer8 = mer8
            fmer12 = mer12
        else:
            fmer8 = False
            fmer12 = False

        return fmer8, fmer12

        # mer12 = midline[(len(midline) - 12):len(midline)].count("X") == 0   # Check the 12mer seed region
        # pam = midline_pam == "X||"  # Check the PAM site
        #
        # if mer12 and pam:
        #     return True
        # else:
        #     return False

    @staticmethod
    def is_perfect_offtarget(midline):
        if midline.count("X") == 1:
            return True
        else:
            return False

    def calc_TM(self):
        for grna in self.gRNAs:
            grna.calc_TM()

    def get_run_Ts(self):
        for grna in self.gRNAs:
            grna.get_run_Ts()

    @staticmethod
    def write_report_header(fp, fp_bed):
        with open(fp, "w") as f:
            f.write(
                "CHR START END STRAND SEQ PAM N_PERFECT_OFFTARGETS N_OFFTARGETS_12MER N_OFFTARGETS_8MER 12_MER 8_MER TM TTTT\n")

        with open(fp_bed, "w") as f:
            f.write("")

    def write_report(self, fp, fp_bed):
        with open(fp, "a") as f:
            for grna in self.gRNAs:
                if grna.strand == "plus":
                    str_strand = "+"
                else:
                    str_strand = "-"

                f.write("{} {} {} {} {} {} {} {} {} {} {} {:.1f} {}\n".format(
                    self.chrom,
                    grna.start,
                    grna.end,
                    str_strand,
                    grna.seq,
                    grna.seq_pam[20:],
                    grna.n_perfect_offtargets,
                    grna.n_offtargets_12mer,
                    grna.n_offtargets_8mer,
                    grna.n_12mers,
                    grna.n_8mers,
                    grna.tm,
                    grna.T_runs,
                ))

        with open(fp_bed, "a") as f:
            for grna in self.gRNAs:
                if grna.strand == "plus":
                    str_strand = "+"
                else:
                    str_strand = "-"

                str_name = "{}_{}_{}_{}_{}_{:.0f}_{}".format(
                    grna.n_perfect_offtargets,
                    grna.n_offtargets_12mer,
                    grna.n_offtargets_8mer,
                    grna.n_12mers,
                    grna.n_8mers,
                    grna.tm,
                    grna.T_runs,
                )

                f.write("{} {} {} {} {} {}\n".format(
                    self.chrom,
                    grna.start,
                    grna.end,
                    str_name,
                    0,
                    str_strand,
                ))


if __name__ == "__main__":
    pass
