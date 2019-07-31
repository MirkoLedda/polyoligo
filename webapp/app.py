import os
from flask import Flask, request, redirect, render_template, url_for, jsonify, send_file
from werkzeug.utils import secure_filename
from os.path import join
from datetime import datetime, timedelta
from uuid import uuid4
import json
import time
import subprocess
from subprocess import STDOUT
from celery import Celery
import yaml

# APP CONFIGS ------------------------------------------------------------------
app = Flask(__name__, instance_relative_config=True)
app.config.from_pyfile('config.py')
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

# GLOBALS ----------------------------------------------------------------------
DEPTH_MAP = {"superficial": 2.5,
             "standard": 10,
             "deep": 25,
             }


# SERVER TO CLIENT MESSAGING ---------------------------------------------------
class Log():
    def __init__(self, dest):
        self.dest = dest
        self.fp = join(dest, "logfile.json")
        self.content = {
            "status": "PENDING",
            "messages": [],
            "n": None,
            "n_total": None,
            "start_time": None,
        }
        self.timestamp()

    def timestamp(self):
        self.content["start_time"] = time.time()

    def update_timer(self):
        self.content["elapsed"] = seconds_to_hms(timedelta(seconds=time.time() - self.content["start_time"]))

    def update_nanobar(self, msg):
        msg = msg.split("-")[2].strip()
        nanobar = msg.split("/")
        self.content["n"] = int(nanobar[0])
        self.content["n_total"] = int(nanobar[1])
        self.write()

    def write(self):
        self.update_timer()
        with open(self.fp, "w") as f:
            json.dump(self.content, f)

    def fetch(self):
        with open(self.fp, "r") as f:
            self.content = json.load(f)

    def get(self):
        self.fetch()
        return self.content

    def message(self, msg):
        self.content["messages"].append(msg)
        self.write()

    def flush(self):
        self.content = {
            "status": "PENDING",
            "messages": [],
            "n": None,
            "n_total": None,
            "start_time": None,
        }
        self.write()

    def update_status(self, status):
        self.content["status"] = status
        self.write()


# FUNCTIONS --------------------------------------------------------------------
def create_task_dest():
    task_id = datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())
    dest_folder = join(app.config['UPLOAD_FOLDER'], task_id)
    return task_id, dest_folder


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


def upload_file_to_server(req_obj, k, dest):
    fp = None

    file = req_obj.files[k]
    filename = secure_filename(file.filename)
    if filename != "":
        fp = join(dest, filename)
        file.save(fp)

    return fp


def parse_range(txt):
    range = txt.strip().split("-")

    if len(range) == 1:
        return int(range[0])
    else:
        return [int(range[0]), int(range[1])]


def parse_blastdb_kwargs(kwargs, request, dest_folder):
    kwargs["include_vcf"] = False

    # BlastDB
    reference = request.form.get("reference")
    if reference == 'Fragaria x ananassa (KnappLab)':
        kwargs["include_vcf"] = True
    kwargs["reference"] = join(app.config["BLASTDB_REPO"], app.config["BLASTDB_FILES"][reference])

    # Populations selection
    if kwargs["include_vcf"]:
        selected_pops = []
        for i in range(9):
            current_pop = "population{}".format(i)
            fetched_pop = request.form.get(current_pop)

            if fetched_pop:
                with open(os.path.abspath(join(app.config["VCF_SUBSETS"], fetched_pop)), "r") as f:
                    for line in f:
                        selected_pops.append(line.strip())
        if len(selected_pops) == 0:
            kwargs["include_vcf"] = False
        else:
            with open(join(dest_folder, "vcf_include.txt"), "w") as f:
                for selected_pop in selected_pops:
                    f.write("{}\n".format(selected_pop))

    strcmd = []
    if kwargs["include_vcf"]:
        strcmd += ["--vcf {}".format(join(app.config["VCF_REPO"], app.config["BLASTDB_FILES"][reference] + ".gz"))]
        strcmd += ["--vcf_include {}".format(join(dest_folder, "vcf_include.txt"))]

    return kwargs, strcmd


# CELERY -----------------------------------------------------------------------
@celery.task()
def run_task(strcmd, log_dest):
    log = Log(log_dest)
    log.fetch()
    log.update_status("RUNNING")

    try:
        process = subprocess.Popen(strcmd, shell=True, stderr=STDOUT, stdout=subprocess.PIPE)
        while True:
            output = process.stdout.readline()
            if process.poll() is not None:
                break
            if output:
                msg = output.decode("utf-8")
                if msg.startswith("INFO - nanobar - "):
                    log.update_nanobar(msg)
                else:
                    log.message(msg)
    except subprocess.CalledProcessError:
        log.update_status("FAILED")

    if not os.path.exists(join(log_dest, "output.txt")):
        log.update_status("FAILED")
    else:
        log.update_status("COMPLETED")

        # Compress result files for download
        strcmd = [
            "cd {};".format(log_dest),
            "ls;",
            "tar -cvf 'output.tar' output.*;",
            "gzip output.tar;",
            "cd ..;",
        ]
        strcmd = " ".join(strcmd)
        subprocess.Popen(strcmd, shell=True)

    return None


# ROUTES -----------------------------------------------------------------------
@app.route('/')
def home():
    return render_template('home.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


@app.route('/pcr', methods=['GET', 'POST'])
def pcr():
    exe = "polyoligo-pcr"
    if request.method == 'POST':

        kwargs_names = ["roi", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth", "primer3"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # ROI
        kwargs["roi"] = join(dest_folder, "rois.txt")
        rois = request.form.get("rois")
        with open(kwargs["roi"], "w") as f:
            f.write(rois)

        # Genome kwargs
        kwargs, nstrcmd = parse_blastdb_kwargs(kwargs=kwargs, request=request, dest_folder=dest_folder)

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))
        kwargs["depth"] = request.form.get("depth")
        kwargs["depth"] = DEPTH_MAP[kwargs["depth"]]
        kwargs["seed"] = int(request.form.get("seed"))
        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        # Primer3 configs
        primer3_yaml = join(dest_folder, "primer3.yaml")
        primer3 = {}
        p_tm = request.form.get("p_tm")
        p_size = request.form.get("p_size")
        p_gc = request.form.get("p_gc")
        p_seq_F = request.form.get("p_seq_F")
        p_seq_R = request.form.get("p_seq_R")

        if len(p_tm) > 0:
            r = parse_range(p_tm)
            if isinstance(r, list):
                primer3["PRIMER_OPT_TM"] = int(r[0] + (r[1] - r[0]) / 2)
                primer3["PRIMER_MIN_TM"] = r[0]
                primer3["PRIMER_MAX_TM"] = r[1]
            else:
                primer3["PRIMER_OPT_TM"] = r
                primer3["PRIMER_MIN_TM"] = r
                primer3["PRIMER_MAX_TM"] = r
        if len(p_size) > 0:
            r = parse_range(p_size)
            if isinstance(r, list):
                primer3["PRIMER_OPT_SIZE"] = int(r[0] + (r[1] - r[0]) / 2)
                primer3["PRIMER_MIN_SIZE"] = r[0]
                primer3["PRIMER_MAX_SIZE"] = r[1]
            else:
                primer3["PRIMER_OPT_SIZE"] = r
                primer3["PRIMER_MIN_SIZE"] = r
                primer3["PRIMER_MAX_SIZE"] = r
        if len(p_gc) > 0:
            r = parse_range(p_gc)
            if isinstance(r, list):
                primer3["PRIMER_MIN_GC"] = r[0]
                primer3["PRIMER_MAX_GC"] = r[1]
            else:
                primer3["PRIMER_MIN_GC"] = r
                primer3["PRIMER_MAX_GC"] = r
        if len(p_seq_F) > 0:
            primer3["SEQUENCE_PRIMER"] = p_seq_F
        if len(p_seq_R) > 0:
            primer3["SEQUENCE_PRIMER_REVCOMP"] = p_seq_R

        with open(primer3_yaml, "w") as f:
            yaml.dump(primer3, f)

        strcmd = [
            exe,
            kwargs["roi"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--primer3 {}".format(primer3_yaml),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        strcmd += nstrcmd
        strcmd = " ".join(strcmd)
        log = Log(dest_folder)
        log.message("File(s) successfully uploaded ...")
        log.message("Job started ...")
        log.write()

        run_task.delay(strcmd, dest_folder)

        return redirect(url_for('processing', task_id=task_id))
    else:
        return render_template('pcr.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


@app.route('/kasp', methods=['GET', 'POST'])
def kasp():
    exe = "polyoligo-kasp"
    if request.method == 'POST':

        kwargs_names = ["markers", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth", "dye1", "dye2"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # Marker file
        kwargs["markers"] = join(dest_folder, "markers.txt")
        markers = request.form.get("markers")
        with open(kwargs["markers"], "w") as f:
            f.write(markers)

        # Genome kwargs
        kwargs, nstrcmd = parse_blastdb_kwargs(kwargs=kwargs, request=request, dest_folder=dest_folder)

        # Dyes
        kwargs["dye1"] = request.form.get("dye1")
        kwargs["dye2"] = request.form.get("dye2")

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))
        kwargs["depth"] = request.form.get("depth")
        kwargs["depth"] = DEPTH_MAP[kwargs["depth"]]
        kwargs["seed"] = int(request.form.get("seed"))

        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        strcmd = [
            exe,
            kwargs["markers"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--dye1 {}".format(kwargs["dye1"]),
            "--dye2 {}".format(kwargs["dye2"]),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        strcmd += nstrcmd
        strcmd = " ".join(strcmd)
        log = Log(dest_folder)
        log.message("File(s) successfully uploaded ...")
        log.message("Job started ...")
        log.write()

        run_task.delay(strcmd, dest_folder)

        return redirect(url_for('processing', task_id=task_id))
    else:
        return render_template('kasp.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


@app.route('/caps', methods=['GET', 'POST'])
def caps():
    exe = "polyoligo-caps"
    if request.method == 'POST':
        kwargs_names = ["markers", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth", "enzymes",
                        "fragment_min_size"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # Marker file
        kwargs["markers"] = join(dest_folder, "markers.txt")
        markers = request.form.get("markers")
        with open(kwargs["markers"], "w") as f:
            f.write(markers)

        # Genome kwargs
        kwargs, nstrcmd = parse_blastdb_kwargs(kwargs=kwargs, request=request, dest_folder=dest_folder)

        # Restriction fragment
        kwargs["fragment_min_size"] = int(request.form.get("fragment_min_size"))

        # Enzymes
        enzymes = request.form.get("enzymes")
        if enzymes:
            kwargs["enzymes"] = join(dest_folder, "enzymes.txt")
            with open(kwargs["enzymes"], "w") as f:
                for enzyme in enzymes.strip().split():
                    f.write("{}\n".format(enzyme))

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))
        kwargs["depth"] = request.form.get("depth")
        kwargs["depth"] = DEPTH_MAP[kwargs["depth"]]
        kwargs["seed"] = int(request.form.get("seed"))

        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        strcmd = [
            exe,
            kwargs["markers"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--fragment_min_size {}".format(kwargs["fragment_min_size"]),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        if kwargs["enzymes"]:
            strcmd += ["--enzymes {}".format(kwargs["enzymes"])]

        strcmd += nstrcmd
        strcmd = " ".join(strcmd)
        log = Log(dest_folder)
        log.message("File(s) successfully uploaded ...")
        log.message("Job started ...")
        log.write()

        run_task.delay(strcmd, dest_folder)

        return redirect(url_for('processing', task_id=task_id))
    else:
        return render_template('caps.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


@app.route('/crispr', methods=['GET', 'POST'])
def crispr():
    exe = "polyoligo-crispr"
    if request.method == 'POST':

        kwargs_names = ["roi", "reference"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # ROI
        kwargs["roi"] = request.form.get("roi")
        if kwargs["roi"] == "":
            kwargs["roi"] = "CHR:0-0"

        # Genome kwargs
        kwargs, nstrcmd = parse_blastdb_kwargs(kwargs=kwargs, request=request, dest_folder=dest_folder)

        strcmd = [
            exe,
            kwargs["roi"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-nt 1",
            "--webapp",
        ]

        strcmd += nstrcmd
        strcmd = " ".join(strcmd)
        log = Log(dest_folder)
        log.message("File(s) successfully uploaded ...")
        log.message("Job started ...")
        log.write()

        run_task.delay(strcmd, dest_folder)

        return redirect(url_for('processing', task_id=task_id))
    else:
        return render_template('crispr.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


# PROCESSING ROUTES ------------------------------------------------------------
@app.route('/tasks/<task_id>')
def processing(task_id):
    return render_template("processing.html", task_id=task_id)


@app.route('/status/<task_id>', methods=['GET', 'POST'])
def get_status(task_id):
    log = Log(join(app.config['UPLOAD_FOLDER'], task_id))
    log.fetch()
    return jsonify(log.content)


@app.route('/downloads/<task_id>/output.txt', methods=['GET', 'POST'])
def download(task_id):
    filename = join(app.config['UPLOAD_FOLDER'], task_id, "output.tar.gz")
    return send_file(filename, as_attachment=True, attachment_filename="output.tar.gz")


# @app.route('/downloads/<task_id>/output.txt', methods=['GET', 'POST'])
# def download(task_id):
#     with open(join(app.config['UPLOAD_FOLDER'], task_id, "output.txt"), "r") as f:
#         content = f.read()
#     return render_template("results.html", content=content)


# OTHER ROUTES ------------------------------------------------------------
@app.route('/contig_list', methods=['GET', 'POST'])
def contig_list():
    with open(join(app.config['CONTIG_LIST']), "r") as f:
        content = f.read()
    return render_template("results.html", content=content)


if __name__ == "__main__":
    app.run()
