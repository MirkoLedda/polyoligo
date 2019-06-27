import os
from flask import Flask, flash, request, redirect, render_template, url_for, Response, jsonify, send_file
from werkzeug.utils import secure_filename
from os.path import join
from datetime import datetime, timedelta
from uuid import uuid4
from threading import Thread, Event
import json
import time
import subprocess
from subprocess import STDOUT, CalledProcessError
from celery import Celery
import yaml

# APP CONFIGS ------------------------------------------------------------------
app = Flask(__name__)
app.secret_key = "secret key"
app.config['UPLOAD_FOLDER'] = os.path.abspath('./uploads')
app.config["BLASTDB_REPO"] = os.path.abspath("../sample_data")
app.config["VCF_REPO"] = os.path.abspath("../sample_data")
app.config["BLASTDB_OPTIONS"] = ['Fragaria ananassa']
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024

app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)


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

# CELERY -----------------------------------------------------------------------
@celery.task()
def run_task(strcmd, log_dest):
    log = Log(log_dest)
    log.fetch()
    log.update_status("RUNNING")

    # log.message(strcmd)

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

    return None

# ROUTES -----------------------------------------------------------------------
@app.route('/')
def home():
    return render_template('home.html')


@app.route('/pcr', methods=['GET', 'POST'])
def pcr():
    exe = "polyoligo-pcr"
    if request.method == 'POST':
        include_vcf = False
        kwargs_names = ["roi", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth",
                        "tm_delta", "primer3"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # ROI
        kwargs["roi"] = request.form.get("roi")

        # BlastDB
        kwargs["reference"] = request.form.get("reference")

        if kwargs["reference"] == 'Fragaria ananassa':
            include_vcf = True

        # Populations selection
        if include_vcf:
            selected_pops = []
            for i in range(9):
                current_pop = "population{}".format(i)
                fetched_pop = request.form.get(current_pop)

                if fetched_pop:
                    with open(os.path.abspath(join("./static/data", fetched_pop)), "r") as f:
                        for line in f:
                            selected_pops.append(line.strip())
            if len(selected_pops) == 0:
                include_vcf = False
            else:
                with open(join(dest_folder, "vcf_include.txt"), "w") as f:
                    for selected_pop in selected_pops:
                        f.write("{}\n".format(selected_pop))

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))

        kwargs["depth"] = request.form.get("depth")
        depth_map = {"superficial": 2.5,
                     "standard": 10,
                     "deep": 25,
                     }
        kwargs["depth"] = depth_map[kwargs["depth"]]

        kwargs["tm_delta"] = float(request.form.get("tm_delta"))
        kwargs["seed"] = int(request.form.get("seed"))

        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        # Primer3 configs
        primer3_yaml = join(dest_folder, "primer3.yaml")
        primer3 = {}
        p_prod_len = request.form.get("p_prod_len")
        p_tm = request.form.get("p_tm")
        p_size = request.form.get("p_size")
        p_gc = request.form.get("p_gc")
        p_include = request.form.get("p_include")
        p_seq_F = request.form.get("p_seq_F")
        p_seq_R = request.form.get("p_seq_R")

        if len(p_prod_len) > 0:
            r = parse_range(p_prod_len)
            if isinstance(r, list):
                primer3["PRIMER_PRODUCT_SIZE_RANGE"] = r
            else:
                primer3["PRIMER_PRODUCT_SIZE"] = r
        if len(p_tm) > 0:
            r = parse_range(p_tm)
            if isinstance(r, list):
                primer3["PRIMER_OPT_TM"] = int(r[0] + (r[1]-r[0])/2)
                primer3["PRIMER_MIN_TM"] = r[0]
                primer3["PRIMER_MAX_TM"] = r[1]
            else:
                primer3["PRIMER_OPT_TM"] = r
                primer3["PRIMER_MIN_TM"] = r
                primer3["PRIMER_MAX_TM"] = r
        if len(p_size) > 0:
            r = parse_range(p_size)
            if isinstance(r, list):
                primer3["PRIMER_OPT_SIZE"] = int(r[0] + (r[1]-r[0])/2)
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
        if len(p_include) > 0:
            r = parse_range(p_include)
            start = int(kwargs["roi"].strip().split(":")[1].split("-")[0])
            if isinstance(r, list):
                primer3["SEQUENCE_TARGET"] = [r[0]-start, r[1]-r[0]]
            else:
                primer3["SEQUENCE_TARGET"] = [r-start, 1]
        if len(p_seq_F) > 0:
            primer3["SEQUENCE_PRIMER"] = p_seq_F
        if len(p_seq_R) > 0:
            primer3["SEQUENCE_PRIMER_REVCOMP"] = p_seq_R

        with open(primer3_yaml, "w") as f:
            yaml.dump(primer3, f)

        # Marker file
        primer3_file = upload_file_to_server(request, "primer3_file", dest_folder)
        if primer3_file is not None:
            primer3_yaml = primer3_file

        # todo rename depth
        kwargs["depth"] = kwargs["depth"]
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            exe,
            kwargs["roi"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--tm_delta {}".format(kwargs["tm_delta"]),
            "--primer3 {}".format(primer3_yaml),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        if include_vcf:
            # TODO: upload proper file
            strcmd += ["--vcf {}".format(join(app.config["VCF_REPO"], "vcf.txt.gz"))]
            strcmd += ["--vcf_include {}".format(join(dest_folder, "vcf_include.txt"))]

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
        include_vcf = False
        kwargs_names = ["markers", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth", "dye1", "dye2",
                        "tm_delta"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # Marker file
        kwargs["markers"] = upload_file_to_server(request, "markers_file", dest_folder)
        if kwargs["markers"] is None:
            kwargs["markers"] = join(dest_folder, "markers.txt")
            markers = request.form.get("markers")
            with open(kwargs["markers"], "w") as f:
                for line in markers:
                    f.write(line)

        # BlastDB
        kwargs["reference"] = request.form.get("ncbi_taxid")
        if kwargs["reference"] == "":
            kwargs["reference"] = request.form.get("reference")

        if kwargs["reference"] == 'Fragaria ananassa':
            include_vcf = True

        # Dyes
        kwargs["dye1"] = request.form.get("dye1")
        kwargs["dye2"] = request.form.get("dye2")

        # Populations selection
        if include_vcf:
            selected_pops = []
            for i in range(9):
                current_pop = "population{}".format(i)
                fetched_pop = request.form.get(current_pop)

                if fetched_pop:
                    with open(os.path.abspath(join("./static/data", fetched_pop)), "r") as f:
                        for line in f:
                            selected_pops.append(line.strip())
            if len(selected_pops) == 0:
                include_vcf = False
            else:
                with open(join(dest_folder, "vcf_include.txt"), "w") as f:
                    for selected_pop in selected_pops:
                        f.write("{}\n".format(selected_pop))

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))

        kwargs["depth"] = request.form.get("depth")
        depth_map = {"superficial": 2.5,
                     "standard": 10,
                     "deep": 25,
                     }
        kwargs["depth"] = depth_map[kwargs["depth"]]

        kwargs["tm_delta"] = float(request.form.get("tm_delta"))
        kwargs["seed"] = int(request.form.get("seed"))

        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        # todo rename depth
        kwargs["depth"] = kwargs["depth"]
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            exe,
            kwargs["markers"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--dye1 {}".format(kwargs["dye1"]),
            "--dye2 {}".format(kwargs["dye2"]),
            "--tm_delta {}".format(kwargs["tm_delta"]),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        if include_vcf:
            # TODO: upload proper file
            strcmd += ["--vcf {}".format(join(app.config["VCF_REPO"], "vcf.txt.gz"))]
            strcmd += ["--vcf_include {}".format(join(dest_folder, "vcf_include.txt"))]

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
        include_vcf = False
        kwargs_names = ["markers", "reference", "vcf", "vcf_include", "vcf_exclude", "n", "depth", "enzymes", "fragment_min_size",
                        "tm_delta"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # Marker file
        kwargs["markers"] = upload_file_to_server(request, "markers_file", dest_folder)
        if kwargs["markers"] is None:
            kwargs["markers"] = join(dest_folder, "markers.txt")
            markers = request.form.get("markers")
            with open(kwargs["markers"], "w") as f:
                for line in markers:
                    f.write(line)

        # BlastDB
        kwargs["reference"] = request.form.get("ncbi_taxid")
        if kwargs["reference"] == "":
            kwargs["reference"] = request.form.get("reference")

        if kwargs["reference"] == 'Fragaria ananassa':
            include_vcf = True

        # Restriction fragment
        kwargs["fragment_min_size"] = int(request.form.get("fragment_min_size"))

        # Enzymes
        enzymes = request.form.get("enzymes")
        if enzymes:
            kwargs["enzymes"] = join(dest_folder, "enzymes.txt")
            with open(kwargs["enzymes"], "w") as f:
                for enzyme in enzymes.strip().split():
                    f.write("{}\n".format(enzyme))

        # Populations selection
        if include_vcf:
            selected_pops = []
            for i in range(9):
                current_pop = "population{}".format(i)
                fetched_pop = request.form.get(current_pop)

                if fetched_pop:
                    with open(os.path.abspath(join("./static/data", fetched_pop)), "r") as f:
                        for line in f:
                            selected_pops.append(line.strip())
            if len(selected_pops) == 0:
                include_vcf = False
            else:
                with open(join(dest_folder, "vcf_include.txt"), "w") as f:
                    for selected_pop in selected_pops:
                        f.write("{}\n".format(selected_pop))

        # Other kwargs
        kwargs["n"] = int(request.form.get("n"))

        kwargs["depth"] = request.form.get("depth")
        depth_map = {"superficial": 2.5,
                     "standard": 10,
                     "deep": 25,
                     }
        kwargs["depth"] = depth_map[kwargs["depth"]]

        kwargs["tm_delta"] = float(request.form.get("tm_delta"))
        kwargs["seed"] = int(request.form.get("seed"))

        offtarget_size = parse_range(request.form.get("offtarget_size"))
        kwargs["offtarget_min_size"] = offtarget_size[0]
        kwargs["offtarget_max_size"] = offtarget_size[1]

        # todo rename depth
        kwargs["depth"] = kwargs["depth"]
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            exe,
            kwargs["markers"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--depth {}".format(kwargs["depth"]),
            "--fragment_min_size {}".format(kwargs["fragment_min_size"]),
            "--tm_delta {}".format(kwargs["tm_delta"]),
            "--seed {}".format(kwargs["seed"]),
            "--offtarget_min_size {}".format(kwargs["offtarget_min_size"]),
            "--offtarget_max_size {}".format(kwargs["offtarget_max_size"]),
            "-nt 1",
            "--webapp",
        ]

        if include_vcf:
            # TODO: upload proper file
            strcmd += ["--vcf {}".format(join(app.config["VCF_REPO"], "vcf.txt.gz"))]
            strcmd += ["--vcf_include {}".format(join(dest_folder, "vcf_include.txt"))]

        if kwargs["enzymes"]:
            strcmd += ["--enzymes {}".format(kwargs["enzymes"])]

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

        # BlastDB
        kwargs["reference"] = request.form.get("reference")

        # todo rename depth
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            exe,
            kwargs["roi"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-nt 1",
            "--webapp",
        ]

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
    filename = join(app.config['UPLOAD_FOLDER'], task_id, "output.txt")
    return send_file(filename)


if __name__ == "__main__":
    app.run(debug=True)
