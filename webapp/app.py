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


# CELERY -----------------------------------------------------------------------
@celery.task()
def task_kasp(strcmd, log_dest):
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


@celery.task()
def task_crispr(strcmd, log_dest):
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


@app.route('/pcr')
def pcr():
    return render_template('pcr.html')


@app.route('/kasp', methods=['GET', 'POST'])
def kasp():
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

        # Selected populations
        if include_vcf:
            selected_pops = []
            for i in range(9):
                current_pop = "population{}".format(i)
                fetched_pop = request.form.get(current_pop)

                if fetched_pop:
                    with open(os.path.abspath(join("./static/data", fetched_pop)), "r") as f:
                        for line in f:
                            selected_pops.append(line.strip())

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

        # todo rename depth
        kwargs["multiplier"] = kwargs["depth"]
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            "polyoligo-kasp",
            kwargs["markers"],
            join(dest_folder, "output"),
            kwargs["reference"],
            "-n {}".format(kwargs["n"]),
            "--multiplier {}".format(kwargs["multiplier"]),
            "--dye1 {}".format(kwargs["dye1"]),
            "--dye2 {}".format(kwargs["dye2"]),
            "--tm_delta {}".format(kwargs["tm_delta"]),
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

        task_kasp.delay(strcmd, dest_folder)

        return redirect(url_for('processing', task_id=task_id))
    else:
        return render_template('kasp.html', blastdb_options=app.config["BLASTDB_OPTIONS"])


@app.route('/crispr', methods=['GET', 'POST'])
def crispr():
    if request.method == 'POST':

        kwargs_names = ["roi", "reference"]
        kwargs = {}
        for kwargs_name in kwargs_names:
            kwargs[kwargs_name] = None

        task_id, dest_folder = create_task_dest()
        os.makedirs(dest_folder)

        # Marker file
        kwargs["roi"] = request.form.get("roi")

        # BlastDB
        kwargs["reference"] = request.form.get("reference")

        # todo rename depth
        kwargs["reference"] = join(app.config["BLASTDB_REPO"], "blastdb")

        strcmd = [
            "polyoligo-crispr",
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

        task_crispr.delay(strcmd, dest_folder)

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
