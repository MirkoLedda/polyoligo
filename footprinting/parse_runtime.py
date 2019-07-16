import datetime
import numpy as np

dat = {}

with open("runtime.txt", "r") as f:

    while True:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()

        if not line3:
            break

        m = int(line1.split("=")[1].strip())
        n = int(line2.split()[5].strip())
        t = line3.split()[3]

        t = datetime.datetime.strptime(t, "%H:%M:%S")
        delta = datetime.timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
        t = int(delta.total_seconds())

        if m not in dat.keys():
            dat[m] = {}

        if n not in dat[m].keys():
            dat[m][n] = [t]
        else:
            dat[m][n].append(t)

with open("runtime.dat", "w") as f:
    for m in dat.keys():
        for n in dat[m].keys():
            for t in dat[m][n]:
                f.write("{} {} {}\n".format(m, n, t))
