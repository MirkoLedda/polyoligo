import numpy as np

N_CPU = 4

max_mem = {}

for i in range(N_CPU):
    j = i + 1
    filename = "footprinting_data/mprof_{}.dat".format(j)
    max_mem[j] = 0

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("MEM"):
                mem = float(line.strip().split()[1])
                if mem > max_mem[j]:
                    max_mem[j] = mem

print(max_mem)

d = []
for i in range(N_CPU-1):
    j = i+1
    d.append(max_mem[j+1] - max_mem[j])

d = np.mean(d)
print(max_mem[1]-d)
print(d)

with open("footprinting_data/mprof.dat", "w") as f:
    for i in range(N_CPU):
        j = i + 1
        f.write("{} {}\n".format(j, max_mem[j]))
