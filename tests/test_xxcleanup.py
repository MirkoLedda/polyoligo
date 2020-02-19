import os
import shutil
from os.path import join, exists

TEMP_DIR = join(os.getcwd(), "temporary")

# Cleanup
try:
    os.remove("out.txt")
except FileNotFoundError:
    pass
try:
    os.remove("out.bed")
except FileNotFoundError:
    pass
try:
    os.remove("out.log")
except FileNotFoundError:
    pass
try:
    os.remove("out_altlist.txt")
except FileNotFoundError:
    pass

if exists(TEMP_DIR):
    shutil.rmtree(TEMP_DIR)

for i in range(100):
    if exists(TEMP_DIR + str(i)):
        shutil.rmtree(TEMP_DIR + str(i))
