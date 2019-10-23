import os
import shutil
from os.path import join

TEMP_DIR = join(os.getcwd(), "temporary")

# Cleanup
os.remove("out.txt")
os.remove("out.bed")
os.remove("out.log")
# os.remove("out_altlist.txt")
shutil.rmtree(TEMP_DIR)
