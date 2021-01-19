import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO

script , binset_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "binset_name.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/binset.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/binning_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)

min_complete = config_file['binsets'][binset_name]['min_complete']
max_contamin = config_file['binsets'][binset_name]['max_contamin']
min_size = config_file['binsets'][binset_name]['min_size']
min_coding = config_file['binsets'][binset_name]['min_coding']

temp_folder = pjoin(config_file['temp_folder'], "binnings", binning_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("copying bins to temp_folder", logfile)

binnings = [pjoin(root_folder, "binnings", binni, "bins") for binni in config_file['binsets'][binset_name]['binnings']]
for binni in binnings:
    call("cp {binni}/* {temp}/bins/".format(binni = binni, temp = temp_folder), shell = True)

title2log("Cleaning up metabat bins", logfile)
tbinfoder = "{temp}/bins".format(temp = temp_folder)
cbinfoder = "{temp}/clean_bins".format(temp = temp_folder)
os.makedirs(cbinfoder, exist_ok = True)
bins = [b for b in os.listdir(tbinfoder) if b.endswith(".fa")]
zeros = len(str(len(bins)))
for b in tqdm(bins):
    clean_bin(pjoin(tbinfoder, b), pjoin(cbinfoder, binning_name + "_bin-" + b.split(".")[-2].zfill(zeros) + ".fna"), binning_name + "_bin-" + b.split(".")[-2].zfill(zeros))

title2log("Cleaning up and moving the bins", logfile)

if os.path.exists(pjoin(out_folder, "binned_assembly.fna")):
    os.remove(pjoin(out_folder, "binned_assembly.fna"))

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file), ass = pjoin(out_folder, "binned_assembly.fna")), shell=True)

title2log("Binning done", logfile)

shutil.rmtree(temp_folder)
