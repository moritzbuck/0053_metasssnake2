import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO
import re

script , binset_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "binset.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/binset.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/binning_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)

min_completeness = config_file['binsets'][binset_name]['min_completeness']
max_contamination = config_file['binsets'][binset_name]['max_contamination']
min_size = config_file['binsets'][binset_name]['min_size']
min_coding = config_file['binsets'][binset_name]['min_coding']
keep_fails = config_file['binsets'][binset_name]['keep_fails']

temp_folder = pjoin(config_file['temp_folder'], "binset", binset_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("copying bins to temp_folder", logfile)

binnings = [pjoin(root_folder, "binnings", binni, "bins") for binni in config_file['binsets'][binset_name]['binnings']]
os.makedirs(pjoin(temp_folder, "bins") , exist_ok = True)
for binni in binnings:
    title2log("copying bins from " + binni.split("/")[-2], logfile)
    call("cp {binni}/* {temp}/bins/".format(binni = binni, temp = temp_folder), shell = True)

tbinfoder = "{temp}/bins".format(temp = temp_folder)
cbinfoder = "{temp}/clean_bins".format(temp = temp_folder)
os.makedirs(cbinfoder, exist_ok = True)

title2log("Handle unbinned", logfile)

append2unkept = lambda f : call("cat {folder}/{f} >> {ofolder}/{binset}_unkept.fna".format(f = f, folder = tbinfoder, ofolder = cbinfoder, binset = binset_name), shell=True)

for f in os.listdir(tbinfoder):
    if f.endswith("-unbinned.fna"):
        if keep_fails:
            append2unkept(f)
        os.remove(pjoin(tbinfoder, f))

title2log("Running checkm", logfile)

#call("checkm lineage_wf -x fna -t 24 {temp}/bins/ {temp}/checkm > {temp}/checkm.txt  2>> {logfile}".format(temp = temp_folder, logfile = logfile), shell = True)

with open(pjoin(temp_folder, "checkm.txt")) as handle:
    all_lines = [l.strip() for l in  handle.readlines() if " INFO:" not in l]


all_lines = [re.sub(r"  +","\t", a).split("\t") for a in all_lines]
header_lines = [i for i,l in enumerate(all_lines) if 'Bin Id' in l and 'Completeness' in l and 'Contamination' in l]
header_lines = header_lines[0]
header_line = all_lines[header_lines]
lines = [l for i,l in enumerate(all_lines) if i != header_lines and len(l) == len(header_line)]
lines = [{a : b if a in ['Marker lineage', 'Bin Id'] else float(b) for a,b in zip(header_line,l) }for l in lines]

checkm_out = {l['Bin Id'] : {k: l[k] for k in ('Completeness', 'Contamination')} for l in lines if l['Completeness'] > min_completeness and l['Contamination'] < max_contamination}

for f in os.listdir(tbinfoder):
    if f[:-4] not in checkm_out:
        if keep_fails:
            append2unkept(f)
        os.remove(pjoin(tbinfoder, f))



syst.exit(0)

bins = [b for b in os.listdir(tbinfoder) if b.endswith(".fa")]
zeros = len(str(len(bins)))
for b in tqdm(bins):
    clean_bin(pjoin(tbinfoder, b), pjoin(cbinfoder, binset_name + "_bin-" + b.split(".")[-2].zfill(zeros) + ".fna"), binset_name + "_bin-" + b.split(".")[-2].zfill(zeros))

title2log("Cleaning up and moving the bins", logfile)

if os.path.exists(pjoin(out_folder, "binned_assembly.fna")):
    os.remove(pjoin(out_folder, "binned_assembly.fna"))

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file), ass = pjoin(out_folder, "binned_assembly.fna")), shell=True)

title2log("Binsetting done", logfile)

shutil.rmtree(temp_folder)
