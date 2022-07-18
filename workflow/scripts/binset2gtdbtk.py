import sys, os
from tqdm import tqdm
from os.path import join as pjoin
sys.path.append(os.getcwd())

from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file, csv2dict
import shutil
from subprocess import call
import json

script , binset_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "gtdbtk.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/gtdbtk.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/gtdbtk_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


temp_folder = pjoin(config_file['temp_folder'], "gtdbtk", binset_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)
os.makedirs(temp_folder, exist_ok=True)

title2log("copying bins to temp_folder", logfile)

os.makedirs(pjoin(temp_folder, "bins") , exist_ok = True)
for bin_ in tqdm(os.listdir(pjoin(out_folder, "bins"))):
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".fna"), pjoin(temp_folder, "bins", bin_ + ".fna"))

call(f"gtdbtk classify_wf --out_dir {temp_folder}/gtdbtk --genome_dir {temp_folder}/bins/ -x fna --cpus {threads}  >> {logfile}", shell = True)

classif = csv2dict(f"{temp_folder}/gtdbtk/gtdbtk.ar53.summary.tsv", sep = "\t")
classif.update(csv2dict(f"{temp_folder}/gtdbtk/gtdbtk.bac120.summary.tsv", sep = "\t"))

binset_stats = csv2dict(pjoin(root_folder, "binsets", binset_name, binset_name + "_basics.csv"))
cleanz = {k : {'gtdbtk_classif' : v['classification'], 'gtdbtk_notes' : ";".join([ field + "=" + v[field].replace(" ","_").replace(",","_") for field in ['note','classification_method','warnings'] if v[field] != "N/A"]), 'translation_table' : v["translation_table"]} for k,v in classif.items()}

for k,v in cleanz.items():
    binset_stats[k].update(v)


title2log("Cleaning up and moving the bins", logfile)
call(f"rm {temp_folder}/gtdbtk/gtdbtk.*", shell=True)
call(f"mv {temp_folder}/gtdbtk/*/gtdbtk.* {temp_folder}/gtdbtk/", shell=True
)
shutil.rmtree(pjoin(temp_folder,"gtdbtk", "align"))
shutil.rmtree(pjoin(temp_folder,"gtdbtk", "identify"))
shutil.rmtree(pjoin(temp_folder,"gtdbtk", "classify"))
shutil.move(pjoin(temp_folder,"gtdbtk"), out_folder)

dict2file(binset_stats, pjoin(out_folder, binset_name + "_basics.csv"))
title2log("Binsetting done", logfile)

shutil.rmtree(temp_folder)
