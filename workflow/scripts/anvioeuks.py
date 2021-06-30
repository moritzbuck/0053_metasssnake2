import sys, os
sys.path.append(os.getcwd())
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO

script , binning_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "binning.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/binning.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/binning_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


title2log("copying assemblies to temp_folder", logfile)

def anvio_completeness(fna):
    os.makedirs("temp")
    shutil.copy(fna, "temp/genome.fna")
    call("""
    anvi-setup-scg-taxonomy -T {threads} 2> /dev/null
    anvi-script-reformat-fasta temp/genome.fna -o temp/contigs4anvio.fa -l 0 --simplify-names 2> /dev/null
    anvi-gen-contigs-database -f temp/contigs4anvio.fa -o temp/contigs.db -T {threads} 2> /dev/null
    anvi-run-hmms -c temp/contigs.db -T {threads} 2> /dev/null
    anvi-run-scg-taxonomy -c temp/contigs.db -T {threads} 2> /dev/null
    anvi-estimate-genome-completeness -c  temp/contigs.db -o temp/completeness.txt 2> /dev/null
    """.format(threads = 24), shell=True)
    with open("temp/completeness.txt") as handle:
        lines = list(handle.readlines())
    shutil.rmtree("temp")
    return {'domain' : lines[1].split()[1], 'completeness' : float(lines[1].split()[3]), 'redundancy' : float(lines[1].split()[4]), 'confidence' : float(lines[1].split()[2]) }

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file), ass = pjoin(out_folder, "binned_assembly.fna")), shell=True)

title2log("Binning done", logfile)

shutil.rmtree(temp_folder)
