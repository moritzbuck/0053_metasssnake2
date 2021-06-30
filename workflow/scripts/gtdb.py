import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file
import shutil
from subprocess import call
from os.path import basename
import json
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

call("checkm lineage_wf -x fna -t 24 {temp}/bins/ {temp}/checkm > {temp}/checkm.txt  2>> {logfile}".format(temp = temp_folder, logfile = logfile), shell = True)

with open(pjoin(temp_folder, "checkm.txt")) as handle:
    all_lines = [l.strip() for l in  handle.readlines() if " INFO:" not in l]

title2log("Parsing and filtering based on checkm", logfile)

all_lines = [re.sub(r"  +","\t", a).split("\t") for a in all_lines]
header_lines = [i for i,l in enumerate(all_lines) if 'Bin Id' in l and 'Completeness' in l and 'Contamination' in l]
header_lines = header_lines[0]
header_line = all_lines[header_lines]
lines = [l for i,l in enumerate(all_lines) if i != header_lines and len(l) == len(header_line)]
lines = [{a : b if a in ['Marker lineage', 'Bin Id'] else float(b) for a,b in zip(header_line,l) }for l in lines]

checkm_out = {l['Bin Id'] : {k: l[k] for k in ('Completeness', 'Contamination', 'Strain heterogeneity')} for l in lines if l['Completeness'] > min_completeness and l['Contamination'] < max_contamination}

for f in os.listdir(tbinfoder):
    if f[:-4] not in checkm_out:
        if keep_fails:
            append2unkept(f)
        os.remove(pjoin(tbinfoder, f))

title2log("Parsing and filtering based faa/fnas", logfile)

for f in os.listdir(tbinfoder):
    fna = list(SeqIO.parse(pjoin(tbinfoder, f), "fasta"))
    faa = list(SeqIO.parse(pjoin(temp_folder, "checkm/bins/", f[:-4], "genes.faa"), "fasta"))
    checkm_out[f[:-4]]['length'] = sum([len(s) for s in fna])
    checkm_out[f[:-4]]['acoding_density'] = 3*sum([len(s) for s in faa])/checkm_out[f[:-4]]['length']

for f in os.listdir(tbinfoder):
    if checkm_out[f[:-4]]['acoding_density'] < min_coding:
        if keep_fails:
            append2unkept(f)
        os.remove(pjoin(tbinfoder, f))
        del checkm_out[f[:-4]]


title2log("Running prokka for the bins", logfile)


call("sed -i 's/my $MAXCONTIGIDLEN = 37/my $MAXCONTIGIDLEN = 250/' `which prokka` ", shell=True)
call("sed -i 's/[^#]tbl2asn -V/#tbl2asn -V/' `which prokka` ", shell=True)


call("prokka --outdir {temp}/clean_bins/{binset}_unkept --prefix {binset}_unkept --locustag {binset}_unkept --metagenome --cpus {threads} {temp}/clean_bins/{binset}_unkept.fna".format(threads= threads, binset = binset_name, temp=temp_folder), shell=True)
os.remove(pjoin(cbinfoder, binset_name + "_unkept.fna"))

call("ls {temp}/bins/ | cut -f1 -d. | parallel -j{threads} prokka --outdir {temp}/clean_bins/{{}} --prefix {{}} --locustag {{}} --cpus 1 {temp}/bins/{{}}.fna".format(threads= threads, temp=temp_folder), shell = True)

for g in checkm_out:
    fna = [s.seq for s in SeqIO.parse(pjoin(cbinfoder, g, g +".fna"), "fasta")]
    faa = [s.seq for s in SeqIO.parse(pjoin(cbinfoder, g, g +".faa"), "fasta")]
    checkm_out[g]['nb_contigs'] = len(fna)
    checkm_out[g]['bin_len'] = sum([len(s) for s in fna])
    checkm_out[g]['nb_cdss'] = len(faa)
    checkm_out[g]['acoding_density'] = sum([len(s)*3 for s in faa])/checkm_out[g]['bin_len']


dict2file(checkm_out, pjoin(out_folder, binset_name + "_basics.csv"))

title2log("Cleaning up and moving the bins", logfile)

if os.path.exists(pjoin(out_folder, binset_name + ".fna")):
    os.remove(pjoin(out_folder,  binset_name +  ".fna"))

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    if file != binset_name + "_unkept":
        call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".fna"), ass = pjoin(out_folder,  binset_name + ".fna")), shell=True)
        call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".faa"), ass = pjoin(out_folder,  binset_name + ".faa")), shell=True)
        call("cat {file} | grep -v '^##' >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".gff"), ass = pjoin(out_folder,  binset_name + ".gff")), shell=True)

title2log("Binsetting done", logfile)

shutil.rmtree(temp_folder)
