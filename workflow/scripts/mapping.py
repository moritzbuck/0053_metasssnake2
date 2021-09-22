import sys, os
sys.path.append(os.getcwd())
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line, csv2dict
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO
import pandas

script , mapping_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "mapping.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/mapping.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/mapping_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


subset = config_file['mappings'][mapping_name]['subset']
method = config_file['mappings'][mapping_name]['mapper']
binset = config_file['mappings'][mapping_name]['binset']
taxfield = config_file['mappings'][mapping_name]['taxfield']

ani = config_file['mappings'][mapping_name]['min_nucleotide_id']
threads = int(threads)

temp_folder = pjoin(config_file['temp_folder'], "mappings", mapping_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("copying binset to temp_folder", logfile)

shutil.copy(pjoin(root_folder, "binsets", binset, binset + ".fna"), pjoin(temp_folder, "binset.fna"))

title2log("indexing binset to temp_folder", logfile)

call("bwa-mem2 index {temp}/binset.fna >> {out_folder}/logs/mapping.log  2>&1".format(temp = temp_folder, threads = threads, out_folder = out_folder), shell=True)

freetxt_line("Starting mappings", logfile)
coverages = {}
total_reads = {}
for lib in config_file['mappings'][mapping_name]['libraries']:
    title2log("copying {subset}{lib} to temp_folder".format(lib = lib, subset = "subset_" if subset else ""), logfile)
    call("""
    unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}_fwd.fastq.gz > {temp}/fwd.fastq 2>> {out_folder}/logs/mapping.log
    unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}_rev.fastq.gz > {temp}/rev.fastq 2>> {out_folder}/logs/mapping.log
    unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}_unp.fastq.gz > {temp}/unp.fastq 2>> {out_folder}/logs/mapping.log
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder, out_folder = out_folder, subset_str = "/subs/subs_" if subset else "/"), shell=True)
    title2log("mapping {lib} to ref".format(lib = lib), logfile)
    call("""
    bwa-mem2 mem -t {threads} {temp}/binset.fna  -o {temp}/mapping.sam {temp}/fwd.fastq {temp}/rev.fastq 2>> {out_folder}/logs/mapping.log
    samtools view  -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_pairs.bam - >> {out_folder}/logs/mapping.log 2>&1
    bwa-mem2 mem -t {threads} {temp}/binset.fna  -o {temp}/mapping.sam {temp}/unp.fastq 2>> {out_folder}/logs/mapping.log
    samtools view -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_unpaired.bam - >> {out_folder}/logs/mapping.log 2>&1
    rm {temp}/mapping.sam 2>> {out_folder}/logs/mapping.log
    samtools merge -f -t {threads} {temp}/mapping.bam  {temp}/mapping_pairs.bam  {temp}/mapping_unpaired.bam >> {out_folder}/logs/mapping.log  2>&1
    rm {temp}/mapping_pairs.bam {temp}/mapping_unpaired.bam 2>> {out_folder}/logs/mapping.log
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder, threads = threads, out_folder = out_folder), shell=True)

    call("""
    coverm filter -b {temp}/mapping.bam -o {temp}/mapping_filtered.bam --min-read-percent-identity {ani} --threads {threads} >> {out_folder}/logs/mapping.log  2>&1
    rm {temp}/mapping.bam >> {out_folder}/logs/mapping.log  2>&1
    coverm contig  --bam-files {temp}/mapping_filtered.bam  --methods count  --threads {threads} > {temp}/mapping.tsv 2>> {out_folder}/logs/mapping.log
    rm {temp}/mapping_filtered.bam
    cat {temp}/fwd.fastq {temp}/rev.fastq {temp}/unp.fastq | wc -l > {temp}/total_reads_x4.txt
    """.format(temp=temp_folder, threads = threads, ani = ani,  out_folder = out_folder), shell=True)

    call("""
rm {temp}/fwd.fastq
rm {temp}/rev.fastq
rm {temp}/unp.fastq
""".format(temp = temp_folder), shell = True)

    with open(pjoin(temp_folder, "mapping.tsv")) as handle:
        handle.readline()
        coverages[lib] = {}
        for l in handle:
            ll = l.strip().split()
            coverages[lib][ll[0]] = float(ll[1])
    with open(pjoin(temp_folder, "total_reads_x4.txt")) as handle:
        total_reads[lib] = int(handle.readline().strip())/4

title2log("Done with the mappings", logfile)

title2log("Making tables", logfile)

for k in total_reads:
    coverages[k]['unmapped'] = total_reads[k] - sum(coverages[k].values())

with open(pjoin(temp_folder, "total_reads_to_map.csv"), "w") as handle:
    handle.writelines(["library,total_reads"] + [f"{k},{v}\n" for k,v in total_reads.items()])

coverages = pandas.DataFrame.from_dict(coverages)
(coverages/coverages.sum()).to_csv(pjoin(temp_folder, "contigs_relative_abundance.csv"), index_label = "contig_name")

coverages['bin'] = [c.split("_ctg")[0] for c in coverages.index]
coverages_by_bin = coverages.groupby("bin").sum()

relab_by_bin = (coverages_by_bin/coverages_by_bin.sum())
relab_by_bin.to_csv(pjoin(temp_folder, "bins_relative_abundance.csv"), index_label = "bin_name")

levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
bin_stats  = csv2dict(pjoin(root_folder, "binsets", binset, binset + "_basics.csv" ))


taxas = [bin_stats.get(bin_, {taxfield : 'unmapped'})[taxfield] for bin_ in relab_by_bin.index]
for l in range(7):
    relabs = relab_by_bin.copy()
    relabs[levels[l]] = [";".join(t.split(";")[0:(l+1)]) for t in taxas]
    relabs = relabs.groupby(levels[l]).sum()
    relabs.to_csv(pjoin(temp_folder, levels[l] + "_relative_abundance.csv"), index_label = levels[l])

relabs = relab_by_bin.copy()
relabs["mOTU"] = [bin_stats.get(bin_, {'mOTU' : 'unmapped'})['mOTU'] for bin_ in relabs.index]
relabs = relabs.groupby("mOTU").sum()
relabs.to_csv(pjoin(temp_folder, "mOTU_relative_abundance.csv"), index_label = levels[l])


call(f"cp {temp_folder}/*.csv {root_folder}/mappings/{mapping_name}/", shell=True)

title2log("Cleaning up and moving the bins", logfile)


shutil.rmtree(temp_folder)
