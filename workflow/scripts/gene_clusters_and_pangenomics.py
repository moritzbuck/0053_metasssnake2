import sys, os
from tqdm import tqdm
from os.path import join as pjoin
sys.path.append(os.getcwd())

from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file, csv2dict
import shutil
from subprocess import call
import json


script , binset_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "gene_clusters.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/gene_clusters.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/gene_clusters_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)

threads = int(threads)
seqid = 0.0
covmode = 0
cov = 0.8
boots = 10
temp_folder = pjoin(config_file['temp_folder'], "binsets", binset_name)
gc_folder = pjoin(temp_folder, "gene_clusters")
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)
os.makedirs(gc_folder, exist_ok=True)

title2log("copying bins to temp_folder", logfile)

os.makedirs(pjoin(temp_folder, "bins") , exist_ok = True)
for bin_ in tqdm(os.listdir(pjoin(out_folder, "bins"))):
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".faa"), pjoin(temp_folder, "bins", bin_ + ".faa"))
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".db"), pjoin(temp_folder, "bins", bin_ + ".db"))
    call(f"cat {pjoin(temp_folder, 'bins', bin_ + '.faa')} >> {temp_folder}/all_proteoms.faa", shell = True)


binset_stats = csv2dict(pjoin(root_folder, "binsets", binset_name, binset_name + "_basics.csv"))
call(f"mmseqs easy-cluster --min-seq-id {seqid} --cov-mode {covmode} -c {cov} --threads {threads} {temp_folder}/all_proteoms.faa {temp_folder}/mmseqs_out {temp_folder}/mmseqs_temp > {logfile}", shell=True)

with open(f"{temp_folder}/mmseqs_out_cluster.tsv") as handle:
    recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

#pretty formating names
fill = len(str(len(set(recs.values()))))

rep2clust = {k : binset_name + "_GC_" + str(i).zfill(fill) for i,k in enumerate(set(recs.values()))}
gene_clusters2rep = {v: k for k,v in rep2clust.items()}

freetxt_line(" ".join(["For", str(len(recs)), "CDSes we got ", str(len(gene_clusters2rep)), " gene-clusters"]), logfile)

recs = {k : rep2clust[v] for k, v in recs.items()}
clst2gene = { clst : [] for clst in gene_clusters2rep}
for gene, clst in recs.items():
    clst2gene[clst].append(gene)

bin2gc = {k : [] for k in binset_stats}
for gc, cdss in clst2gene.items():
    for cds in cdss:
        bin_ = "_".join(cds.split("_")[:-1])
        bin2gc[bin_] += [gc]

with open(pjoin(gc_folder, "gene_cluster2representative_cds.json") , "w") as handle:
    json.dump(gene_clusters2rep, handle, indent = 2, sort_keys = True)

with open(pjoin(gc_folder, "bin2gene_clusters.json") , "w") as handle:
    json.dump(bin2gc, handle, indent = 2, sort_keys = True)

with open(f"{temp_folder}/checkm.txt", "w") as handle :
    handle.writelines(["Bin Id\tCompleteness\tContamination\n"] + [ f"{k}\t{v['percent_completion']}\t{v['percent_redundancy']}\n" for k,v in binset_stats.items()])

all_motus = {v['mOTU'] : [] for v in binset_stats.values() if v['mOTU']}
for k,v in binset_stats.items():
    if v['mOTU']:
        all_motus[v['mOTU']] += [k]

os.makedirs(f"{gc_folder}/mOTUs", exist_ok = True)
for motu, bins_ in tqdm(all_motus.items()):
    bins_ = set(bins_)
    if len(bins_) > 2:
        with open(pjoin(temp_folder, "bin2gc.json") , "w") as handle:
            json.dump({k : v for k,v in bin2gc.items() if k in bins_}, handle, indent = 2, sort_keys = True)
        call(f"mOTUpan.py --name {motu} --gene_clusters_file {temp_folder}/bin2gc.json  --checkm {temp_folder}/checkm.txt  -o {gc_folder}/mOTUs/{motu}.tsv --boots {boots}  >>{logfile} 2>&1 ", shell = True)

title2log("pulling anvi'o dataSets for annot if needed", logfile)

call(f"""
anvi-setup-kegg-kofams --quiet >> {logfile}  2>&1
""", shell = True)


title2log("Cleaning up and moving the bins", logfile)


shutil.rmtree(temp_folder)
