import sys, os
from tqdm import tqdm
from os.path import join as pjoin
sys.path.append(os.getcwd())

from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file, csv2dict
import shutil
from subprocess import call
import json
from collections import OrderedDict
import anvio
import pandas

script , binset_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "gene_clusters.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/gene_clusters.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/gene_clusters_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)

threads = int(threads)
seqid = 0.7
covmode = 1
cov = 0.8
boots = 10
min_genomes = 10
remove_singleton_gcs = True
module_cutoff = 0.75
temp_folder = pjoin(config_file['temp_folder'], "binsets", binset_name)
gc_folder = pjoin(temp_folder, "gene_clusters")
kegg_folder = pjoin(temp_folder, "KEGGs")

freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)
os.makedirs(gc_folder, exist_ok=True)

title2log("copying bins to temp_folder", logfile)

os.makedirs(pjoin(temp_folder, "bins") , exist_ok = True)
for bin_ in tqdm(os.listdir(pjoin(out_folder, "bins"))):
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".gff"), pjoin(temp_folder, "bins", bin_ + ".gff"))
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".faa"), pjoin(temp_folder, "bins", bin_ + ".faa"))
    shutil.copyfile(pjoin(out_folder, "bins", bin_, bin_ + ".db"), pjoin(temp_folder, "bins", bin_ + ".db"))

title2log("copying big faa to temp_folder", logfile)
shutil.copyfile(pjoin(out_folder, binset_name + ".faa"), pjoin(temp_folder, "all_proteoms.faa"))

title2log("Running mmseqs easy-cluster for gene clusters", logfile)
freetxt_line(f" --min-seq-id {seqid} --cov-mode {covmode} -c {cov}", logfile)

binset_stats = csv2dict(pjoin(root_folder, "binsets", binset_name, binset_name + "_basics.csv"))
call(f"mmseqs easy-cluster --min-seq-id {seqid} --cov-mode {covmode} -c {cov} --threads {threads} {temp_folder}/all_proteoms.faa {temp_folder}/mmseqs_out {temp_folder}/mmseqs_temp > {logfile}", shell=True)

title2log("Parsing the output", logfile)

with open(f"{temp_folder}/mmseqs_out_cluster.tsv") as handle:
    recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

#pretty formating names
fill = len(str(len(set(recs.values()))))

rep2clust = {k : binset_name + "_GC_" + str(i).zfill(fill) for i,k in enumerate(set(recs.values()))}
gene_clusters2rep = {v: k for k,v in rep2clust.items()}

freetxt_line(" ".join(["For", str(len(recs)), "CDSes we got ", str(len(gene_clusters2rep)), " gene-clusters"]), logfile)

clst2gene = { clst : [] for clst in gene_clusters2rep}
recs = {k : rep2clust[v] for k, v in recs.items()}
for gene, clst in recs.items():
    clst2gene[clst].append(gene)

bin2gc = {k : [] for k in binset_stats}
for gc, cdss in clst2gene.items():
    for cds in cdss:
        bin_ = "_".join(cds.split("_")[:-1])
        bin2gc[bin_] += [gc]

title2log("Editing the gff-files and dumping data", logfile)
for bin_ in tqdm([f for f in os.listdir(pjoin(temp_folder, "bins")) if f.endswith(".gff")]):
    with open(pjoin(temp_folder,"bins", bin_)) as handle:
        lines = handle.readlines()
    for i in range(len(lines)):
        ll = lines[i].strip()
        if ll=="##FASTA":
            break
        if ll[0] != "#" and "CDS" in ll:
            ll = ll.split("\t")
            tail = ll[-1]
            tail = OrderedDict([tuple(kv.split("=")) for kv in tail.strip().split(";")])
            tail['gene_cluster']  =  recs[tail['ID']]
            tail = ";".join([k + "=" + v for k,v in tail.items()])
            ll[-1] = tail
            ll = "\t".join(ll)+ "\n"
            lines[i] = ll
    with open(pjoin(temp_folder,"bins", bin_), "w") as handle:
        handle.writelines(lines)
    call("sed '/##FASTA/q' {file} | grep -v '^# ' >> {ass}".format(file = pjoin(temp_folder,"bins", bin_), ass = pjoin(temp_folder,  binset_name + ".gff")), shell=True)


with open(pjoin(gc_folder, "gene_cluster2representative_cds.json") , "w") as handle:
    json.dump(gene_clusters2rep, handle, indent = 2, sort_keys = True)

with open(pjoin(gc_folder, "bin2gene_clusters.json") , "w") as handle:
    json.dump(bin2gc, handle, indent = 2, sort_keys = True)

with open(pjoin(gc_folder, "gene2gene_clusters.json") , "w") as handle:
    json.dump(recs, handle, indent = 2, sort_keys = True)


with open(f"{temp_folder}/checkm.txt", "w") as handle :
    handle.writelines(["Bin Id\tCompleteness\tContamination\n"] + [ f"{k}\t{v['percent_completion']}\t{v['percent_redundancy']}\n" for k,v in binset_stats.items()])

all_motus = {v['mOTU'] : [] for v in binset_stats.values() if v['mOTU']}
for k,v in binset_stats.items():
    if v['mOTU']:
        all_motus[v['mOTU']] += [k]

title2log(f"Running mOTUpan on mOTUs >{min_genomes} genomes", logfile)
freetxt_line(f" and parsing mOTUpan-outs, removing singleton GCs : {remove_singleton_gcs}", logfile)

gc_pangenome_stats = {}
os.makedirs(f"{gc_folder}/mOTUs", exist_ok = True)
for motu, bins_ in tqdm(all_motus.items()):
    bins_ = set(bins_)
    if len(bins_) > min_genomes:
        with open(pjoin(temp_folder, "bin2gc.json") , "w") as handle:
            json.dump({k : v for k,v in bin2gc.items() if k in bins_}, handle, indent = 2, sort_keys = True)
        call(f"mOTUpan.py --name {motu} --gene_clusters_file {temp_folder}/bin2gc.json  --checkm {temp_folder}/checkm.txt  -o {gc_folder}/mOTUs/{motu}.tsv --boots {boots}  >>{logfile} 2>&1 ", shell = True)
        with open(f"{gc_folder}/mOTUs/{motu}.tsv") as handle:
            header = {}
            main = []
            for l in handle:
                l = l.strip()
                if l != "#" and not l.startswith("#mOTUlizer"):
                    if l.startswith("#"):
                        header[l.split("=")[0][1:]] = "=".join(l.split("=")[1:])
                    else :
                        main += [l]
            cols = main[0].split()
            gcs = [{k : ll for k,ll in zip(cols, l.split())} for l in main[1:] if l != '']
            if remove_singleton_gcs:
                gcs = [d for d in gcs if int(d['genome_occurences']) > 1]
            core = []
            accessory = []
            for d in gcs:
                if d['type'] == 'core':
                    core.append(d['trait_name'])
                else :
                    accessory.append(d['trait_name'])
                    for dat in header['genomes'].split(";"):
                        dat = dat.split(":")
                        bin_ = dat[0]
                        new_comp = [dd for dd in dat[1:] if dd.startswith("posterior_complete=")][0]
                        binset_stats[bin_]['mOTUpan_completeness'] = float(new_comp.replace("posterior_complete=",""))
                        gc_pangenome_stats[motu] = {
                        'genome_count' : int(header['genome_count']),
                        'core_length': int(header['core_length']),
                        'mean_prior_completeness' : float(header['mean_prior_completeness']),
                        'mean_posterior_completeness' : float(header['mean_posterior_completeness']),
                        "nb_bootstraps" : 10 if "bootstrapped_mean_false_positive_rate" not in header else int(header['bootstrapped_nb_reps']),
                        "fpr" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[0]),
                        "sd_fpr" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[1].split("=")[1]),
                        "lowest_false" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[0]),
                        "sd_lowest_false" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[1].split("=")[1]),
                        'core' : ";".join(core),
                        'accessory' : ";".join(accessory)
                        }

title2log("pulling anvi'o dataSets for annot if needed", logfile)

call(f"""
anvi-setup-kegg-kofams --quiet >> {logfile}  2>&1
""", shell = True)

with open(f"{temp_folder}/bin_list.tsv", "w") as handle :
    handle.writelines(["name\tcontigs_db_path\n"] + [ f"{k}\t{temp_folder}/bins/{k}.db\n" for k,v in binset_stats.items()])


title2log(f"Running KEGG module estimation for all genoems", logfile)
bin2ko = {}
for bin_ in tqdm([f for f in os.listdir(pjoin(temp_folder, "bins")) if f.endswith(".db")]):
    out = bin_.replace(".db", ".kofam")
    if bin_.replace(".db", "") not in bin2ko:
        call(f"anvi-export-functions -c {temp_folder}/bins/{bin_} --annotation-sources KOfam -o {temp_folder}/bins/{out} --quiet", shell = True)
        with open(f"{temp_folder}/bins/{out}") as handle:
            kos = set([l.split()[2] for l in handle.readlines()])
            bin2ko[bin_.replace(".db", "")] = kos
call(f"anvi-estimate-metabolism -e{temp_folder}/bin_list.tsv --output-file-prefix {temp_folder}/KEGG-metabolism --matrix-format --quiet", shell = True)



module2bin = {}
with open(f"{temp_folder}/KEGG-metabolism-presence-MATRIX.txt") as handle:
    header = handle.readline()
    header = header.strip().split()[1:]
    for l in handle:
        dat = l.strip().split()
        module2bin[dat[0]] = [bin_ for bin_, val in zip(header,dat[1:]) if val == '1' ]

bin2module = {bin_ : []for bin_ in binset_stats}
for module, bins_ in module2bin.items():
    for bin_ in bins_:
        bin2module[bin_].append(module)

os.makedirs(kegg_folder, exist_ok = True)

with open(pjoin(kegg_folder, "bin2module.json") , "w") as handle:
    json.dump(bin2module, handle, indent = 2, sort_keys = True)

with open(pjoin(kegg_folder, "bin2ko.json") , "w") as handle:
    json.dump({k: list(v) for k,v in bin2ko.items()}, handle, indent = 2, sort_keys = True)


title2log(f"Running KEGG mOTUpan on mOTUs >{min_genomes} genomes", logfile)
freetxt_line(f" and parsing mOTUpan-outs", logfile)

kegg_pangenome_stats = {}
os.makedirs(f"{kegg_folder}/mOTUs", exist_ok = True)
for motu, bins_ in tqdm(all_motus.items()):
    bins_ = set(bins_)
    if len(bins_) > min_genomes:
        with open(pjoin(temp_folder, "bin2ko.json") , "w") as handle:
            json.dump({k : list(v) for k,v in bin2ko.items() if k in bins_}, handle, indent = 2, sort_keys = True)
        call(f"mOTUpan.py --name {motu} --gene_clusters_file {temp_folder}/bin2ko.json  --checkm {temp_folder}/checkm.txt  -o {kegg_folder}/mOTUs/{motu}_ko.tsv --boots {boots}  >>{logfile} 2>&1 ", shell = True)
        with open(f"{kegg_folder}/mOTUs/{motu}_ko.tsv") as handle:
            header = {}
            main = []
            for l in handle:
                l = l.strip()
                if l != "#" and not l.startswith("#mOTUlizer"):
                    if l.startswith("#"):
                        header[l.split("=")[0][1:]] = "=".join(l.split("=")[1:])
                    else :
                        main += [l]
            cols = main[0].split()
            gcs = [{k : ll for k,ll in zip(cols, l.split())} for l in main[1:] if l != '']
            core = []
            accessory = []
            for d in gcs:
                if d['type'] == 'core':
                    core.append(d['trait_name'])
                else :
                    accessory.append(d['trait_name'])
                    for dat in header['genomes'].split(";"):
                        dat = dat.split(":")
                        bin_ = dat[0]
                        new_comp = [dd for dd in dat[1:] if dd.startswith("posterior_complete=")][0]
                        binset_stats[bin_]['mOTUpan_completeness'] = float(new_comp.replace("posterior_complete=",""))
                        kegg_pangenome_stats[motu] = {
                        'genome_count' : int(header['genome_count']),
                        'core_length': int(header['core_length']),
                        'mean_prior_completeness' : float(header['mean_prior_completeness']),
                        'mean_posterior_completeness' : float(header['mean_posterior_completeness']),
                        "nb_bootstraps" : 10 if "bootstrapped_mean_false_positive_rate" not in header else int(header['bootstrapped_nb_reps']),
                        "fpr" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[0]),
                        "sd_fpr" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[1].split("=")[1]),
                        "lowest_false" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[0]),
                        "sd_lowest_false" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[1].split("=")[1]),
                        'core' : ";".join(core),
                        'accessory' : ";".join(accessory)
                        }

module_path = pjoin(os.path.dirname(anvio.__file__), "data/misc/KEGG/modules/")
title2log("Loading kegg_paths",logfile)
module2name = {}
module2def = {}
module2ko = {}

def s(keggs, def_line, andline = True):
    def_line = def_line.split() if andline else def_line.split(",")
    i = 0
    blocks = []
    fwd = None
    rev = None
    for v in def_line:
        blocks += [i]
        if v.startswith("("):
            fwd = v.count("(")
            rev = v.count(")")
        elif fwd:
            fwd += v.count("(")
            rev += v.count(")")
        if fwd == rev:
            i += 1
            fwd = None
            rev = None
    operation = " " if andline else ","
    big_operation = lambda x: sum(x)/len(x) if andline else max(x)
    and_block = [operation.join([b for i,b in enumerate(def_line) if blocks[i] == block]) for block in set(blocks)]
    and_block = [b[1:-1] if b.startswith("(") else b for b in and_block]
    return big_operation([s(keggs, b, not andline) if ("," in b or " " in b) else int((b in keggs) if not b.startswith("-") else True) for b in and_block])


for m in os.listdir(module_path):
    with open(pjoin(module_path, m)) as handle:
        lines = handle.readlines()
        for l in lines:
            if l.startswith("DEFINITION"):
                module2def[m] = " ".join(l.split()[1:])
            if l.startswith("NAME"):
                module2name[m] = " ".join(l.split()[1:])

motu2modules = {}
for motu, dd in tqdm(kegg_pangenome_stats.items()):
    kos = set(dd['core'].split(";"))
    completes = { module : s(kos, defline)for module, defline in module2def.items()}
    completes = { k : v for k, v in completes.items() if v >= module_cutoff}
    motu2modules[motu] = completes


motu2core  = {k :set(v)  for k,v in motu2modules.items()}
union_sizes = {}
inter_sizes = {}
for g1, gc1 in tqdm(motu2core.items()):
    for g2, gc2 in motu2core.items():
        pair = frozenset((g1,g2))
        if pair not in union_sizes:
            inti = len(gc1.intersection(gc2))
            inter_sizes[pair] = inti
            union_sizes[pair] = len(gc1) + len(gc2) - inti

union_sizes = {k : v for k,v in union_sizes.items() if len(k) >1}
inter_sizes = {k : v for k,v in inter_sizes.items() if len(k) >1}

jaccard_index = {motu : {} for motu in kegg_pangenome_stats}
inter_by_small_size = {motu : {} for motu in kegg_pangenome_stats}

for k, inter in inter_sizes.items():
    motu1, motu2 = tuple(k)
    if motu1 not in jaccard_index:
        jaccard_index[motu1] = {}
        inter_by_small_size[motu1] = {}
    if motu2 not in jaccard_index:
        jaccard_index[motu2] = {}
        inter_by_small_size[motu2] = {}

    if union_sizes[k] == 0 :
        jaccard_index[motu1][motu2] = 0
        jaccard_index[motu2][motu1] = 0
    else :
        jaccard_index[motu2][motu1] = inter/union_sizes[k]
        jaccard_index[motu1][motu2] = inter/union_sizes[k]
        inter_by_small_size[motu1][motu2] =  0 if not min(len(motu2modules[motu1]),len(motu2modules[motu2])) else inter/min(len(motu2modules[motu1]),len(motu2modules[motu2]))
        inter_by_small_size[motu2][motu1] = 0 if not min(len(motu2modules[motu1]),len(motu2modules[motu2])) else inter/min(len(motu2modules[motu1]),len(motu2modules[motu2]))

motu_kegg_jaccard = pandas.DataFrame.from_records(jaccard_index).fillna(1)
motu_kegg_inter_by_small = pandas.DataFrame.from_records(inter_by_small_size).fillna(1)
orderz = sorted(motu_kegg_inter_by_small.columns)

motu_kegg_inter_by_small.loc[orderz,orderz].to_csv(pjoin(temp_folder, "motu_kegg_inter_by_small.csv"))
motu_kegg_jaccard.loc[orderz,orderz].to_csv(pjoin(temp_folder, "motu_kegg_jaccard.csv"))
motu_stats = {v['mOTU'] : v for k,v in bindat.items() if k == v['representative']}
for k in motu_stats:
    if k in gc_pangenome_stats:
        keys = ['core_length','nb_bootstraps', 'fpr', 'sd_fpr', 'lowest_false', 'sd_lowest_false', 'core', 'accessory']
        gc_stats = { "gc_" + kk : v for kk, v in gc_pangenome_stats[k].items() if kk in keys}
        kegg_stats = { "kegg_" + kk : v for kk, v in kegg_pangenome_stats[k].items() if kk in keys}
        motu_stats[k].update(gc_stats)
        motu_stats[k].update(kegg_stats)
        motu_stats[k]['ko_modules'] = ";".join(motu2modules[k])

pandas.DataFrame.from_dict(motu_stats, orient = "index").to_csv(pjoin(temp_folder, "motu_dat.csv"), index_label="mOTU")

title2log("Cleaning up and moving the bins", logfile)


shutil.rmtree(temp_folder)
