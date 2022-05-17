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
import numpy

script , mapping_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "mapping.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/mapping.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/mapping_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


subset = config_file['mappings'][mapping_name]['subset']
is_rna = config_file['mappings'][mapping_name]['is_rna']
method = config_file['mappings'][mapping_name]['mapper']
binset = config_file['mappings'][mapping_name]['binset']
taxfield = config_file['mappings'][mapping_name]['taxfield']
precluster = config_file['mappings'][mapping_name]['precluster']
keep_mapped = config_file['mappings'][mapping_name]['keep_mapped']

alternate_root = config_file['mappings'][mapping_name]['alternate_root']
ani = config_file['mappings'][mapping_name]['min_nucleotide_id']
min_len = config_file['mappings'][mapping_name]['min_len']

threads = int(threads)
mrna_flag = "_mrna" if is_rna else ""

if not alternate_root :
    alternate_root = binset

temp_folder = pjoin(config_file['temp_folder'], "mappings", mapping_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("copying binset to temp_folder", logfile)

shutil.copy(pjoin(root_folder, "binsets", binset, alternate_root + ".fna"), pjoin(temp_folder, "binset.fna"))
seqid = 0.95
cov = 0.9
covmode = 2
if precluster:
    call(f"mmseqs easy-cluster --min-seq-id {seqid} --cov-mode {covmode} -c {cov} --threads {threads} {temp_folder}/binset.fna {temp_folder}/binset {temp_folder}/mmseqs_temp #> {logfile}", shell=True)


title2log("indexing binset to temp_folder", logfile)
if method == "bwa-mem2":
    call("bwa-mem2 index {temp}/binset.fna >> {out_folder}/logs/mapping.log  2>&1".format(temp = temp_folder, threads = threads, out_folder = out_folder), shell=True)
if method == "bowtie2":
    call(f"bowtie2-build --threads {threads} {temp_folder}/binset.fna {temp_folder}/binset.fna", shell=True)
if method == "bbmap.sh":
    pass
if method == "minimap2":
       call(f"minimap2 -I 30G -t {threads} -d {temp_folder}/binset.idx {temp_folder}/binset.fna", shell=True)

freetxt_line("Starting mappings", logfile)
if os.path.exists(f"{temp_folder}/coverages.json"):
    with open(f"{temp_folder}/coverages.json") as handle:
        tt = json.load(handle)
    coverages = tt['coverages']
    total_reads = tt['total_reads']
    del tt
else :
    coverages = {}
    total_reads = {}

for lib in config_file['mappings'][mapping_name]['libraries']:
    if lib not in coverages:
        title2log("copying {subset}{lib}{mrna_flag} to temp_folder".format(lib = lib, subset = "subset_" if subset else "", mrna_flag = mrna_flag), logfile)
        call("""
        unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}{mrna_flag}_fwd.fastq.gz > {temp}/fwd.fastq 2>> {out_folder}/logs/mapping.log
        unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}{mrna_flag}_rev.fastq.gz > {temp}/rev.fastq 2>> {out_folder}/logs/mapping.log
        unpigz -kc {root_folder}/libraries/{lname}{subset_str}{lname}{mrna_flag}_unp.fastq.gz > {temp}/unp.fastq 2>> {out_folder}/logs/mapping.log
        """.format(root_folder = root_folder, lname = lib, temp=temp_folder, out_folder = out_folder, subset_str = "/subs/subs_" if subset else "/", mrna_flag = mrna_flag), shell=True)
        title2log("mapping {lib} to ref with {method}".format(lib = lib, method = method), logfile)
        if method == "bwa-mem2":
            call("""
            bwa-mem2 mem -t {threads} {temp}/binset.fna  -o {temp}/mapping.sam {temp}/fwd.fastq {temp}/rev.fastq 2>> {out_folder}/logs/mapping.log
            samtools view  -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_pairs.bam - >> {out_folder}/logs/mapping.log 2>&1
            bwa-mem2 mem -t {threads} {temp}/binset.fna  -o {temp}/mapping.sam {temp}/unp.fastq 2>> {out_folder}/logs/mapping.log
            samtools view -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_unpaired.bam - >> {out_folder}/logs/mapping.log 2>&1
            rm {temp}/mapping.sam 2>> {out_folder}/logs/mapping.log
            samtools merge -f -t {threads} {temp}/mapping.bam  {temp}/mapping_pairs.bam  {temp}/mapping_unpaired.bam >> {out_folder}/logs/mapping.log  2>&1
            rm {temp}/mapping_pairs.bam {temp}/mapping_unpaired.bam 2>> {out_folder}/logs/mapping.log
            """.format(root_folder = root_folder, lname = lib, temp=temp_folder, threads = threads, out_folder = out_folder), shell=True)
        if method == "bowtie2" :
            call("""
            bowtie2 -p {threads} -x  {temp}/binset.fna --very-sensitive -S {temp}/mapping.sam  -1 {temp}/fwd.fastq -2 {temp}/rev.fastq >> {out_folder}/logs/mapping.log 2>&1
            samtools view  -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_pairs.bam - >> {out_folder}/logs/mapping.log 2>&1
            bowtie2 -p {threads} -x  {temp}/binset.fna --very-sensitive -S {temp}/mapping.sam  -U {temp}/unp.fastq >> {out_folder}/logs/mapping.log 2>&1
            samtools view -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/mapping_unpaired.bam - >> {out_folder}/logs/mapping.log 2>&1
            rm {temp}/mapping.sam 2>> {out_folder}/logs/mapping.log
            samtools merge -f -t {threads} {temp}/mapping.bam  {temp}/mapping_pairs.bam  {temp}/mapping_unpaired.bam >> {out_folder}/logs/mapping.log  2>&1
            rm {temp}/mapping_pairs.bam {temp}/mapping_unpaired.bam 2>> {out_folder}/logs/mapping.log
            """.format(root_folder = root_folder, lname = lib, temp=temp_folder, threads = threads, out_folder = out_folder), shell=True)
        if method == "minimap2":
            call(f"""
            minimap2 -x sr --secondary=no -t 24  {temp_folder}/binset.idx {temp_folder}/fwd.fastq {temp_folder}/rev.fastq -a -o {temp_folder}/mapping.sam --MD  >> {out_folder}/logs/mapping.log 2>&1
            samtools view  --reference {temp_folder}/binset.fna -F0x900 -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {out_folder}/logs/mapping.log 2>&1
            minimap2 -x sr --secondary=no -t 24  {temp_folder}/binset.idx {temp_folder}/unp.fastq -a -o {temp_folder}/mapping.sam --MD  >> {out_folder}/logs/mapping.log 2>&1
            samtools view  --reference {temp_folder}/binset.fna -F0x900 -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_unpaired.bam - >> {out_folder}/logs/mapping.log 2>&1
            rm {temp_folder}/mapping.sam 2>> {out_folder}/logs/mapping.log
            samtools merge -f -t {threads} {temp_folder}/mapping.bam  {temp_folder}/mapping_pairs.bam  {temp_folder}/mapping_unpaired.bam >> {out_folder}/logs/mapping.log  2>&1
            rm {temp_folder}/mapping_pairs.bam {temp_folder}/mapping_unpaired.bam 2>> {out_folder}/logs/mapping.log
            """, shell = True)


        call("""
        coverm filter -b {temp}/mapping.bam -o {temp}/mapping_filtered.bam --min-read-percent-identity {ani} --min-read-aligned-length {min_len} --threads {threads} >> {out_folder}/logs/mapping.log  2>&1
        rm {temp}/mapping.bam >> {out_folder}/logs/mapping.log  2>&1
        coverm contig  --bam-files {temp}/mapping_filtered.bam  --methods count  --threads {threads} > {temp}/mapping.tsv 2>> {out_folder}/logs/mapping.log
        """.format(temp=temp_folder, threads = threads, ani = ani,  out_folder = out_folder,  min_len = min_len), shell=True)
        if keep_mapped:
            title2log("extracting mapped reads from  {lib}".format(lib = lib), logfile)
            call(f"samtools fastq -@ {threads} {temp_folder}/mapping_filtered.bam -o {temp_folder}/{lib}.fastq  2>> {out_folder}/logs/mapping.log", shell=True)
        call("""
        rm {temp}/mapping_filtered.bam
        cat {temp}/fwd.fastq {temp}/rev.fastq {temp}/unp.fastq | wc -l > {temp}/total_reads_x4.txt
        """.format(temp=temp_folder, threads = threads, ani = ani,  out_folder = out_folder,  min_len = min_len), shell=True)

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

        with open(f"{temp_folder}/coverages.json", "w") as handle:
            json.dump({"coverages" : coverages, "total_reads" : total_reads}, handle)



title2log("Done with the mappings", logfile)

title2log("Making tables", logfile)

if not is_rna :
    for k in total_reads:
        coverages[k]['unmapped'] = total_reads[k] - sum(coverages[k].values())

with open(pjoin(temp_folder, "total_reads_to_map.csv"), "w") as handle:
    handle.writelines(["library,total_reads"] + [f"{k},{v}\n" for k,v in total_reads.items()])

coverages = pandas.DataFrame.from_dict(coverages)
coverages.to_csv(pjoin(temp_folder, "contigs_mapped_reads.csv"), index_label = "contig_name")
normed_cov = coverages/coverages.sum()
normed_cov.to_csv(pjoin(temp_folder, "contigs_relative_abundance.csv"), index_label = "contig_name")

if is_rna :
    gene_stats  = pandas.read_csv(pjoin(root_folder, "binsets", binset, alternate_root + "_basics.csv" ), index_col=0)
    lens = [len(gene_stats.loc[g, 'representative_nucls']) for g in normed_cov.index]
    rpk = 10_000*(coverages.transpose()/lens).transpose()
    norm_facts = rpk.sum()/1_000_000
    tpm = rpk/norm_facts
    tpm.to_csv(pjoin(temp_folder, "tpm.csv"), index_label = "derep_gene")

    for l in ['root_eggNOG', 'COG_category', 'symbol', 'KO']:
        relabs = tpm.copy()
        relabs[l] = [t for t in gene_stats.loc[tpm.index,l] ]
        if l in ['KO', 'COG_category']:
            if l == "KO":
                relabs[l] = [ "" if v != v else v.split(",") for v in relabs[l]]
            else :
                relabs[l] = [ "" if v != v else  list(v) for v in relabs[l]]
            multiples = [i for i,n in zip(relabs.index,relabs[l]) if n == n and len(n) >1 ]
            lines = []
            for gc in tqdm(multiples):
                line = relabs.loc[gc]
                values = line[l]
                del line[l]
                line = list(line/len(values))
                lines += [ line + [values[j]] for j in range(len(values))]
            relabs = relabs.loc[[i for i,n in zip(relabs.index,relabs[l]) if len(n) == 1 or n == ""]]
            relabs[l] = ["NA" if len(v) == 0 else v[0] for v in relabs[l]]
            lines =  pandas.DataFrame(lines)
            lines.columns = relabs.columns
            relabs = pandas.concat([relabs,lines])

        relabs = relabs.groupby(l).sum()
        relabs.to_csv(pjoin(temp_folder, l + "_relative_abundance.csv"), index_label = l)


else :
    coverages['bin'] = ["_".join(c.split("_")[:-1]) if c != "unmapped"  else "unmapped" for c in coverages.index]
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
call(f"mkdir -p {root_folder}/mappings/{mapping_name}/mapped_reads/; cp {temp_folder}/*.fastq {root_folder}/mappings/{mapping_name}/mapped_reads/", shell=True)

title2log("Cleaning up and moving the bins", logfile)


#shutil.rmtree(temp_folder)
