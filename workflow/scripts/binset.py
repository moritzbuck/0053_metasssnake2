import sys, os
from tqdm import tqdm
from os.path import join as pjoin
sys.path.append(os.getcwd())

from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file, gff2anvio
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO
import re
from anvio.summarizer import ContigSummarizer
from tempfile import NamedTemporaryFile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from anvio.dbops import ContigsSuperclass

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

title2log("Parsing and filtering based on faa/fnas", logfile)

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
title2log("Live hacking prokka script to allow for long contig names", logfile)


call("sed -i 's/my $MAXCONTIGIDLEN = 37/my $MAXCONTIGIDLEN = 250/' `which prokka` ", shell=True)
call("sed -i 's/[^#]tbl2asn -V/#tbl2asn -V/' `which prokka` ", shell=True)

if keep_fails:
    nb_contigs = len([ None for s in tqdm(SeqIO.parse(pjoin(temp_folder, "clean_bins", binset_name + "_unkept.fna"), "fasta"))])
    zeros = len(str(nb_contigs))
    max_buffer_size = config_file['seqio_buffer_size']
    buffer = []
    with open(pjoin(temp_folder, "clean_bins", binset_name + "_unkept2.fna"), "w") as handle:
        for i,s in tqdm(enumerate(SeqIO.parse(pjoin(temp_folder, "clean_bins", binset_name + "_unkept.fna"), "fasta"))):
            s.id = binset_name + "-kept_" + str(i+1).zfill(zeros) + "__" + s.id
            s.description = ""
            buffer += [s]
            if len(buffer) > max_buffer_size:
                SeqIO.write(buffer, handle, "fasta")
                buffer = []
        SeqIO.write(buffer, handle, "fasta")

shutil.move(pjoin(temp_folder, "clean_bins", binset_name + "_unkept2.fna"),pjoin(temp_folder, "clean_bins", binset_name + "_unkept.fna"))
call("prokka --outdir {temp}/clean_bins/{binset}_unkept --prefix {binset}_unkept --locustag {binset}_unkept --metagenome --cpus {threads} {temp}/clean_bins/{binset}_unkept.fna >> {logfile}  2>&1".format(threads= threads, binset = binset_name, temp=temp_folder, logfile = logfile), shell=True)
os.remove(pjoin(cbinfoder, binset_name + "_unkept.fna"))

call("ls {temp}/bins/ | cut -f1 -d. | parallel -j{threads} prokka --outdir {temp}/clean_bins/{{}} --prefix {{}} --locustag {{}} --cpus 1 {temp}/bins/{{}}.fna >> {logfile}  2>&1".format(logfile = logfile, threads= threads, temp=temp_folder), shell = True)

title2log("Creating anvi'o databases", logfile)

formating_dat = {
'out_folder' : out_folder,
'temp_folder' : temp_folder,
'threads' : threads,
'binset_name' : binset_name
}
make_dbs_line = """ls {temp_folder}/clean_bins/ | grep _bin- | parallel -j{threads} anvi-gen-contigs-database --ignore-internal-stop-codons --quiet -n {binset_name} -f {temp_folder}/clean_bins/{{}}/{{}}.fna -o {temp_folder}/clean_bins/{{}}/{{}}.db -T 1 --external-gene-calls {temp_folder}/clean_bins/{{}}/{{}}.cdss
ls {temp_folder}/clean_bins/ | grep _bin- | parallel -j{threads} anvi-import-functions --quiet -c {temp_folder}/clean_bins/{{}}/{{}}.db -i {temp_folder}/clean_bins/{{}}/{{}}.annot"""
for bin_id in tqdm(os.listdir(pjoin(temp_folder, "clean_bins"))):
    if "_bin-" in bin_id:
        prokka_out = gff2anvio(pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".gff"))
        with open(pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".annot"), "w") as handle:
            handle.writelines(prokka_out['annot'])
        with open(pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".cdss"), "w") as handle:
            handle.writelines(prokka_out['cdss'])

call(make_dbs_line.format(**formating_dat), shell = True)
title2log("done with anvi'o databases", logfile)


anvi_pipe = ['anvi-run-hmms', 'anvi-run-kegg-kofams', 'anvi-run-scg-taxonomy']#, 'anvi-run-pfams']
anvi_line = "ls {temp_folder}/clean_bins/ | grep _bin- | parallel -j{threads} {program} --quiet --tmp-dir  {temp_folder} -c {temp_folder}/clean_bins/{{}}/{{}}.db -T {threads}"

title2log("Starting anvio pipe".format(**formating_dat), logfile)

for program in anvi_pipe:
    formating_dat['program'] = program
    title2log("running {program} on anvi'o databases".format(**formating_dat), logfile)
    call(anvi_line.format(**formating_dat), shell = True)

class mock:
     def __init__(self):
         self.__dict__ = {}
params = mock()

stats = {}
fields = ['project_name', 'contigs_db_hash', 'num_contigs', 'total_length', 'creation_date', 'gc_content', 'num_genes', 'percent_completion', 'percent_redundancy', 'scg_domain', 'scg_domain_confidence']
get_taxo_line = "anvi-estimate-scg-taxonomy --quiet -T {threads} -c {temp_folder}/clean_bins/{bin_id}/{bin_id}.db -o {tempfile}"
head = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
for bin_id in tqdm(os.listdir(pjoin(temp_folder, "clean_bins"))):
    if os.path.isdir(pjoin(temp_folder, "clean_bins", bin_id)):
        tt = ContigSummarizer(pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".db")).get_contigs_db_info_dict(gene_caller_to_use="Prodigal")
        t_file = NamedTemporaryFile()
        formating_dat['bin_id'] = bin_id
        formating_dat['tempfile'] = t_file.name
        call(get_taxo_line.format(**formating_dat), shell = True)
        with open(t_file.name) as handle:
            handle.readline().split()
            scg_taxo = handle.readline().strip().split("\t")
        params.__dict__['contigs_db'] = pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".db")
        c = ContigsSuperclass(params)
        calls = c.get_sequences_for_gene_callers_ids(simple_headers=False)[1]
        with open(pjoin(temp_folder, binset_name + ".faa"), "a") as handle:
            seqs = []
            for k,v in calls.items():
                ss = Seq(v['sequence']).translate()
                seqs.append(SeqRecord(ss, id = bin_id + ";" + str(k), description = ""))
            SeqIO.write(seqs, handle, "fasta")

        utils.export_sequences_from_contigs_db(pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".db"), t_file.name)
        call("cat {tempfile} >> {temp_folder}/{binset_name}.fna".format(**formating_dat), shell = True)
        t_file.close()
        est_coding = tt['avg_gene_length']*tt['num_genes']/tt['total_length']
        tt = {k : v for k,v in tt.items() if k in fields}
        stats[bin_id] = tt
        stats[bin_id]['scg_taxo_support'] = scg_taxo[1] + "/" + scg_taxo[2]
        stats[bin_id]['scg_taxo'] = ";".join([h+t for h,t in  zip(head,scg_taxo[3:])])
        stats[bin_id]['approx_coding_density'] = est_coding
        pjoin(temp_folder, "clean_bins", bin_id, bin_id + ".db")

dict2file(stats, pjoin(out_folder, binset_name + "_basics.csv"))

title2log("Cleaning up and moving the bins", logfile)

if os.path.exists(pjoin(out_folder, binset_name + ".fna")):
    os.remove(pjoin(out_folder,  binset_name +  ".fna"))

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".fna"), ass = pjoin(out_folder,  binset_name + ".fna")), shell=True)
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".faa"), ass = pjoin(out_folder,  binset_name + ".faa")), shell=True)
    call("sed '/##FASTA/q' {file} | grep -v '^# ' >> {ass}".format(file = pjoin(out_folder, "bins", file, file +".gff"), ass = pjoin(out_folder,  binset_name + ".gff")), shell=True)

title2log("Binsetting done", logfile)

shutil.rmtree(temp_folder)
