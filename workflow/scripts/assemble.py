import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import generate_config, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO

script , ass_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', ass_name + ".log")

config_file = generate_config(config_file)


call("conda env export > {out_folder}/logs/assembly.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/assembly_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


temp_folder = pjoin(config_file['temp_folder'], "assemblies", ass_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("Ungzipping and cat-ing libs to temp_folder".format(ass_name = ass_name), logfile)

for lib in config_file['assemblies'][ass_name]['libraries']:
    call("""
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_fwd.fastq.gz >> {temp}/fwd.fastq
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_rev.fastq.gz >> {temp}/rev.fastq
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_unp.fastq.gz >> {temp}/unp.fastq
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder), shell=True)

if config_file['assemblies'][ass_name]['preprocess'] == 'none':
    if config_file['assemblies'][ass_name]['keep_unpaired'] :
        unp = "-r  {temp}/unp.fastq "
    else :
        unp = ""
    megahit_line = "megahit -m 0.9 -1 {temp}/fwd.fastq -2 {temp}/rev.fastq " + unp + "-t {threads} -o {temp}/assembly --min-contig-len {min_len} 2> {log}"
else :
    print("Other preprocesssing then 'none' not implemented yet")
    system.exit(0)

title2log("assembling {ass_name}".format(ass_name = ass_name), logfile)

call(megahit_line.format(temp = temp_folder, threads = threads, min_len = config_file['assemblies'][ass_name]['length_cutoff'], log = "{out_folder}/logs/megahit.log".format(out_folder = out_folder)), shell = True)


title2log("Cleaning up and moving {ass_name}".format(ass_name = ass_name), logfile)

nb_contigs = len([ None for s in tqdm(SeqIO.parse(pjoin(temp_folder, "assembly", "final.contigs.fa"), "fasta"))])
zeros = len(str(nb_contigs))

max_buffer_size = config_file['seqio_buffer_size']
buffer = []
with open(pjoin(out_folder, "assembly.fna"), "w") as handle:
    for i,s in tqdm(enumerate(SeqIO.parse(pjoin(temp_folder, "assembly", "final.contigs.fa"), "fasta"))):
        s.id = ass_name + "_" + str(i+1).zfill(zeros)
        s.description = ""
        buffer += [s]
        if len(buffer) > max_buffer_size:
            SeqIO.write(buffer, handle, "fasta")
            buffer = []
    SeqIO.write(buffer, handle, "fasta")


shutil.rmtree(temp_folder)
