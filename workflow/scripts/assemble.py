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

script , ass_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', ass_name + ".log")

config_file = generate_config(config_file)


call("conda env export > {out_folder}/logs/assembly.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/assembly_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


temp_folder = pjoin(config_file['temp_folder'], "assemblies", ass_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("cat-ing libs to temp_folder".format(ass_name = ass_name), logfile)

for lib in config_file['assemblies'][ass_name]['libraries']:
    call("""
    cat {root_folder}/libraries/{lname}/{lname}_fwd.fastq.gz >> {temp}/fwd.fastq.gz
    cat {root_folder}/libraries/{lname}/{lname}_rev.fastq.gz >> {temp}/rev.fastq.gz
    cat {root_folder}/libraries/{lname}/{lname}_unp.fastq.gz >> {temp}/unp.fastq.gz
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder), shell=True)



if config_file['assemblies'][ass_name]['preprocess'] == 'none':
    pass
elif config_file['assemblies'][ass_name]['preprocess'] == 'bbnorm':
    title2log("Running bbnorm diginorm".format(ass_name = ass_name), logfile)

    call("""
    bbnorm.sh -Xmx100g in={temp}/fwd.fastq.gz in2={temp}/rev.fastq.gz out={temp}/tt_fwd.fastq.gz out2={temp}/tt_rev.fastq.gz t={threads} 2>> {log_file}
    bbnorm.sh -Xmx100g in={temp}/unp.fastq.gz out={temp}/tt_unp.fastq.gz t={threads} 2>> {log_file}
    mv {temp}/tt_fwd.fastq.gz {temp}/fwd.fastq.gz
    mv {temp}/tt_rev.fastq.gz {temp}/rev.fastq.gz
    mv {temp}/tt_unp.fastq.gz {temp}/unp.fastq.gz
    """.format(temp = temp_folder, threads = threads, log_file = logfile), shell = True)
else :
    print("Other preprocesssing then 'none' or 'bbnorm' not implemented yet")
    system.exit(0)

if config_file['assemblies'][ass_name]['assembler'] == 'megahit':
    if config_file['assemblies'][ass_name]['keep_unpaired'] :
        unp = "-r  {temp}/unp.fastq.gz "
    else :
        unp = ""
    title2log("Running megahit".format(ass_name = ass_name), logfile)
    megahit_line = "megahit -m 0.9 -1 {temp}/fwd.fastq.gz -2 {temp}/rev.fastq.gz " + unp + "-t {threads} -o {temp}/assembly --min-contig-len {min_len} 2> {log}"


    title2log("assembling {ass_name}".format(ass_name = ass_name), logfile)
    call(megahit_line.format(temp = temp_folder, threads = threads, min_len = config_file['assemblies'][ass_name]['length_cutoff'], log = "{out_folder}/logs/megahit.log".format(out_folder = out_folder)), shell = True)

    title2log("Cleaning up and moving {ass_name}".format(ass_name = ass_name), logfile)
    nb_contigs = len([ None for s in tqdm(SeqIO.parse(pjoin(temp_folder, "assembly", "final.contigs.fa"), "fasta"))])
    max_buffer_size = config_file['seqio_buffer_size']
    zeros = len(str(nb_contigs))


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
elif config_file['assemblies'][ass_name]['assembler'] == 'spades':
        if config_file['assemblies'][ass_name]['keep_unpaired'] :
            unp = "-s  {temp}/unp.fastq.gz "
        else :
            unp = ""
        title2log("Running spades".format(ass_name = ass_name), logfile)
        megahit_line = "spades.py --meta -1 {temp}/fwd.fastq.gz -2 {temp}/rev.fastq.gz " + unp + "-t {threads} -o {temp}/assembly > {log}"


        title2log("assembling {ass_name}".format(ass_name = ass_name), logfile)
        call(megahit_line.format(temp = temp_folder, threads = threads, log = "{out_folder}/logs/spades.log".format(out_folder = out_folder)), shell = True)

        title2log("Cleaning up and moving {ass_name}".format(ass_name = ass_name), logfile)
        nb_contigs = len([ None for s in tqdm(SeqIO.parse(pjoin(temp_folder, "assembly", "scaffolds.fasta"), "fasta"))])
        zeros = len(str(nb_contigs))

        max_buffer_size = config_file['seqio_buffer_size']
        buffer = []
        shutil.copy(pjoin(temp_folder, "assembly",'assembly_graph.fastg'), pjoin(out_folder, 'assembly.fastg') )
        with open(pjoin(out_folder, "assembly.fna"), "w") as handle:
            for i,s in tqdm(enumerate(SeqIO.parse(pjoin(temp_folder, "assembly", "scaffolds.fasta"), "fasta"))):
                if len(s) > config_file['assemblies'][ass_name]['length_cutoff']:
                    s.id = ass_name + "_" + str(i+1).zfill(zeros)
                    s.description = ""
                    buffer += [s]
                    if len(buffer) > max_buffer_size:
                        SeqIO.write(buffer, handle, "fasta")
                        buffer = []
            SeqIO.write(buffer, handle, "fasta")
else :
    print("Other assembler than 'megahit' and 'spades' not implemented yet")
    system.exit(0)


shutil.rmtree(temp_folder)
