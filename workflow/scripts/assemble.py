import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import validate_description_json, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json

script , ass_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', ass_name + ".log")

config_file = validate_description_json(config_file)


call("conda env export > {out_folder}/logs/assembly.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/assembly.json".format(out_folder = out_folder), "w") as handle:
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
    megahit_line = "megahit -1 {temp}/fwd.fastq -2 {temp}/rev.fastq -r  {temp}/unp.fastq -t {threads} -o {temp}/assembly --min-contig-len {min_len} > {log}.log"
else :
    print("Other preprocesssing then 'none' not implemented yet")
    system.exit(0)

title2log("assembling {ass_name}".format(ass_name = ass_name), logfile)

call(megahit_line.format(temp = temp_folder, threads = threads, min_len = config_file['assemblies'][ass_name]['length_cutoff'], log = "{out_folder}/logs/megahit.log".format(out_folder = out_folder)), shell = True)



title2log("Cleaning up and moving".format(ass_name = ass_name), logfile)

shutil.move(pjoin(temp_folder, "smrna_paired", "out", "aligned_fwd.fastq") , pjoin(temp_folder, "rrna_fwd.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","aligned_rev.fastq") , pjoin(temp_folder, "rrna_rev.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","other_fwd.fastq") , pjoin(temp_folder, "mrna_fwd.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","other_rev.fastq") , pjoin(temp_folder, "mrna_rev.fastq") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","aligned.fastq") , pjoin(temp_folder, "rrna_unp.fastq") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","other.fastq") , pjoin(temp_folder, "mrna_unp.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","aligned.log") , pjoin(out_folder, "stats", "sortmerna_" + ass_name + "_paired.stats") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","aligned.log") , pjoin(out_folder, "stats", "sortmerna_" + ass_name + "_unpaired.stats") )



to_gz = ['{temp}/rrna_fwd.fastq',
         '{temp}/mrna_fwd.fastq',
         '{temp}/rrna_rev.fastq',
         '{temp}/mrna_rev.fastq',
         '{temp}/rrna_unp.fastq',
         '{temp}/mrna_unp.fastq',
         ]

to_gz = [g.format(temp = temp_folder, lname = ass_name) for g in to_gz]

for f in to_gz:
    call("pigz {}".format(f), shell=True)

to_rename = ['rrna_fwd.fastq.gz',
         'mrna_fwd.fastq.gz',
         'rrna_rev.fastq.gz',
         'mrna_rev.fastq.gz',
         'rrna_unp.fastq.gz',
         'mrna_unp.fastq.gz',
         ]

for f in to_rename:
    shutil.move(pjoin(temp_folder, f), pjoin(temp_folder, ass_name + "_" + f))

to_move = ['{temp}/{lname}_rrna_fwd.fastq.gz',
         '{temp}/{lname}_mrna_fwd.fastq.gz',
         '{temp}/{lname}_rrna_rev.fastq.gz',
         '{temp}/{lname}_mrna_rev.fastq.gz',
         '{temp}/{lname}_rrna_unp.fastq.gz',
         '{temp}/{lname}_mrna_unp.fastq.gz',
         ]

to_move = [g.format(temp = temp_folder, lname = ass_name) for g in to_move]

for f in to_move:
    shutil.move(f, out_folder )

shutil.rmtree(temp_folder)
