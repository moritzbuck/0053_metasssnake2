import sys, os
from tqdm import tqdm
from os.path import join as pjoin
from workflow.scripts.utils import validate_description_json, title2log, freetxt_line
import shutil
from subprocess import call
from os.path import basename
import json

script , lib_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', lib_name + ".log")

config_file = validate_description_json(config_file)


call("conda env export > {out_folder}/logs/library_rrna_spliting.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/library_rrna_spliting_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


temp_folder = pjoin(config_file['temp_folder'], "library_processing", lib_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)



os.makedirs(temp_folder, exist_ok=True)

title2log("Ungzipping {lib_name}'s reads to temp_folder".format(lib_name = lib_name), logfile)

call("""
unpigz -kc {out_folder}/{lname}_fwd.fastq.gz > {temp}/fwd.fastq
unpigz -kc {out_folder}/{lname}_rev.fastq.gz > {temp}/rev.fastq
unpigz -kc {out_folder}/{lname}_unp.fastq.gz > {temp}/unp.fastq
""".format(out_folder = out_folder, lname = lib_name, temp=temp_folder), shell=True)

refs = "".join([" --ref " + f for f in config_file['libraries_config']['sortmerna_refs']])

title2log("Running sortmeRNA on {lib_name}".format(lib_name = lib_name), logfile)


call("""
sortmerna --task 3 --out2 --threads {threads} {refs}  --reads {temp}/fwd.fastq --reads {temp}/rev.fastq --workdir {temp}/smrna_paired/  -out2 -num_alignments 1 -v -fastx --threads {threads} --other > {log} 2>&1
echo "=========== unpaired reads =============" >> {log} 2>&1
sortmerna --task 3  --threads {threads} {refs}  --reads {temp}/unp.fastq --workdir {temp}/smrna_unpaired/  -num_alignments 1 -v -fastx --threads {threads} --other >> {log} 2>&1
""".format(threads = threads, refs = refs, temp=temp_folder, log = pjoin(out_folder, "logs", "sortmerna.log")), shell=True)

title2log("Cleaning up and moving".format(lib_name = lib_name), logfile)

shutil.move(pjoin(temp_folder, "smrna_paired", "out", "aligned_fwd.fastq") , pjoin(temp_folder, "rrna_fwd.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","aligned_rev.fastq") , pjoin(temp_folder, "rrna_rev.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","other_fwd.fastq") , pjoin(temp_folder, "mrna_fwd.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","other_rev.fastq") , pjoin(temp_folder, "mrna_rev.fastq") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","aligned.fastq") , pjoin(temp_folder, "rrna_unp.fastq") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","other.fastq") , pjoin(temp_folder, "mrna_unp.fastq") )
shutil.move(pjoin(temp_folder, "smrna_paired","out","aligned.log") , pjoin(out_folder, "stats", "sortmerna_" + lib_name + "_paired.stats") )
shutil.move(pjoin(temp_folder, "smrna_unpaired","out","aligned.log") , pjoin(out_folder, "stats", "sortmerna_" + lib_name + "_unpaired.stats") )



to_gz = ['{temp}/rrna_fwd.fastq',
         '{temp}/mrna_fwd.fastq',
         '{temp}/rrna_rev.fastq',
         '{temp}/mrna_rev.fastq',
         '{temp}/rrna_unp.fastq',
         '{temp}/mrna_unp.fastq',
         ]

to_gz = [g.format(temp = temp_folder, lname = lib_name) for g in to_gz]

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
    shutil.move(pjoin(temp_folder, f), pjoin(temp_folder, lib_name + "_" + f))

to_move = ['{temp}/{lname}_rrna_fwd.fastq.gz',
         '{temp}/{lname}_mrna_fwd.fastq.gz',
         '{temp}/{lname}_rrna_rev.fastq.gz',
         '{temp}/{lname}_mrna_rev.fastq.gz',
         '{temp}/{lname}_rrna_unp.fastq.gz',
         '{temp}/{lname}_mrna_unp.fastq.gz',
         ]

to_move = [g.format(temp = temp_folder, lname = lib_name) for g in to_move]

for f in to_move:
    shutil.move(f, out_folder )

shutil.rmtree(temp_folder)
