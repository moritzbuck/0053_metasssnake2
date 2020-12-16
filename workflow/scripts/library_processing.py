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

call("conda env export > {out_folder}/logs/library_processing.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/library_processing_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


title2log("Starting processing library {}".format(lib_name), logfile)

temp_folder = pjoin(config_file['temp_folder'], "library_processing", lib_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)
os.makedirs(pjoin(out_folder,  "logs/fastp_logs/"))

paired_fastp_line = "fastp -h /dev/null -j {temp}/{lib1}.json  --in1 {temp}/{lib1} --in2 {temp}/{lib2} --out1 {temp}/{lib1}_clean.fastq --out2 {temp}/{lib2}_clean.fastq --unpaired1 {temp}/{lib1}_unp.fastq --unpaired2 {temp}/{lib2}_unp.fastq  -w {threads}  >> {log} 2>&1"
single_fastp_line = "fastp -h /dev/null -j {temp}/{lib}.json  --in1 {temp}/{lib}  --out1 {temp}/{lib}_clean.fastq  -w {threads}  >> {log} 2>&1"
qc_log = {'paired' : dict(), 'unpaired' : dict()}

for fwd, rev in zip(config_file['libraries'][lib_name]["fwd"], config_file['libraries'][lib_name]["rev"]):
    title2log("QCing paired reads_library {} and {}, ".format(os.path.basename(fwd), os.path.basename(rev)), logfile)
    freetxt_line("first copying reads to {temp_folder} ".format(temp_folder = temp_folder), logfile)
    shutil.copy(pjoin(config_file['root_folder'], fwd), temp_folder)
    shutil.copy(pjoin(config_file['root_folder'], rev), temp_folder)

    freetxt_line("Doing actual QC", logfile)
    call(paired_fastp_line.format(lib1 = basename(fwd), lib2 = basename(rev), out1 = basename(fwd) + "clean.fastq", out2 = basename(rev) + "clean.fastq",
        threads = threads,  temp=temp_folder, log = pjoin(out_folder,  "logs/fastp_logs/{}.log".format(basename(fwd)))
    ), shell = True)
    freetxt_line("Concatenanting stats and reads", logfile)
    with open(pjoin(temp_folder,  basename(fwd) + ".json")) as handle:
        qc_log['paired'][fwd] = json.load(handle)
    shutil.move(pjoin(temp_folder,  basename(fwd) + ".json"), pjoin(out_folder,  "logs/fastp_logs/"))
    cat_line = """
    cat {temp}/{lib1}_clean.fastq >> {temp}/fwd.fastq 2>> {log}
    cat {temp}/{lib2}_clean.fastq >> {temp}/rev.fastq 2>> {log}
    cat {temp}/{lib1}_unp.fastq >> {temp}/unp.fastq 2>> {log}
    cat {temp}/{lib2}_unp.fastq >> {temp}/unp.fastq 2>> {log}
    """
    call(cat_line.format(lib1 = basename(fwd),lib2 = basename(rev), temp=temp_folder, log = logfile), shell = True)
    freetxt_line("Cleaning up", logfile)
    rm_line = """
    rm {temp}/{lib1} 2>> {log}
    rm {temp}/{lib2} 2>> {log}
    rm {temp}/{lib1}_clean.fastq 2>> {log}
    rm {temp}/{lib2}_clean.fastq 2>> {log}
    rm {temp}/{lib1}_unp.fastq 2>> {log}
    rm {temp}/{lib2}_unp.fastq 2>> {log}
    rm {temp}/{lib1}.json 2>> {log}
    """
    call(rm_line.format(lib1 = basename(fwd),lib2 = basename(rev), temp=temp_folder, log = logfile), shell = True)


for unp in config_file['libraries'][lib_name]["unp"]:
    title2log("QCing unpaired reads_library {}".format(os.path.basename(unp)), logfile)
    freetxt_line("first copying reads to {temp_folder} ".format(temp_folder = temp_folder), logfile)
    shutil.copy(pjoin(config_file['root_folder'], unp), temp_folder)

    freetxt_line("Doing actual QC", logfile)
    call(single_fastp_line.format(lib = basename(unp), out1 = basename(unp) + "clean.fastq",
        threads = threads,  temp=temp_folder, log = logfile
    ), shell = True)
    freetxt_line("Concatenanting stats and reads", logfile)
    with open(pjoin(temp_folder,  basename(unp) + ".json")) as handle:
        qc_log['unpaired'][unp] = json.load(handle)
    shutil.move(pjoin(temp_folder,  basename(unp) + ".json"), pjoin(out_folder,  "logs/fastp_logs/"))
    cat_line = """
    cat {temp}/{lib}_unp.fastq >> {temp}/unp.fastq 2>> {log}
    """
    call(cat_line.format(lib = basename(unp), temp=temp_folder, log = logfile), shell = True)
    freetxt_line("Cleaning up", logfile)
    rm_line = """
    rm {temp}/{lib} 2>> {log}
    rm {temp}/{lib}_clean.fastq 2>> {log}
    rm {temp}/{lib}.json 2>> {log}
    """
    call(rm_line.format(lib = basename(unp), temp=temp_folder, log = logfile), shell = True)


with open(pjoin(out_folder, "stats/fastp_{lib_name}.stats".format(lib_name = lib_name)), "w") as handle:
   json.dump(qc_log, handle, indent = 2, sort_keys = True)

title2log("Read subsetting", logfile)


if not os.path.exists(pjoin(temp_folder, "fwd.fastq")):
    only_singles = True
    paired_size = 0
else :
    only_singles = False
    paired_size = os.path.getsize(pjoin(temp_folder, "fwd.fastq"))+os.path.getsize(pjoin(temp_folder, "rev.fastq"))

unpaired_size = os.path.getsize(pjoin(temp_folder, "unp.fastq"))

read_proc_line = """
reformat.sh {in_reads} {out_reads} samplereadstarget={subcount} sampleseed=42 t={threads}  2>> {log}
"""

if only_singles:
    inreads = "in={temp}/unp.fastq".format(temp = temp_folder)
    outreads = "out={temp}/subs/subs_{lname}_unp.fastq".format(temp = temp_folder, lname = lib_name)
    call(read_proc_line.format(in_reads = inreads, out_reads = outreads, subcount = config_file['libraries_config']['read_subset'], threads = threads, log = pjoin(out_folder,  "logs/bbtools_reformat.log")), shell=True)
    with open("{temp}/subs/subs_{lname}_fwd.fastq".format(temp = temp_folder, lname = lib_name), "w") as handle:
        pass
    with open("{temp}/subs/subs_{lname}_rev.fastq".format(temp = temp_folder, lname = lib_name), "w") as handle:
        pass
else :
    unp_ratio = unpaired_size/(unpaired_size+paired_size)
    inreads = "in={temp}/unp.fastq".format(temp = temp_folder)
    outreads = "out={temp}/subs/subs_{lname}_unp.fastq".format(temp = temp_folder, lname = lib_name)
    call(read_proc_line.format(in_reads = inreads, out_reads = outreads, subcount = int(config_file['libraries_config']['read_subset']*unp_ratio), threads = threads, log = pjoin(out_folder,  "logs/bbtools_reformat.log")), shell=True)

    inreads = "in={temp}/fwd.fastq in2={temp}/rev.fastq ".format(temp = temp_folder)
    outreads = "out={temp}/subs/subs_{lname}_fwd.fastq out2={temp}/subs/subs_{lname}_rev.fastq ".format(temp = temp_folder, lname = lib_name )
    call(read_proc_line.format(in_reads = inreads, out_reads = outreads, subcount = int(config_file['libraries_config']['read_subset']*(1-unp_ratio)), threads = threads, log = pjoin(out_folder,  "logs/bbtools_reformat.log")), shell=True)

title2log("Read sketching", logfile)

sourmash_line = """
sourmash compute --track-abundance --merge {lname} -k{k} --scaled {scale} {temp}/fwd.fastq {temp}/rev.fastq {temp}/unp.fastq -o {temp}/{lname}.sig -p {threads}  2>> {log}
"""

call(sourmash_line.format(k = config_file['libraries_config']['sourmash_k'], scale = config_file['libraries_config']['sourmash_scaled'], lname = lib_name, temp = temp_folder, threads = threads, log = pjoin(out_folder,  "logs/sourmash_compute.log")), shell=True)
with open("{temp}/{lname}.sig".format(temp = temp_folder, lname = lib_name) , "w") as handle:
    pass

title2log("zipping things and moving back", logfile)

to_gz = ['{temp}/fwd.fastq',
         '{temp}/rev.fastq',
         '{temp}/unp.fastq',
         '{temp}/{lname}.sig',
         '{temp}/subs/subs_{lname}_unp.fastq',
         '{temp}/subs/subs_{lname}_fwd.fastq',
         '{temp}/subs/subs_{lname}_rev.fastq',
         ]
to_gz = [g.format(temp = temp_folder, lname = lib_name) for g in to_gz]

for f in to_gz:
    call("pigz {}".format(f), shell=True)

to_rename = ['fwd.fastq.gz',
         'rev.fastq.gz',
         'unp.fastq.gz']

for f in to_rename:
    shutil.move(pjoin(temp_folder, f), pjoin(temp_folder, lib_name + "_" + f))

to_move = ['{temp}/{lname}_fwd.fastq.gz',
         '{temp}/{lname}_rev.fastq.gz',
         '{temp}/{lname}_unp.fastq.gz',
         '{temp}/{lname}.sig.gz'
         ]

to_move = [g.format(temp = temp_folder, lname = lib_name) for g in to_move]

for f in to_move:
    shutil.move(f, out_folder )

for f in os.listdir(pjoin(temp_folder, "subs")):
    shutil.move(pjoin(temp_folder, "subs", f), pjoin(out_folder, "subs", f))

shutil.rmtree(temp_folder)
