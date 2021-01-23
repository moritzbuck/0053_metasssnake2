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

script , binning_name, config_file, root_folder, out_folder, threads = sys.argv

logfile = pjoin(out_folder, 'logs', "binning.log")
config_file = generate_config(config_file)

call("conda env export > {out_folder}/logs/binning.yaml".format(out_folder = out_folder), shell=True)
with open("{out_folder}/logs/binning_settings.json".format(out_folder = out_folder), "w") as handle:
    json.dump(config_file, handle, indent = 2, sort_keys = True)


def clean_bin(inbin, outbin, name, buffer_size = 1000 ):
    lens = [ len(s) for s in SeqIO.parse(inbin, "fasta")]
    nb_contigs = len(lens)
    if sum(lens) > min_bin_size:
        zeros = len(str(nb_contigs))

        max_buffer_size = buffer_size
        buffer = []
        with open(outbin, "w") as handle:
            for i,s in enumerate(SeqIO.parse(inbin, "fasta")):
                s.id = name + "_ctg-" + str(i+1).zfill(zeros)
                s.description = ""
                buffer += [s]
                if len(buffer) > max_buffer_size:
                    SeqIO.write(buffer, handle, "fasta")
                    buffer = []
            SeqIO.write(buffer, handle, "fasta")
        return None
    else :
        return (inbin, len(lens))

min_bin_size = config_file['binnings'][binning_name]['min_bin_size']
method = config_file['binnings'][binning_name]['binner']

cuda = config_file['binnings'][binning_name]['other_parameters'].get('cuda')
if cuda :
    cuda = bool(cuda)

temp_folder = pjoin(config_file['temp_folder'], "binnings", binning_name)
freetxt_line("Creating temp folder: " + temp_folder, logfile)

os.makedirs(temp_folder, exist_ok=True)

title2log("copying assemblies to temp_folder", logfile)

assemblies = [pjoin(root_folder, "assemblies", ass, "assembly.fna") for ass in config_file['binnings'][binning_name]['assemblies']]
for ass in assemblies:
    call("cat {ass} >> {temp}/assembly.fna".format(ass = ass, temp = temp_folder), shell = True)

title2log("indexing assembly to temp_folder", logfile)

call("bwa-mem2 index {temp}/assembly.fna >> {out_folder}/logs/binning.log  2>&1".format(temp = temp_folder, threads = threads, out_folder = out_folder), shell=True)

freetxt_line("Starting mappings", logfile)

for lib in config_file['binnings'][binning_name]['libraries']:
    title2log("copying {lib} to temp_folder".format(lib = lib), logfile)
    call("""
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_fwd.fastq.gz >> {temp}/fwd.fastq 2>> {out_folder}/logs/binning.log
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_rev.fastq.gz >> {temp}/rev.fastq 2>> {out_folder}/logs/binning.log
    unpigz -kc {root_folder}/libraries/{lname}/{lname}_unp.fastq.gz >> {temp}/unp.fastq 2>> {out_folder}/logs/binning.log
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder, out_folder = out_folder), shell=True)
    title2log("mapping {lib} to ref".format(lib = lib), logfile)
    call("""
    bwa-mem2 mem -t {threads} {temp}/assembly.fna  -o {temp}/mapping.sam {temp}/fwd.fastq {temp}/rev.fastq 2>> {out_folder}/logs/binning.log
    samtools view -F 3588 -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/{lname}_pairs.bam - >> {out_folder}/logs/binning.log 2>&1
    bwa-mem2 mem -t {threads} {temp}/assembly.fna  -o {temp}/mapping.sam {temp}/unp.fastq 2>> {out_folder}/logs/binning.log
    samtools view -F 3588 -b -S -@{threads}  {temp}/mapping.sam | samtools sort -@ 24 -o {temp}/{lname}_unpaired.bam - >> {out_folder}/logs/binning.log 2>&1
    rm {temp}/mapping.sam 2>> {out_folder}/logs/binning.log
    samtools merge -t {threads} {temp}/{lname}.bam  {temp}/{lname}_pairs.bam  {temp}/{lname}_unpaired.bam >> {out_folder}/logs/binning.log  2>&1
    rm {temp}/{lname}_pairs.bam {temp}/{lname}_unpaired.bam 2>> {out_folder}/logs/binning.log
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder, threads = threads, out_folder = out_folder), shell=True)

call("""
rm {temp}/fwd.fastq
rm {temp}/rev.fastq
rm {temp}/unp.fastq
""".format(temp = temp_folder), shell = True)


title2log("Done with the mappings", logfile)

if method == "metabat" :
    title2log("Making a coverage file for metabat", logfile)
    call("""jgi_summarize_bam_contig_depths --outputDepth {temp}/coverage_table.dat {temp}/*.bam 2>> {out_folder}/logs/binning.log
    """.format(temp=temp_folder, threads = threads, out_folder = out_folder), shell = True)
    title2log("Binning with metabat", logfile)
    call("metabat2 -i {temp}/assembly.fna -o {temp}/bins/bin --unbinned -a {temp}/coverage_table.dat -s {min_bin_size} -t {threads}   >> {out_folder}/logs/binning.log 2>&1".format(temp=temp_folder, threads = threads, out_folder = out_folder, min_bin_size = min_bin_size), shell=True)
    title2log("Cleaning up metabat bins", logfile)
    tbinfoder = "{temp}/bins".format(temp = temp_folder)
    cbinfoder = "{temp}/clean_bins".format(temp = temp_folder)
    os.makedirs(cbinfoder, exist_ok = True)
    bins = [b for b in os.listdir(tbinfoder) if b.endswith(".fa")]
    zeros = len(str(len(bins)))
    for b in tqdm(bins):
        clean_bin(pjoin(tbinfoder, b), pjoin(cbinfoder, binning_name + "_bin-" + b.split(".")[-2].zfill(zeros) + ".fna"), binning_name + "_bin-" + b.split(".")[-2].zfill(zeros))
elif method == "vamb" :
    title2log("Making a coverage file for vamb", logfile)
    call("jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {temp}/coverage_table.dat {temp}/*.bam 2>> {out_folder}/logs/binning.log".format(temp=temp_folder, threads = threads, out_folder = out_folder), shell = True)
    title2log("Binning with vamb", logfile)
    call('vamb --outdir {temp}/bins/ --fasta {temp}/assembly.fna  --jgi {temp}/coverage_table.dat --minfasta {min_bin_size} {cuda} -o "_"'.format(temp=temp_folder, threads = threads, out_folder = out_folder, cuda = "--cuda" if cuda else "", min_bin_size = min_bin_size), shell= True)
    title2log("Cleaning up vamb bins", logfile)
    tbinfoder = "{temp}/bins/bins/".format(temp = temp_folder)
    cbinfoder = "{temp}/clean_bins".format(temp = temp_folder)
    os.makedirs(cbinfoder, exist_ok = True)
    bins = [b for b in os.listdir(tbinfoder) if b.endswith(".fna")]
    zeros = len(str(len(bins)))
    unbinned = [clean_bin(pjoin(tbinfoder, b), pjoin(cbinfoder, binning_name + "_bin-" + b.split("_")[-1].zfill(zeros)), binning_name + "_bin-" + b.split("_")[-1][:-4].zfill(zeros)) for b in tqdm(bins)]
    buffer = []
    zeros = len(str(sum([v[1] for v in unbinned if v])))
    i = 0
    with open(pjoin(cbinfoder, binning_name + "_unbinned.fna"), "w") as handle:
        for file in unbinned:
            if file:
                for s in SeqIO.parse(file[0], "fasta"):
                    s.id = binning_name + "_unbinned_ctg-" + str(i+1).zfill(zeros)
                    i += 1
                    s.description = ""
                    buffer += [s]
                    if len(buffer) > max_buffer_size:
                        SeqIO.write(buffer, handle, "fasta")
                        buffer = []
        SeqIO.write(buffer, handle, "fasta")
    nb_bins = len(os.listdir(cbinfoder)) -1
    zeros = len(str(nb_bins))
    for i,file in enumerate([f for f in os.listdir(cbinfoder) if "unbinned" not in f]):
        new_file = pjoin(cbinfoder,file.split('_bin-')[0] + "_bin-" + str(i+1).zfill(zeros) + ".fna")
        shutil.move(pjoin(cbinfoder,file), new_file)
        call("sed -i 's/{oname}/{nname}/' {file}".format(oname=file[:-4].split('_bin-')[1], nname = str(i+1).zfill(zeros), file = new_file), shell=True)
else :
    print("the binner you want is not implemented yet")
    sys.exit(0)



title2log("Cleaning up and moving the bins", logfile)

if os.path.exists(pjoin(out_folder, "binned_assembly.fna")):
    os.remove(pjoin(out_folder, "binned_assembly.fna"))

os.makedirs(pjoin(out_folder, "bins") , exist_ok = True)
for file in os.listdir(cbinfoder):
    shutil.move(pjoin(cbinfoder,file), pjoin(out_folder, "bins"))
    call("cat {file} >> {ass}".format(file = pjoin(out_folder, "bins", file), ass = pjoin(out_folder, "binned_assembly.fna")), shell=True)

title2log("Binning done", logfile)

shutil.rmtree(temp_folder)
