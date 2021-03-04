import pandas

"ls | parallel -j 20 sourmash compute --track-abundance --merge {}  -k 31 --scaled 1000 {}/{}_fwd.fastq.gz {}/{}_rev.fastq.gz {}/{}_unp.fastq.gz -o {}/{}.sig -p 1"
run("sourmash compare --ignore-abundance -p 20 --csv ../compared_noabs.tsv --traverse-directory  .")

with open("../compared_noabs.tsv") as handle:
    hh = [handle.readline().strip().split(",")][0]
    chosenones = {k : sorted([(kk,float(ll)) for kk,ll in zip(hh,l.split(","))], key = lambda t : t[1], reverse = True)[0:10] for k,l in zip(hh,handle.readlines())}
    chosenones = {k.replace("Sample_","") : ";".join([vv[0] for vv in v]) for k,v in chosenones.items()}
    chosenones['assemblies'] = "libraries"

with open("/home/moritz/uppmax/projects/0056_anderspolar/anderspolar_binnings_fixed.csv") as handle:
            lines = [l[:-1] + (",\n" if l.split(",")[1] not in chosenones else ("," + chosenones[l.split(",")[1]] + "\n") )    for l in handle]
with open("/home/moritz/uppmax/projects/0056_anderspolar/anderspolar_binnings_fixed.csv","w") as handle:
    handle.writelines(lines)
