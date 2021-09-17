from datetime import datetime
from math import floor, ceil
from sys import stderr, stdout
import subprocess
import json
import os
import sys

sys.path.append(os.getcwd())
from os.path import join as pjoin
from workflow.scripts.hard_config import *
from workflow._version import __version__

def trystr2int(i):
    try:
        return float(i)
    except ValueError:
        return i

def csv2dict(file, sep = ","):
    with open(file) as handle:
        lines = [l.strip() for l in handle.readlines() if not l.startswith("#")]
        header = lines[0].split(sep)[1:]
        data = {l.split(sep)[0] : l.split(sep)[1:]  for l in lines[1:]}
        data = {k : {kk : trystr2int(vv) for kk,vv in zip(header,v)} for k,v in data.items()}
    return data



def dict2file(dd, file):
    cols = list(list(dd.values())[0].keys())
    with open(file, "w") as handle:
        handle.writelines([",".join([""] + cols) + "\n"] )
        handle.writelines([ k + "," + ",".join([str(v[kk]) for kk in cols]) + "\n"  for k,v in dd.items()])

def is_type(s, typ):
    try:
        typ(s)
        return True
    except ValueError:
        return False

def validate_field(value, validator, name):
    if not validator:
        return value
    if not value or value == '':
        return validator['default']
    if 'possibles' in validator:
        assert value in validator['possibles'] or all([vv in validator['possibles'] for vv in value.split(";")]), 'The value(s) of the field "' + name + '" is(are) '+ str(value) + ', and should be in {}'.format(validator['possibles'])
    if 'type' in validator:
        assert is_type(value,validator['type']), 'The value of the field "' + name + '" should be convertible to {}'.format(validator['type'].__name__)
        if   validator['type'] == bool:
            assert value in ["True", "False"],  'The value of the field "' + name + '" should be either "True" or "False"'
            value = value == "True"
        else :
            value = validator['type'](value)
        if 'min' in validator:
            assert value > validator['min'], 'The value of the field "' + name + '" should be larger then ' + str(validator['min'])
        if 'max' in validator:
            assert value < validator['max'], 'The value of the field "' + name + '" should be smaller then ' + str(validator['max'])
    if name == 'other_parameters':
        pairs = [ tuple(v.split("=")) for v in value.split(";")]
        assert all([len(l) == 2 for l in pairs]), "the other_parameters field is not properly formated"
        value = dict(pairs)
    return value


def generate_config(file_or_dict):
    if type(file_or_dict) != dict:
        try :
            with open(file_or_dict) as handle:
                config_dat = json.load(handle)
        except Exception as err:
            print("ERROR : {file} not valid\nERROR Your description/config json-file is not a valid json file... Check the doc for formating advice".format(file = file_or_dict, err = err), file = sys.stderr)
            return None
    else :
        config_dat = file_or_dict
    try :
        if config_dat['temp_folder'].startswith("$"):
            assert config_dat['temp_folder'][1:] in os.environ, "The env-variable " + config_dat['temp_folder'] + " you picked for temp folder does not exist"
            config_dat['temp_folder'] = os.environ[config_dat['temp_folder'][1:]]
        for k in general_fields:
            config_dat[k] = validate_field(config_dat.get(k), general_fields[k], k)
        for param in necessary_paths:
            assert param in config_dat, param + " should be in your config_file"
            assert os.path.exists(config_dat[param]), "The " + param + " file in your config does not exist, you have " + config_dat[param]

        libraries_dat = csv2dict(config_dat['libraries_file'])
        for k,v in libraries_dat.items():
            v.update({kk : validate_field(libraries_dat[k].get(kk), libraries_fields[kk], kk) for kk in libraries_fields})
            v['fwd'] = v['fwd'].split(";")
            v['rev'] = v['rev'].split(";")
            v['unp'] = v['unp'].split(";")
            for l in v['fwd'] + v['rev'] + v['unp'] :
                assert os.path.exists(pjoin(config_dat['raw_folder'], l)), pjoin(config_dat['raw_folder'], l) + " does not exist"
            assert len(v['fwd']) == len(v['rev']), "The library " + k + " does not have the same number of fwd and reverse libraries"
            assert len(v['fwd']) > 0 or len(v['unp']) > 0, "A library should have at least some paired or an unpaired fastq"
            libraries_dat[k] = v
        config_dat['libraries'] = libraries_dat

        assemblies_dat = csv2dict(config_dat['assemblies_file'])
        for k,v in assemblies_dat.items():
            for field in necessary_assemblies_fields:
                assert field in v, field + " needs to be in the assemblies file"
            v.update({kk : validate_field(assemblies_dat[k].get(kk), assemblies_fields[kk], kk) for kk in assemblies_fields})
            v['libraries'] = v['libraries'].split(";")
            for lib in v['libraries']:
                assert lib in config_dat['libraries'], "no such library as " + lib + " for your assembly " + k
            assemblies_dat[k] = v

        config_dat['assemblies'] = assemblies_dat

        binnings_dat = csv2dict(config_dat['binnings_file'])
        for k,v in binnings_dat.items():
            for field in necessary_binnings_fields:
                assert field in v, field + " needs to be in the binnings file"
            v.update({kk : validate_field(binnings_dat[k].get(kk), binnings_fields[kk], kk) for kk in binnings_fields})
            v['assemblies'] = v['assemblies'].split(";")
            v['libraries'] = v['libraries'].split(";")

            for lib in v['libraries']:
                assert lib in config_dat['libraries'], "no such library as " + lib + " for your binning " + k
            for ass in v['assemblies']:
                assert ass in config_dat['assemblies'], "no such assembly as " + ass + " for your binning " + k

            binnings_dat[k] = v

        config_dat['binnings'] = binnings_dat

        binsets_dat = csv2dict(config_dat['binsets_file'])
        for k,v in binsets_dat.items():
            for field in necessary_binsets_fields:
                assert field in v, field + " needs to be in the binsets file"
            v.update({kk : validate_field(binsets_dat[k].get(kk), binsets_fields[kk], kk) for kk in binsets_fields})
            v['binnings'] = v['binnings'].split(";")
            for bini in v['binnings']:
                assert bini in config_dat['binnings'],  "no such binning as " + bini + " for your binset " + k

            binsets_dat[k] = v
        config_dat['binsets'] = binsets_dat

        mappings_dat = csv2dict(config_dat['mappings_file'])
        for k,v in mappings_dat.items():
            for field in necessary_mappings_fields:
                assert field in v, field + " needs to be in the mappings file"
            v.update({kk : validate_field(mappings_dat[k].get(kk), mappings_fields[kk], kk) for kk in mappings_fields})
            v['libraries'] = v['libraries'].split(";")
            for lib in v['libraries']:
                assert lib in config_dat['libraries'],  "no such library as " + lib + " for your mapping " + k
            assert v['binset'] in config_dat['binsets'],  "no such binset as " + v['binset'] + " for your mapping " + k

            mappings_dat[k] = v
        config_dat['mappings'] = mappings_dat
    except AssertionError as err:
        print("ERROR : config not valid\nERROR {err}\nERROR Check the doc for formating advice".format(err = err), file = sys.stderr)
        return None
    except Exception as err:
        print("ERROR : Something unexpected went wrong, so your probably not valid\nERROR {err}\nERROR Check the doc for formating advice".format(err = err), file = sys.stderr)
        return None

    config_dat['version'] = __version__
    if os.path.exists(".git"):
        config_dat['git_commit'] = subprocess.check_output(['git', 'log']).decode().split()[1]

    return config_dat


def title2log(title, logfile, llen = 90, also_stderr = True) :
    text_insert = "{title} started at : {time}".format(title = title, time = datetime.now())
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    if also_stderr:
        print(text, file = stderr, flush = True)

def freetxt_line(text, logfile, llen = 90, also_stderr = True) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")

    if also_stderr:
        print(text, file = stderr, flush = True)

def folder2csvs(folder, oprefix):
    fastqs = []
    for v in os.walk(folder):
        for vv in v[2]:
            if vv.endswith(".fastq.gz"):
                fastqs += [pjoin(folder,v[0], vv)]
    libraries = {os.path.basename(os.path.dirname(f)) for f in fastqs}
    libraries_dat = dict()
    assemblies_dat = dict()
    binnings_dat = dict()
    for lname in libraries:
        fwds = [l for l in fastqs if (lname + "/") in l  and "_R1_" in l]
        revs = [l for l in fastqs if (lname + "/") in l  and "_R1_" in l]
        libraries_dat[lname] = { 'fwd' : ";".join(fwds), 'rev' : ";".join(revs) }
        assemblies_dat[lname.replace("Sample_", "")] = { 'libraries' : lname }
        binnings_dat[lname.replace("Sample_", "binning-")] = { 'assemblies' : lname.replace("Sample_", "") , 'libraries' : lname}

    binsets_dat = { "all-single-samples" : { 'binnings' : ";".join(list(binnings_dat.keys()))}}


    oass = oprefix + "_assemblies.csv"
    olib = oprefix + "_libraries.csv"
    obin = oprefix + "_binnings.csv"
    oset = oprefix + "_binsets.csv"
    ojson = oprefix + "_config.json"

    cfg = """
    {{
        "root_folder"       : "INSERT_DATA_OUTFOLDER",
        "temp_folder"       : "INSERT_TMP_FOLDER_CAN_BE_ENVVARIABLE_IF_START_WITH_DOLLAR",
        "raw_folder"        : "/",
        "config_file"       : "{ojson}",
        "libraries_file"    : "{olib}",
        "assemblies_file"   : "{oass}",
        "binnings_file"     : "{obin}",
        "binsets_file"      : "{oset}"
    }}
    """.format(ojson = ojson, olib = olib, oass = oass, obin = obin, oset = oset)

    with open(ojson, "w") as handle:
        handle.writelines(cfg)

    dict2file(assemblies_dat, oass)
    dict2file(libraries_dat, olib)
    dict2file(binnings_dat, obin)
    dict2file(binsets_dat, oset)


    return [oass, olib,obin, oset, ojson]

def folder2csvs2(folder, oprefix):
    fastqs = []
    for v in os.walk(folder):
        for vv in v[2]:
            if vv.endswith(".fastq.gz"):
                fastqs += [pjoin(folder,v[0], vv)]
    libraries = {os.path.basename(f).split("_")[1] for f in fastqs}
    libraries_dat = dict()
    assemblies_dat = dict()
    binnings_dat = dict()
    for lname in libraries:
        fwds = [l for l in fastqs if "_" + lname + "_" in os.path.basename(l)  and "_R1_" in l]
        revs = [l for l in fastqs if "_" + lname + "_" in os.path.basename(l)  and "_R2_" in l]
        libraries_dat[lname] = { 'fwd' : ";".join(fwds), 'rev' : ";".join(revs) }
        assemblies_dat[lname.replace("Sample_", "")] = { 'libraries' : lname }
        binnings_dat[lname.replace("Sample_", "binning-")] = { 'assemblies' : lname.replace("Sample_", "") , 'libraries' : lname}

    binsets_dat = { "all-single-samples" : { 'binnings' : ";".join(list(binnings_dat.keys()))}}


    oass = oprefix + "_assemblies.csv"
    olib = oprefix + "_libraries.csv"
    obin = oprefix + "_binnings.csv"
    oset = oprefix + "_binsets.csv"
    ojson = oprefix + "_config.json"

    cfg = """
    {{
        "root_folder"       : "INSERT_DATA_OUTFOLDER",
        "temp_folder"       : "INSERT_TMP_FOLDER_CAN_BE_ENVVARIABLE_IF_START_WITH_DOLLAR",
        "raw_folder"        : "/",
        "config_file"       : "{ojson}",
        "libraries_file"    : "{olib}",
        "assemblies_file"   : "{oass}",
        "binnings_file"     : "{obin}",
        "binsets_file"      : "{oset}"
    }}
    """.format(ojson = ojson, olib = olib, oass = oass, obin = obin, oset = oset)

    with open(ojson, "w") as handle:
        handle.writelines(cfg)

    dict2file(assemblies_dat, oass)
    dict2file(libraries_dat, olib)
    dict2file(binnings_dat, obin)
    dict2file(binsets_dat, oset)


    return [oass, olib,obin, oset, ojson]

def gff2anvio(file):
    import gffutils
    from collections import Counter

    #inspired by https://github.com/karkman/gff_parser

    #Input file and output files
    GFF = file
    OUT_CDS = ""
    OUT_ANNO = ""
    SEP = ':'

    #load in the GFF3 file
    db = gffutils.create_db(GFF, ':memory:')

    #Print headers for anvi'o
    OUT_CDS += "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\n"
    OUT_ANNO += "gene_callers_id\tsource\taccession\tfunction\te_value\n"

    #running gene ID and a trumped-up e-value for the gene calls.
    gene_id = 1
    e_value = "0"

    # keping track of things we haven't processed
    feature_types = Counter()
    call_types = Counter()
    total_num_features = 0
    features_missing_product_or_note = 0

    #parse the GFF3 file and write results to output files
    for feature in db.all_features():
        total_num_features += 1
        # determine source
        source, version = feature.source.split(SEP, 1)

        start = feature.start - 1
        stop = feature.stop

        feature_types[feature.featuretype] += 1
        if feature.featuretype == 'CDS':
            call_type = 1
            call_types['CDS'] += 1
        elif 'RNA' in feature.featuretype:
            call_type = 2
            call_types['RNA'] += 1
        else:
            call_type = 3
            call_types['unknown'] += 1

        if (float(start - stop)/float(3)).is_integer() == True:
            partial = str(0)
        else:
            partial = str(1)

        try:
            gene_acc = feature.attributes['gene'][0]
        except KeyError:
            gene_acc = ""

        # if a feature is missing both, move on.
        if 'product' not in feature.attributes.keys() and 'note' not in feature.attributes.keys():
            features_missing_product_or_note += 1
            continue

        try:
            product = feature.attributes['product'][0]
        except KeyError:
            product = feature.attributes['note'][0]

        # skip if hypotethical proiten:
        if product == 'hypothetical protein':
            product = ""
            gene_acc = ""

        # determine direction
        if feature.featuretype=='repeat_region':
            direction='f'
        else:
            if feature.strand=='+':
                direction='f'
            else:
                direction='r'

        OUT_CDS += '%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\n' % (gene_id, feature.seqid, start, stop, direction, partial, call_type, source, version)
        OUT_ANNO += '%d\t%s:%s\t%s\t%s\t%s\n' % (gene_id, "prokka", source, gene_acc, product, e_value)

        gene_id = gene_id + 1
    return {'cdss' : OUT_CDS, 'annot' : OUT_ANNO}


def main():
    import sys
    cline = sys.argv

    if cline[1] == "validate_descriptor":
        test = generate_config(cline[2])
        if test:
            print("File " + cline[2] + " is valid")
            #json.dump(test[cline[3]], stdout, sort_keys = True, indent = 2)
    if cline[1] == "csv_generator":
        test = folder2csvs(cline[2], cline[3])
        if test:
            print("Files " + ", ".join(test) + " generated")
    if cline[1] == "csv_generator2":
        test = folder2csvs2(cline[2], cline[3])
        if test:
            print("Files " + ", ".join(test) + " generated")
if __name__ == "__main__":
    main()
