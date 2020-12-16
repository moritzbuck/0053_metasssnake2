from workflow.scripts.hard_config import assembly_fields, loglinelen, libraries_fields
from datetime import datetime
from math import floor, ceil
from sys import stderr
import subprocess

def is_type(s, typ):
    try:
        typ(s)
        return True
    except ValueError:
        return False

def validate_field(value, validator, name):

    if not validator:
        return value
    if not value:
        return validator['default']
    if 'possibles' in validator:
        assert value in validator['possibles'], 'The value of the field "' + name + '" is '+ str(value) + ' should be in {}'.format(validator['possibles'])
    if 'type' in validator:
        assert is_type(value,validator['type']), 'The value of the field "' + name + '" should be convertible to {}'.format(validator['type'].__name__)
        value = validator['type'](value)
        if 'min' in validator:
            assert value > validator['min'], 'The value of the field "' + name + '" should be larger then ' + str(validator['min'])
        if 'max' in validator:
            assert value < validator['max'], 'The value of the field "' + name + '" should be smaller then ' + str(validator['max'])
    return value

def validate_description_json(file_or_dict):
    import json
    import os
    import sys
    from collections import OrderedDict
    from os.path import join as pjoin

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
        assert os.path.exists(config_dat.get('root_folder', "")), 'No "root_folder", or "root_folder" is an invalid path'
        assert os.path.exists(config_dat.get('temp_folder', "")), 'No "temp_folder", or "temp_folder" is an invalid path'
        assert os.path.exists(pjoin(config_dat['root_folder'],config_dat['config_file'])), 'path to "config" file is wrong it is relative to "root"'
        assert 'libraries' in config_dat, 'no "libraries" in the config-json'
        assert type(config_dat['libraries']) == dict or type(config_dat['libraries']) == OrderedDict, "'libraries' should be a dictionary"
        assert all([k in libraries_fields for k in config_dat['libraries_config'].keys()]), 'All keys of "libraries_config"' + k + '" should be in {}'.format(list(libraries_fields.keys()))
        for f in libraries_fields:
            config_dat['libraries_config'][f] = validate_field(config_dat['libraries'].get(f), libraries_fields[f], f)
        for k, v in config_dat['libraries'].items():
            assert ('fwd' in v and 'rev' in v) or 'unp' in v, 'Library ' + k + ' needs either "fwd" and "rev" or "unp" or all three'
            assert all([ type(v.get(l, [])) == list for l in ['fwd','rev','unp']]) , 'In library ' + k + ': "fwd", "rev", and "unp" need to be lists'
            if 'fwd' in v:
                assert len(v['fwd']) == len(v['rev']), 'In library ' + k + ': "fwd" and "rev" must be of same length'
            for f in sum([ v.get(l, []) for l in ['fwd','rev','unp']], []):
                assert os.path.exists(os.path.join(config_dat['root_folder'], f)), "File " + os.path.join(config_dat['root_folder'], f) + " in Library " + k + " does not exist, all paths are relative to 'root_folder'"
        assert 'assemblies' in config_dat, 'no "assemblies" in the config-json'
        assert type(config_dat['assemblies']) == dict  or type(config_dat['libraries']) == OrderedDict, "'assemblies' should be a dictionary"
        for k, v in config_dat['assemblies'].items():
            assert type(v) == dict  or type(config_dat['libraries']) == OrderedDict, 'Assembly "' + k + '" should be a dictionary'
            assert all([k in assembly_fields for k in v.keys()]), 'All keys of Assembly "' + k + '" should be in {}'.format(list(assembly_fields.keys()))
            for f in assembly_fields:
                v[f] = validate_field(v.get(f), assembly_fields[f], f)
            assert all([l in config_dat['libraries'] for l in v['libraries']]), 'All "libraries" of Assembly "' + k + '" have to be in your set of libraries: {}'.format(list(config_dat['libraries']))
            if len(v['bin_mapping']) == 0:
                v['bin_mapping'] = v['libraries']
            assert all([l in config_dat['libraries'] for l in v['bin_mapping']]), 'All "bin_mapping" of Assembly "' + k + '" have to be in your set of libraries: {}'.format(list(config_dat['libraries']))

    except AssertionError as err:
        print("ERROR : config not valid\nERROR {err}\nERROR Check the doc for formating advice".format(err = err), file = sys.stderr)
        return None
    except Exception as err:
        print("ERROR : Something unexpected went wrong, so your probably not valid\nERROR {err}\nERROR Check the doc for formating advice".format(err = err), file = sys.stderr)
        return None

    if os.path.exists(".git"):
        config_dat['git_commit'] = subprocess.check_output(['git', 'log']).decode().split()[1]

    return config_dat

def title2log(title, logfile, llen = loglinelen, also_stderr = True) :
    text_insert = "{title} started at : {time}".format(title = title, time = datetime.now())
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    if also_stderr:
        print(text, file = stderr, flush = True)

def freetxt_line(text, logfile, llen = loglinelen, also_stderr = True) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")

    if also_stderr:
        print(text, file = stderr, flush = True)



def main():
    import sys

    cline = sys.argv
    if cline[1] == "validate_descriptor":
        test = validate_description_json(cline[2])
        if test:
            print("File " + cline[2] + " is valid")

if __name__ == "__main__":
    main()
