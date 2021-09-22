from os.path import join as pjoin
import json
import sys
from workflow.scripts.utils import generate_config

config = generate_config(config)

if not config:
    print("You need to provied a config-file!")
    sys.exit(0)

include : "workflow/rules/rules.smk"


rule all :
    input : expand(pjoin(config['root_folder'], "binsets/{binni_name}/{binni_name}.fna"), binni_name = config['binsets'].keys())

rule all_libs :
    input : expand(pjoin(config['root_folder'], "libraries/{lib_name}/{lib_name}_fwd.fastq.gz"), lib_name = config['libraries'].keys())

rule all_asses :
    input : expand(pjoin(config['root_folder'], "assemblies/{ass_name}/assembly.fna"), ass_name = config['assemblies'].keys())

rule all_binnings :
    input : expand(pjoin(config['root_folder'], "binnings/{binni_name}/binned_assembly.fna"), binni_name = config['binnings'].keys())

rule all_binsets :
    input : expand(pjoin(config['root_folder'], "binsets/{binni_name}/{binni_name}.fna"), binni_name = config['binsets'].keys())
