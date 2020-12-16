from os.path import join as pjoin
import json
import sys
from workflow.scripts.utils import validate_description_json

config = validate_description_json(config)

if not config:
    print("You need to provied a config-file!")
    sys.exit(0)

include : "workflow/rules/rules.smk"

rule all :
    input : expand(config['root_folder'] + "/libraries/{lib_name}/{lib_name}_fwd.fastq.gz", lib_name = config['libraries'].keys())
